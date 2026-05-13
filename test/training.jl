using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DiscreteCurves          
using Zygote           
using Optim
using LinearAlgebra
using Random                 
using Statistics 
using Plots
using Printf
Random.seed!(42)


struct RobustCurveParams{T <: Real}
    p23 :: T  
    q11 :: T 
    q2  :: T
    q33 :: T
end

RobustCurvParams(v::AbstractVector{T}) where T =
    RobustCurvParams{T}(v[1], v[2], v[3], v[4])
to_vector(p::RobustCurveParams) = [p.p23, p.q11, p.q2, p.q33]

pinkall_params() = RobustCurveParams(0.0, 0.0, 0.0, 0.0)


function compute_features(γ::AbstractMatrix{T}) where T
    ℓ    = edge_lengths(γ)
    φ    = turning_angles(γ)
    lstar = dual_lengths(γ)
    N    = normals(γ)

    n = size(γ, 1)

    
    ℓ_prev = circshift(ℓ, 1)
    ε_clamp = T(1e-12)
    x1 = log.(max.(ℓ_prev, ε_clamp) ./ max.(ℓ, ε_clamp))
    x2 = cos.(φ ./ 2)
    x3 = sin.(φ ./ 2)
    return x1, x2, x3, lstar, N
end

function robust_curvature(γ::AbstractMatrix{T}, p::RobustCurveParams) where T
    x1, x2, x3, lstar, _ = compute_features(γ)

    
    num = @. 2x3 + p.p23 * x2 * x3

    δ = T(1e-6)
    denom = @. max(1 + p.q11 * x1^2 + p.q2 * x2 + p.q33 * x3^2, δ)

    κ = @. num / (lstar * denom)

    return κ
end
function pinkall_curvature(γ::AbstractMatrix{T}) where T
    ℓ    = edge_lengths(γ)
    φ    = turning_angles(γ)
    lstar = dual_lengths(γ)
    return @. 2sin(φ / 2) / lstar
end


function flow_step(γ::AbstractMatrix{T}, κ::AbstractVector{T}, Δt::T) where T
    N = normals(γ)

    κ_centered = κ .- mean(κ)

    
    γ_new = γ .+ Δt .* (κ_centered .* N)

    return γ_new
end

function isoperimetric_ratio(γ::AbstractMatrix{T}) where T
    L = arc_length(γ)    
    A = abs(signed_area(γ))
    ε = T(1e-10)
    return L^2 / (4T(π) * max(A, ε))
end

function snapshot_loss(γ::AbstractMatrix{T}, p::RobustCurveParams, Δt::T) where T
    κ     = robust_curvature(γ, p)
    γ_new = flow_step(γ, κ, Δt)
    return isoperimetric_ratio(γ_new)
end

function make_circle(n::Int; R::Float64=1.0)
    θ = range(0, 2π, length=n+1)[1:end-1]
    hcat(R .* cos.(θ), R .* sin.(θ))
end

function make_ellipse(n::Int; a::Float64=1.5, b::Float64=0.7)
    θ = range(0, 2π, length=n+1)[1:end-1]
    hcat(a .* cos.(θ), b .* sin.(θ))
end

function make_noisy_circle(n::Int; R::Float64=1.0, σ::Float64=0.15)
    γ_clean = make_circle(n; R)
    noise   = σ .* randn(n, 2)
    return γ_clean .+ noise
end

function make_epsilon_edge_curve(n::Int; R::Float64=1.0, ε_frac::Float64=0.05,
                                  pos::Union{Int,Nothing}=nothing)
    γ = make_circle(n; R)
    k = isnothing(pos) ? rand(1:n) : pos

    k_next = mod1(k + 1, n)

    edge_vec = γ[k_next, :] .- γ[k, :]
    new_pos  = γ[k_next, :] .- ε_frac .* edge_vec

    γ_new = copy(γ)
    γ_new[k, :] = new_pos
    return γ_new
end

function generate_dataset(n_verts::Int; n_samples::Int=30)
    dataset = Vector{Matrix{Float64}}()

    n_smooth  = round(Int, 0.30 * n_samples)
    n_noisy   = round(Int, 0.40 * n_samples)
    n_epsilon = n_samples - n_smooth - n_noisy

    # ── Smooth curves ───────────────────────────────────────────────────────
    for i in 1:n_smooth
        if isodd(i)
            push!(dataset, make_circle(n_verts; R=0.5 + 0.5*rand()))
        else
            push!(dataset, make_ellipse(n_verts; a=1.0+0.5*rand(), b=0.5+0.3*rand()))
        end
    end

    # ── Noisy circles ───────────────────────────────────────────────────────
    for _ in 1:n_noisy
        σ = 0.05 + 0.15*rand()   # noise level between 5% and 20%
        push!(dataset, make_noisy_circle(n_verts; σ))
    end

    # ── ε-edge curves ───────────────────────────────────────────────────────
    for _ in 1:n_epsilon
        ε = 0.02 + 0.08*rand()   # ε-fraction between 2% and 10%
        push!(dataset, make_epsilon_edge_curve(n_verts; ε_frac=ε))
    end

    shuffle!(dataset)
    return dataset
end


function train_robust_operator(
    dataset :: Vector{Matrix{Float64}};
    Δt      :: Float64 = 0.01,
    n_epochs:: Int     = 50,
    verbose :: Bool    = true,
)
    
    θ_init = [0.01, 0.1, 0.0, 0.05]

    softplus(x) = log1p(exp(x))

    function θ_to_params(θ)
        RobustCurveParams(
            θ[1],             # p23  : unconstrained
            softplus(θ[2]),   # q11  : must be ≥ 0
            θ[3],             # q2   : unconstrained
            softplus(θ[4]),   # q33  : must be ≥ 0
        )
    end

    
    function objective(θ::AbstractVector)
        p    = θ_to_params(θ)
        losses = map(γ -> snapshot_loss(γ, p, Float64(Δt)), dataset)
        return mean(losses)
    end

    function objective_and_grad!(G, θ)
        val, back = Zygote.withgradient(objective, θ)
        ∇θ = back[1]
        if !isnothing(G) && !isnothing(∇θ)
            G .= ∇θ
        end
        return val
    end

    loss_history = Float64[]
    epoch_hook   = Optim.Options(
        iterations       = n_epochs,
        show_trace       = verbose,
        show_every       = 5,
        callback         = state -> begin
            push!(loss_history, state.value)
            false
        end,
    )
    result = Optim.optimize(
        Optim.only_fg!(objective_and_grad!),
        θ_init,
        LBFGS(),
        epoch_hook,
    )

    θ_opt = Optim.minimizer(result)
    p_opt = θ_to_params(θ_opt)

    if verbose
        println("\n━━━ Training complete ━━━")
        println("  Final mean IQ loss : $(round(Optim.minimum(result), digits=6))")
        println("  Trained parameters :")
        println("    p23 = $(round(p_opt.p23, digits=5))")
        println("    q11 = $(round(p_opt.q11, digits=5))   ← ε-edge suppressor")
        println("    q2  = $(round(p_opt.q2,  digits=5))")
        println("    q33 = $(round(p_opt.q33, digits=5))   ← angle self-saturation")
    end

    return p_opt, loss_history
end



function simulate_flow(
    γ₀          :: Matrix{Float64},
    curvature_fn :: Function,
    n_steps     :: Int,
    Δt          :: Float64,
)
    trajectory = Vector{Matrix{Float64}}(undef, n_steps + 1)
    trajectory[1] = γ₀

    # Reference diameter for blow-up detection
    diam = maximum(norm(γ₀[i,:] - γ₀[j,:]) for i in 1:size(γ₀,1) for j in i+1:size(γ₀,1))

    γ = copy(γ₀)
    for t in 1:n_steps
        κ = curvature_fn(γ)
        γ = flow_step(γ, κ, Float64(Δt))

        # ── Blow-up guard ──────────────────────────────────────────────────
        if any(!isfinite, γ) || maximum(abs, γ) > 10 * diam
            @warn "Flow diverged at step $t — truncating trajectory"
            return trajectory[1:t]
        end

        trajectory[t + 1] = γ
    end

    return trajectory
end


# ============================================================
# 9.  VISUALISATION
# ============================================================

"""
    plot_flow_comparison(γ₀, p_trained; n_steps, Δt, n_show)

Generate a side-by-side (or grid) comparison:
  Left  : flow with classical Pinkall curvature  (unstable on ε-edge input)
  Right : flow with robust trained operator       (stable)

Shows the initial curve (gray), intermediate snapshots (colored by time),
and the final curve (red).
"""
function plot_flow_comparison(
    γ₀        :: Matrix{Float64},
    p_trained :: RobustCurveParams;
    n_steps   :: Int    = 200,
    Δt        :: Float64 = 0.005,
    n_show    :: Int    = 5,       # number of intermediate curves to draw
)
    println("Simulating Pinkall flow...")
    traj_pinkall = simulate_flow(γ₀, pinkall_curvature, n_steps, Δt)

    println("Simulating Robust flow...")
    traj_robust  = simulate_flow(γ₀, γ -> robust_curvature(γ, p_trained), n_steps, Δt)

    # Compute IQ evolution for both
    iq_pinkall = [isoperimetric_ratio(γ) for γ in traj_pinkall]
    iq_robust  = [isoperimetric_ratio(γ) for γ in traj_robust]

    # ── Colour scheme: early = blue, late = red ────────────────────────────
    colors = [RGB(t, 0.2, 1-t) for t in range(0, 1, length=n_show)]

    # ── Helper: plot a trajectory ──────────────────────────────────────────
    function draw_traj!(plt, traj, title_str)
        n_frames = length(traj)
        snap_idx = round.(Int, range(1, n_frames, length=n_show))

        # Initial curve
        γ_init = traj[1]
        plot!(plt, [γ_init[:,1]; γ_init[1,1]], [γ_init[:,2]; γ_init[1,2]],
              color=:gray70, lw=1.5, label="t=0", ls=:dash)

        for (j, idx) in enumerate(snap_idx)
            γ = traj[idx]
            label_str = j == n_show ? "t=$(n_frames-1) (final)" : ""
            lw = j == n_show ? 2.5 : 1.0
            plot!(plt, [γ[:,1]; γ[1,1]], [γ[:,2]; γ[1,2]],
                  color=colors[j], lw=lw, label=label_str)
        end

        # Mark ε-edge vertex (initial) with a star
        κ_init  = pinkall_curvature(traj[1])
        k_worst = argmax(abs.(κ_init))
        scatter!(plt, [traj[1][k_worst,1]], [traj[1][k_worst,2]],
                 marker=:star5, ms=8, color=:orange, label="ε-vertex")

        title!(plt, title_str)
        xlabel!(plt, "x"); ylabel!(plt, "y")
        aspect_ratio!(plt, 1)
    end

    p1 = plot(legend=:outertopright, title="Pinkall (classical)")
    draw_traj!(p1, traj_pinkall, "Classical Pinkall Flow\n(unstable on ε-edge)")

    p2 = plot(legend=:outertopright)
    draw_traj!(p2, traj_robust, "Robust Trained Operator\n(stable)")

    # ── IQ evolution plot ──────────────────────────────────────────────────
    p3 = plot(iq_pinkall, label="Pinkall IQ",  color=:crimson, lw=2,
              xlabel="Step", ylabel="IQ = L²/(4πA)",
              title="Isoperimetric Ratio Evolution", yscale=:log10)
    plot!(p3, iq_robust,  label="Robust IQ",   color=:steelblue, lw=2)
    hline!(p3, [1.0], color=:black, ls=:dot, label="IQ=1 (circle limit)")

    # ── Combine ───────────────────────────────────────────────────────────
    fig = plot(p1, p2, p3, layout=@layout([a b; c]), size=(1000, 800),
               plot_title="Robust vs Classical Discrete Curvature Flow")

    return fig
end

"""
    plot_feature_landscape(p_trained, n_verts)

Visualise how the trained operator modulates κ compared to Pinkall
as a function of the log-length-ratio x₁ (the ε-edge detector).
Shows a 1D slice through feature space with x₂=0, x₃ fixed.
"""
function plot_feature_landscape(p_trained::RobustCurveParams; x3_val::Float64=0.3)
    x1_range = range(-3, 3, length=200)
    x2_val   = 0.0   # neutral stretch

    # Pinkall response (no denominator correction)
    κ_pinkall = @. 2x3_val   # 2sin(φ/2) without the 1/l* (l*=1 normalisation)

    # Robust response along x₁ axis
    function robust_response(x1)
        num   = 2x3_val + p_trained.p23 * x2_val * x3_val
        denom = max(1 + p_trained.q11 * x1^2 + p_trained.q2 * x2_val +
                    p_trained.q33 * x3_val^2, 1e-6)
        return num / denom
    end

    κ_robust = robust_response.(x1_range)

    p = plot(x1_range, fill(κ_pinkall, length(x1_range)),
             label="Pinkall (constant)", color=:crimson, lw=2, ls=:dash,
             xlabel="x₁ = ln(ℓᵢ₋₁/ℓᵢ)  [log-length ratio]",
             ylabel="Curvature response (× l*)",
             title="Curvature as a function of ε-edge indicator x₁\n(x₂=0, x₃=$(x3_val))")
    plot!(p, x1_range, κ_robust, label="Robust (trained)", color=:steelblue, lw=2.5)
    vspan!(p, [-0.5, 0.5], alpha=0.1, color=:green, label="Uniform mesh zone")
    vspan!(p, [2.0, 3.0],  alpha=0.15, color=:red,   label="ε-edge zone")
    vspan!(p, [-3.0, -2.0], alpha=0.15, color=:red,   label="")

    return p
end



"""
    diagnose_curve(γ)

Print a summary of geometric quantities for a discrete curve.
Useful during development to understand the feature distribution.
"""
function diagnose_curve(γ::Matrix{Float64}; label::String="")
    ℓ     = edge_lengths(γ)
    φ     = turning_angles(γ)
    lstar = dual_lengths(γ)
    x1, x2, x3, _, _ = compute_features(γ)

    println("━━━ Curve diagnostics$(isempty(label) ? "" : " [$label]") ━━━")
    println("  n_verts : $(size(γ,1))")
    @printf("  Edge lengths : min=%.4f  max=%.4f  ratio=%.2f\n",
            minimum(ℓ), maximum(ℓ), maximum(ℓ)/minimum(ℓ))
    @printf("  |x₁| (log-ratio) : max=%.3f  (ε-edge if > 2)\n", maximum(abs, x1))
    @printf("  |x₂| (stretch)   : max=%.3f\n", maximum(abs, x2))
    @printf("  x₃ range         : [%.3f, %.3f]\n", minimum(x3), maximum(x3))
    @printf("  IQ               : %.5f  (1.0 = circle)\n", isoperimetric_ratio(γ))
    println()
end


function main()
    println("═══════════════════════════════════════════════════════")
    println("  Robust Discrete Curvature Operator — Training Pipeline")
    println("═══════════════════════════════════════════════════════\n")

    n_verts = 32   
    Δt      = 0.005

    
    println("Generating dataset...")
    dataset = generate_dataset(n_verts; n_samples=40)
    println("  Generated $(length(dataset)) training curves.\n")

    # ── Step 2: Diagnostics on a representative ε-edge curve ───────────────
    γ_demo = make_epsilon_edge_curve(n_verts; ε_frac=0.04, pos=8)
    diagnose_curve(γ_demo; label="ε-edge demo curve")

    # ── Step 3: Train ─────────────────────────────────────────────────────
    println("Training robust curvature operator via Zygote + L-BFGS...\n")
    p_trained, loss_history = train_robust_operator(
        dataset; Δt=Δt, n_epochs=80, verbose=true
    )

    # ── Step 4: Visualise comparison flow ─────────────────────────────────
    println("\nGenerating comparison flow visualisation...")
    # Use an ε-edge curve as the stress-test input
    γ_test = make_epsilon_edge_curve(n_verts; ε_frac=0.04, pos=8)
    fig_flow = plot_flow_comparison(γ_test, p_trained;
                                    n_steps=300, Δt=Δt, n_show=6)
    savefig(fig_flow, "flow_comparison.png")
    println("  Saved: flow_comparison.png")

    # ── Step 5: Feature landscape ─────────────────────────────────────────
    println("Generating feature landscape plot...")
    fig_feat = plot_feature_landscape(p_trained; x3_val=0.3)
    savefig(fig_feat, "feature_landscape.png")
    println("  Saved: feature_landscape.png")

    # ── Step 6: Training loss curve ────────────────────────────────────────
    fig_loss = plot(loss_history, xlabel="Iteration", ylabel="Mean IQ Loss",
                    title="Training Loss (Snapshot Approach)\nMinimise IQ → regular n-gon",
                    color=:steelblue, lw=2, legend=false, yscale=:log10)
    hline!(fig_loss, [1.0], color=:black, ls=:dot, label="IQ=1")
    savefig(fig_loss, "training_loss.png")
    println("  Saved: training_loss.png")

    # ── Summary ────────────────────────────────────────────────────────────
    println("\n═══════════════════════════════════════════════════════")
    println("  Final trained parameters:")
    println("    p23 = $(round(p_trained.p23, digits=6))")
    println("    q11 = $(round(p_trained.q11, digits=6))  ← ε-edge suppressor")
    println("    q2  = $(round(p_trained.q2,  digits=6))")
    println("    q33 = $(round(p_trained.q33, digits=6))  ← angle self-saturation")
    println()
    println("  Interpretation:")
    q11 = p_trained.q11
    println("    q11=$(round(q11,digits=3)): at an ε-edge with ratio 10:1,")
    println("    x₁=ln(10)≈2.3, denominator≈$(round(1+q11*2.3^2, digits=2))×,")
    println("    so curvature is suppressed by factor $(round(1/(1+q11*2.3^2), digits=3)).")
    println("═══════════════════════════════════════════════════════")

    return p_trained, loss_history
end

p_trained, loss_history = main()
