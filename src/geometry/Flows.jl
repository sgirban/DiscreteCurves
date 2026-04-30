module Flows

using LinearAlgebra
using StaticArrays
using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators, ..Lengths, ..Orientation
using ..Curvature: _turning_angle_at, _kn_B, _curvature_at,
	TurningAngle, SteinerCurvature, TangentCurvature, OsculatingCircle,
	SignedCurvature2D, LengthWeightedCurvature, CustomCurvature


export CurvatureFlow, PointwiseFlow, BulkFlow, GradientFlow
export ScaledFlow, SumFlow, ConstrainedFlow, TangentialCorrection
export compute_velocities!
export cfl_timestep
export FlowResult, evolve, evolve_step!
export curve_shortening_flow

"""
	CurvatureFlow(model=SteinerCurvature(); tangential=false)
 
Curvature-driven flow.  Normal velocity: `vᵢ = κᵢ(model) · N̂ᵢ`.
 
**`tangential=false` (default):** pure normal flow. Vertices move only inward/outward.
  - Efficient, no reparametrisation overhead.
  - Edge lengths become unequal during long flows → use `resample_every` to fix.
 
**`tangential=true`:** adds a tangential velocity component that keeps vertices
  equally spaced along arc length.  This is the *DeTurck trick* (a.k.a.
  arc-length-equalising tangential reparametrisation):
  - Adds `vᵢᵀ = (ℓᵢ − ℓᵢ₋₁)/(ℓᵢ + ℓᵢ₋₁) · Tᵢ` where Tᵢ is the unit tangent.
  - Does NOT change the geometric shape of the curve (same orbits, different speed).
  - Makes the scheme better conditioned for long flows (no mesh degeneration).
  - Slightly more expensive: one extra normalisation per vertex.
 
```julia
CurvatureFlow()                               # κ^B, normal only
CurvatureFlow(TurningAngle())                 # κ^A, normal only
CurvatureFlow(SteinerCurvature(); tangential=true)  # κ^B + arc-length reparam
```
"""
struct CurvatureFlow{M}
	model      :: M
	tangential :: Bool
end
CurvatureFlow(model = SteinerCurvature(); tangential::Bool = false) =
	CurvatureFlow{typeof(model)}(model, tangential)

# ── PointwiseFlow ─────────────────────────────────────────────────────────────
 
"""
    PointwiseFlow(f)
 
Flow where `f(curve, i) → SVector{N,T}` gives the velocity at vertex `i`.
`f` receives the full curve, so it may reference any vertex or global quantity.
 
```julia
# Attraction toward centroid
PointwiseFlow((c, i) -> centroid(c) - vertex(c, i))
 
# Repulsion from the origin
PointwiseFlow((c, i) -> vertex(c, i) / norm(vertex(c, i))^2)
```
"""
struct PointwiseFlow{F}
    f :: F 
end

"""
    BulkFlow(f)
 
Flow where `f(vel, curve) → nothing` fills the entire velocity buffer at once.
Preferred when all velocities are computed in a single coordinated pass.
 
```julia
BulkFlow() do vel, c
    κs = curvature_vectors(c)
    for i in eachindex(κs)
        vel[i] = κs[i]
    end
end
```
"""
struct BulkFlow{F}
    f :: F
end

"""
    GradientFlow(f; h=1e-6)
 
Gradient-descent flow minimising scalar functional `f(curve) → T`.
Velocity: `vᵢ = −∇_γᵢ f`  via central finite differences with step `h`.
 
Cost: O(2·N·n) functional evaluations per step. Fine for n ≤ 1000 with cheap f.
For arc-length gradient, use `CurvatureFlow(SteinerCurvature())` instead —
they are theoretically equivalent but CurvatureFlow is exact and O(n).
 
```julia
# Minimise bending energy
GradientFlow(total_squared_curvature)
 
# Minimise Hausdorff distance to a target curve
GradientFlow(c -> hausdorff_distance(c, target))
 
# Minimise enclosed area (collapses the curve)
GradientFlow(enclosed_area; h=1e-5)
```
"""
struct GradientFlow{F, T <: AbstractFloat}
    f :: F
    h :: T
end
GradientFlow(f; h::Real=1e-6) = GradientFlow(f, Float64(h))

"""
    ScaledFlow(law, α)
 
Multiply all velocities by constant `α`.
 
```julia
ScaledFlow(CurvatureFlow(), 0.5)   # half-speed curve shortening
ScaledFlow(CurvatureFlow(), -1.0)  # curve expansion (reverse)
```
"""
struct ScaledFlow{L, S <: Real}
    law   :: L
    scale :: S
end
 
"""
    SumFlow(law1, law2)
 
Additive superposition: `vᵢ = v1ᵢ + v2ᵢ`.
 
```julia
SumFlow(
    CurvatureFlow(),
    ScaledFlow(PointwiseFlow((c,i) -> centroid(c) - vertex(c,i)), 0.1)
)
```
"""
struct SumFlow{L1, L2}
    law1 :: L1
    law2 :: L2
end
 
"""
    ConstrainedFlow(law, g!)
 
Apply law, then project velocities via `g!(vel, c) → nothing`.
Use for hard constraints (fixed endpoints, prescribed total velocity, etc.).
 
```julia
# Fix both endpoints of an open curve
ConstrainedFlow(CurvatureFlow(), (vel, c) -> begin
    vel[1] = zero(eltype(vel))
    vel[end] = zero(eltype(vel))
end)
 
# Zero-mean constraint (remove drift)
ConstrainedFlow(CurvatureFlow(), (vel, c) -> begin
    drift = sum(vel) / length(vel)
    for i in eachindex(vel); vel[i] -= drift; end
end)
```
"""
struct ConstrainedFlow{L, G}
    law        :: L
    constraint :: G 
end

"""
    TangentialCorrection(law)
 
Wrap any law with an arc-length-equalising tangential component.
This is equivalent to `CurvatureFlow(model; tangential=true)` for curvature laws
but works with ANY velocity law.
 
The tangential component:
  `vᵢᵀ = ½(ℓᵢ − ℓᵢ₋₁) · Tᵢ`
moves vertex `i` along the curve toward the midpoint of its two neighbours,
keeping edge lengths equal without changing the geometric shape.
 
```julia
TangentialCorrection(CurvatureFlow(TurningAngle()))
TangentialCorrection(GradientFlow(enclosed_area))
TangentialCorrection(BulkFlow(my_flow_fn))
```
"""
struct TangentialCorrection{L}
    law :: L
end
"""
    compute_velocities!(vel, c, law)
 
Fill the velocity buffer `vel` from curve `c` according to `law`.
 
**This is the single extensibility point.** Define this method for your custom law
and it integrates automatically with `evolve`, `cfl_timestep`, tracking, resampling,
and all compositors (`SumFlow`, `ScaledFlow`, etc.).
 
Contract:
- Fill `vel[i]` for every `i` in `1:nvertices(c)`.
- Read `c` only (do not modify it).
- May not allocate per-vertex (pre-allocate outside if needed).
- Fixed boundary for open curves: set `vel[1] = vel[end] = zero(SVector{N,T})`.
 
```julia
# Example: uniform inward contraction
struct ContractLaw end
function DiscreteCurves.Flow.compute_velocities!(vel, c, ::ContractLaw)
    ctr = centroid(c)
    for i in 1:nvertices(c)
        vel[i] = ctr - vertex(c, i)   # point toward centroid
    end
end
 
r = evolve(my_curve, ContractLaw(); nsteps=100, track=[:L => arc_length])
```
"""
function compute_velocities! end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::CurvatureFlow{SteinerCurvature}) where {N,T}
    n      = nvertices(c)
    closed = isclosed(c)
    @inbounds for i in 1:n
        if !closed && (i == 1 || i == n)
            vel[i] = zero(SVector{N,T}); continue
        end
        ip = _adj(i-1, n, closed);  iq = _adj(i+1, n, closed)
        d1 = c.points[i] - c.points[ip];  n1 = norm(d1)
        d2 = c.points[iq] - c.points[i];  n2 = norm(d2)
        if n1 < eps(T) || n2 < eps(T)
            vel[i] = zero(SVector{N,T}); continue
        end
        v_normal = d2/n2 - d1/n1
        vel[i]   = law.tangential ? v_normal + _tangential_v(d1, d2, n1, n2) : v_normal
    end
end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::CurvatureFlow{M}) where {N,T,M}
    n      = nvertices(c)
    closed = isclosed(c)
    @inbounds for i in 1:n
        if !closed && (i == 1 || i == n)
            vel[i] = zero(SVector{N,T}); continue
        end
        ip = _adj(i-1, n, closed);  iq = _adj(i+1, n, closed)
        p, q, r = c.points[ip], c.points[i], c.points[iq]
        d1 = q - p;  n1 = norm(d1)
        d2 = r - q;  n2 = norm(d2)
        if n1 < eps(T) || n2 < eps(T)
            vel[i] = zero(SVector{N,T}); continue
        end
        kn      = d2/n2 - d1/n1
        kn_norm = norm(kn)
        if kn_norm < eps(T)
            vel[i] = zero(SVector{N,T}); continue
        end
        κ        = _curvature_at(p, q, r, law.model)
        v_normal = (κ / kn_norm) * kn
        vel[i]   = law.tangential ? v_normal + _tangential_v(d1, d2, n1, n2) : v_normal
    end
end

function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::PointwiseFlow) where {N,T}
    @inbounds for i in 1:nvertices(c)
        vel[i] = SVector{N,T}(law.f(c, i))
    end
end

function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::BulkFlow) where {N,T}
    law.f(vel, c)
end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::GradientFlow) where {N,T}
    n        = nvertices(c)
    h        = T(law.h)
    pts_work = copy(c.points)
 
    @inbounds for i in 1:n
        grad = zero(MVector{N,T})
        p_i  = c.points[i]
        for d in 1:N
            pts_work[i] = _svec_setindex(p_i, p_i[d] + h, d)
            f_fwd = law.f(DiscreteCurve(pts_work; closed=isclosed(c)))
            pts_work[i] = _svec_setindex(p_i, p_i[d] - h, d)
            f_bwd = law.f(DiscreteCurve(pts_work; closed=isclosed(c)))
            grad[d] = (f_fwd - f_bwd) / (2h)
            pts_work[i] = p_i
        end
        vel[i] = -SVector{N,T}(grad)
    end
end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::ScaledFlow) where {N,T}
    compute_velocities!(vel, c, law.law)
    α = T(law.scale)
    @inbounds for i in eachindex(vel); vel[i] = α * vel[i]; end
end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::SumFlow) where {N,T}
    vel2 = similar(vel)  
    compute_velocities!(vel,  c, law.law1)
    compute_velocities!(vel2, c, law.law2)
    @inbounds for i in eachindex(vel); vel[i] = vel[i] + vel2[i]; end
end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::ConstrainedFlow) where {N,T}
    compute_velocities!(vel, c, law.law)
    law.constraint(vel, c)
end
function compute_velocities!(vel::Vector{SVector{N,T}},
                              c::AbstractDiscreteCurve{N,T},
                              law::TangentialCorrection) where {N,T}
    compute_velocities!(vel, c, law.law)
    n      = nvertices(c)
    closed = isclosed(c)
    @inbounds for i in 1:n
        ((!closed) && (i == 1 || i == n)) && continue
        ip = _adj(i-1, n, closed);  iq = _adj(i+1, n, closed)
        d1 = c.points[i] - c.points[ip]
        d2 = c.points[iq] - c.points[i]
        n1 = norm(d1);  n2 = norm(d2)
        vel[i] += _tangential_v(d1, d2, n1, n2)
    end
end
#─────────────────────────────────────────────────────────────────────────────

"""
    cfl_timestep(c; safety=0.25)  →  T
 
CFL-stable timestep τ for curvature flows.
 
τ has units of [length²]. A curve of edge length 0.1 needs τ ≤ 0.0025.
If your curve has very short edges and flow seems to explode, reduce safety.
 
```julia
τ = cfl_timestep(c)              # τ = 0.25 × min_edge²
τ = cfl_timestep(c; safety=0.5) # τ = 0.50 × min_edge² (faster, less stable)
```
"""
function cfl_timestep(c::AbstractDiscreteCurve{N,T}; safety::Real=T(0.25)) where {N,T}
    min_ℓ² = typemax(T)
    n, closed = nvertices(c), isclosed(c)
    ne = closed ? n : n-1
    @inbounds for i in 1:ne
        j  = closed ? mod1(i+1, n) : i+1
        ℓ² = sum(abs2, c.points[j] - c.points[i])
        ℓ² < min_ℓ² && (min_ℓ² = ℓ²)
    end
    T(safety) * min_ℓ²
end
#
# ── Flow Result ───────────────────────────────────────────────────────────────
"""
    FlowResult{N,T}
 
Output of `evolve`. Carries the final curve, time series of tracked functionals,
and optional snapshots of intermediate curves.
 
### Field access
```julia
r.curve        # final DiscreteCurve after all steps
r.nsteps       # actual number of steps taken (≤ requested nsteps)
r.time         # total physical time elapsed: ∑ τₖ
r.max_v        # Vector{T}: max vertex speed at each step (auto-tracked)
r.curves       # Vector of saved intermediate curves (empty if save_every=0)
r.times        # Vector{T}: physical times of each saved curve snapshot
 
r[:L]          # Vector of arc lengths at each tracked step
r[:A]          # Vector of areas at each tracked step
r[:κmax]       # Vector of max curvatures, etc.
r[σ]           # DiscreteCurve at step σ  (requires save_every > 0)
r[σ => :L]     # tracked :L at the saved snapshot closest to step σ
```
 
### σ vs τ
- `σ` (sigma) is the step index: `r[5]` gives the curve after 5 Euler steps.
- `τ` (tau) is the physical timestep: `r.time` = ∑ τ gives total flow time.
- Physical time at step σ: `r.times[σ]` (if saved).
- These are independent: 100 small steps and 10 large steps can give the same
  flow time but very different accuracy and curve shape.
"""
struct FlowResult{N, T <: AbstractFloat}
    curve   :: AbstractDiscreteCurve{N,T}
    nsteps  :: Int
    time    :: T
    max_v   :: Vector{T}
    history :: Dict{Symbol, AbstractVector}
    curves  :: Vector{<:AbstractDiscreteCurve{N,T}}
    times   :: Vector{T}
end

Base.getindex(r::FlowResult, k::Symbol) = r.history[k]
Base.haskey(r::FlowResult, k::Symbol)   = haskey(r.history, k)

function Base.getindex(r::FlowResult, σ::Int)
    isempty(r.curves) &&
        error("No curve snapshots saved. Run evolve with save_every > 0 to save curves.\n" *
              "Example: evolve(c, law; save_every=10)")
    idx = σ + 1
    checkbounds(r.curves, idx)
    r.curves[idx]
end
function Base.getindex(r::FlowResult, p::Pair{Int,Symbol})
    σ, sym = p
    isempty(r.curves) && error("save_every=0; no curve snapshots")
    r[sym][σ + 1]
end
 
function Base.show(io::IO, r::FlowResult{N,T}) where {N,T}
    tracked = filter(!=(:max_v), collect(keys(r.history)))
    snap_str = isempty(r.curves) ? "" : ", $(length(r.curves)) snapshots"
    print(io, "FlowResult{$N,$T}  σ=$(r.nsteps) steps  t=$(round(r.time; sigdigits=4))" *
              "  tracked=$(tracked)$snap_str")
end


"""
    evolve_step!(pts, vel, c_view, law, τ) → max_v::T
 
One explicit Euler step: compute velocities from OLD positions, then update.
Returns maximum velocity magnitude (for convergence monitoring).
```julia
pts    = copy(c.points)
vel    = similar(pts)
c_view = DiscreteCurve{N,T}(pts, isclosed(c))
τ      = cfl_timestep(c)
law    = CurvatureFlow()
 
for σ in 1:10_000
    max_v = evolve_step!(pts, vel, c_view, law, τ)
    max_v < 1e-9 && break
end
final = DiscreteCurve(copy(pts); closed=isclosed(c))
```
"""
function evolve_step!(pts::Vector{SVector{N,T}},
                      vel::Vector{SVector{N,T}},
                      c::AbstractDiscreteCurve{N,T},
                      law, τ::T) where {N,T}
    compute_velocities!(vel, c, law)
 
    max_v² = zero(T)
    @inbounds for i in eachindex(pts)
        v² = sum(abs2, vel[i])
        v² > max_v² && (max_v² = v²)
        pts[i] = pts[i] + τ * vel[i]
    end
    sqrt(max_v²)
end

"""
    evolve(c, law; kwargs...)  →  FlowResult
 
Evolve curve `c` under velocity `law` using explicit Euler time stepping.
 
**τ and σ**
 
- `τ` (timestep) controls the *size* of each step in physical time.
  Too large → numerical blow-up. Auto-CFL picks a safe value automatically.
- `σ` (step count, 0-indexed) counts *how many* steps were taken.
  `r[σ]` retrieves the curve after σ steps (if `save_every > 0`).
- Physical time = ∑ τₖ (variable if τ adapts). Access via `r.time`.
 
**Keyword arguments**
 
| kwarg             | default            | meaning                                             |
|-------------------|--------------------|-----------------------------------------------------|
| `τ`               | auto CFL           | Fixed timestep (nothing = CFL auto-selection)       |
| `cfl_safety`      | 0.25               | CFL safety factor. 0.5 = boundary of stability.     |
| `nsteps`          | 200                | Max steps σ. Actual may be less if converged.       |
| `tol`             | eps^(1/3)          | Stop when max\\|vᵢ\\| < tol.                        |
| `resample_every`  | 0                  | Resample to uniform arc-length every k steps.       |
| `adapt_τ_every`   | 50                 | Shrink τ if min edge collapses (0 = never).         |
| `save_every`      | 0                  | Save curve snapshot every k steps (0 = never).      |
| `track`           | []                 | Named functionals: `[:L => arc_length, ...]`        |
| `track_every`     | 1                  | Evaluate tracked functionals every k steps.         |
| `callback`        | nothing            | `f(σ, c, max_v) → Bool`. Return true to stop early.|
 
**Examples**
 
```julia
c = @curve (cos(t), sin(t)) t ∈ 0..2π  n=128  closed
 
# Minimal
r = evolve(c, CurvatureFlow())
r.curve          # final curve
 
# Full tracking + snapshots
r = evolve(c, CurvatureFlow();
    nsteps      = 500,
    save_every  = 10,        # save curve every 10 steps
    track       = [
        :L    => arc_length,
        :A    => enclosed_area,
        :κmax => c -> maximum(curvatures(c)),
    ])
 
r[0]        # initial curve
r[10]       # curve after 10 steps
r[100]      # curve after 100 steps
r[:L]       # Vector{Float64} of arc lengths at each tracked step
r[:A]       # Vector{Float64} of areas
 
# With tangential correction (prevents mesh degeneration)
r = evolve(c, CurvatureFlow(SteinerCurvature(); tangential=true); nsteps=2000)
 
# Custom law: gradient flow minimising bending energy
r = evolve(c, GradientFlow(total_squared_curvature); nsteps=300,
    track = [:E => total_squared_curvature])
 
# Composed law: shortening + attraction to a target
target = SVector(0.5, 0.0)
r = evolve(c, SumFlow(
    CurvatureFlow(),
    ScaledFlow(PointwiseFlow((c,i) -> target - vertex(c,i)), 0.05)
); nsteps=200)
 
# Your own custom law (just define compute_velocities!)
struct MyLaw end
DiscreteCurves.Flow.compute_velocities!(vel, c, ::MyLaw) =
    for i in 1:nvertices(c); vel[i] = curvature_normal_B(c, i); end
r = evolve(c, MyLaw(); nsteps=100)
 
# Manual loop for maximum control
pts    = copy(c.points)
vel    = similar(pts)
c_view = DiscreteCurve{2,Float64}(pts, true)
τ      = cfl_timestep(c)
for σ in 1:10_000
    max_v = evolve_step!(pts, vel, c_view, CurvatureFlow(), τ)
    max_v < 1e-9 && break
end
```
"""
function evolve(c0::AbstractDiscreteCurve{N,T},
                law;
                τ              :: Union{Nothing,Real} = nothing,
                cfl_safety     :: Real                = T(0.25),
                nsteps         :: Int                 = 200,
                tol            :: Real                = T(eps(T)^(1/3)),
                resample_every :: Int                 = 0,
                adapt_τ_every  :: Int                 = 50,
                save_every     :: Int                 = 0,
                track          :: AbstractVector      = Pair{Symbol}[],
                track_every    :: Int                 = 1,
                callback                             = nothing) where {N,T}
 
    closed = isclosed(c0)
    n      = nvertices(c0)
    tol_T  = T(tol)

    τ_cur  = τ === nothing ? cfl_timestep(c0; safety=T(cfl_safety)) : T(τ)
    τ_cur > zero(T) || throw(ArgumentError("τ must be positive, got $τ_cur"))

    pts    = copy(c0.points)
    vel    = Vector{SVector{N,T}}(undef, n)
    c_view = DiscreteCurve{N,T}(pts, closed)
 
    trackers  = [(Symbol(k), f) for (k,f) in track]
    history   = Dict{Symbol, AbstractVector}()
    for (sym, f) in trackers
        sample       = f(c0)
        V            = typeof(sample)
        history[sym] = Vector{V}([sample])
    end
    max_v_hist = T[]
 
    curves_saved = AbstractDiscreteCurve{N,T}[]
    times_saved  = T[]
    if save_every > 0
        push!(curves_saved, DiscreteCurve(copy(pts); closed))
        push!(times_saved,  zero(T))
    end
 
    total_time   = zero(T)
    actual_steps = 0
 
    @inbounds for σ in 1:nsteps
        actual_steps = σ

        max_v = evolve_step!(pts, vel, c_view, law, τ_cur)
 
        total_time += τ_cur
        push!(max_v_hist, max_v)
        if save_every > 0 && σ % save_every == 0
            push!(curves_saved, DiscreteCurve(copy(pts); closed))
            push!(times_saved, total_time)
        end
 
        if !isempty(trackers) && σ % track_every == 0
            for (sym, f) in trackers
                push!(history[sym], f(c_view))
            end
        end
        (callback !== nothing && callback(σ, c_view, max_v)) && break

        max_v < tol_T && break

        if resample_every > 0 && σ % resample_every == 0
            pts .= resample(c_view, n).points
        end
 
        if adapt_τ_every > 0 && σ % adapt_τ_every == 0
            τ_cur = min(τ_cur, _cfl_pts(pts, n, closed, T(cfl_safety)))
        end
    end

    history[:max_v] = max_v_hist
    final = DiscreteCurve(copy(pts); closed)
    FlowResult{N,T}(final, actual_steps, total_time,
                    max_v_hist, history, curves_saved, times_saved)
end

"""
    curve_shortening_flow(c, model=SteinerCurvature(); tangential=false, kwargs...)
 
Curvature-driven curve shortening. Forwards all kwargs to `evolve`.
 
```julia
curve_shortening_flow(c)                                   # κ^B default
curve_shortening_flow(c, TurningAngle())                   # κ^A
curve_shortening_flow(c, OsculatingCircle())               # κ^D
curve_shortening_flow(c; tangential=true, nsteps=2000)     # with reparam
curve_shortening_flow(c; save_every=10, track=[:L=>arc_length])
```
"""
function curve_shortening_flow(c, model=SteinerCurvature();
                                tangential::Bool=false, kw...)
    evolve(c, CurvatureFlow(model; tangential); kw...)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers
# ─────────────────────────────────────────────────────────────────────────────
 
@inline _adj(i::Int, n::Int, closed::Bool) = closed ? mod1(i, n) : clamp(i, 1, n)
 
@inline function _cfl_pts(pts::Vector{SVector{N,T}}, n::Int, closed::Bool, safety::T) where {N,T}
    min_ℓ² = typemax(T)
    ne = closed ? n : n-1
    @inbounds for i in 1:ne
        j  = closed ? mod1(i+1, n) : i+1
        ℓ² = sum(abs2, pts[j] - pts[i])
        ℓ² < min_ℓ² && (min_ℓ² = ℓ²)
    end
    safety * min_ℓ²
end
 
@inline function _tangential_v(d1::SVector{N,T}, d2::SVector{N,T}, n1::T, n2::T) where {N,T}
    # T̂ = unit forward tangent = d2/n2
    (iszero(n2)) && return zero(SVector{N,T})
    (n2 - n1) / 2 * (d2 / n2)
end

@inline _svec_setindex(v::SVector{N,T}, val::T, d::Int) where {N,T} =
    SVector{N,T}(ntuple(i -> i == d ? val : v[i], Val(N)))
 
end
