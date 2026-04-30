using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DiscreteCurves
using LinearAlgebra
using StaticArrays
using Statistics

println("=" ^ 80)
println("FLOW COMPUTATION TESTS: Verification with known curves")
println("=" ^ 80)

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: Unit Circle Curvature Flow
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 1: UNIT CIRCLE - Curvature Flow Velocity Verification")
println("=" ^ 80)

# Create unit circle with 32 vertices
n_circle = 32
θs = range(0, 2π, length=n_circle+1)[1:end-1]
circle_points = [SVector{2,Float64}(cos(θ), sin(θ)) for θ in θs]
circle = ClosedCurve(circle_points)

println("\nCircle properties:")
println("  Vertices: $(nvertices(circle))")
println("  Closed: $(isclosed(circle))")
println("  Arc length: $(arc_length(circle)) (expected: ≈ 2π = $(2π))")

# For a unit circle, curvature should be 1 everywhere
κs_steiner = curvatures(circle, SteinerCurvature())
κs_turning = curvatures(circle, TurningAngle())

println("\nCurvatures:")
println("  SteinerCurvature:   mean=$(round(mean(κs_steiner); digits=6)), std=$(round(std(κs_steiner); digits=6))")
println("  TurningAngle:       mean=$(round(mean(κs_turning); digits=6)), std=$(round(std(κs_turning); digits=6))")
println("  Expected (≈1.0):    1.0")

# Compute velocities for curve shortening flow (should point inward)
vel_steiner = similar(circle_points)
compute_velocities!(vel_steiner, circle, CurvatureFlow(SteinerCurvature()))

println("\nVelocity verification (should point toward center):")
for i in 1:min(4, nvertices(circle))
    v = vel_steiner[i]
    v_norm = norm(v)
    # For unit circle, velocity should point toward origin (inward)
    p = circle_points[i]
    direction_to_center = -p  # From point to origin
    direction_to_center_norm = normalize(direction_to_center)
    v_normalized = v_norm > 0 ? normalize(v) : v
    
    # Dot product should be close to 1 (same direction)
    if v_norm > 1e-10
        dot_prod = dot(v_normalized, direction_to_center_norm)
        println("  Vertex $i: |v|=$(round(v_norm; digits=6)), direction dot=$(round(dot_prod; digits=6)) (expected ≈ 1)")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Straight Line (Zero Curvature)
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 2: STRAIGHT LINE - Should have zero curvature")
println("=" ^ 80)

# Create a straight line (open curve)
line_points = [SVector{2,Float64}(i/10, 0) for i in 0:10]
line = OpenCurve(line_points)

println("\nLine properties:")
println("  Vertices: $(nvertices(line))")
println("  Closed: $(isclosed(line))")
println("  Arc length: $(arc_length(line)) (expected: ≈ 1.0)")

κs_line = curvatures(line, SteinerCurvature())
println("\nCurvatures (should be ≈ 0):")
println("  Interior vertices: $(κs_line)")
    println("  Mean: $(round(mean(κs_line); digits=8))")
    println("  Max: $(round(maximum(abs.(κs_line)); digits=8))")
vel_line = similar(line_points)
compute_velocities!(vel_line, line, CurvatureFlow(SteinerCurvature()))

println("\nVelocities (should be ≈ 0):")
for i in 1:nvertices(line)
    println("  Vertex $i: $(vel_line[i]) (norm=$(norm(vel_line[i])))")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: One Euler Step Evolution
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 3: ONE EULER STEP - Check convergence toward equilibrium")
println("=" ^ 80)

# Start with circle
c0 = ClosedCurve(circle_points)
pts_0 = copy(c0.points)

# Compute one small Euler step
τ = 0.001  # Small timestep
pts_1 = copy(pts_0)
vel = Vector{SVector{2,Float64}}(undef, nvertices(c0))
c_view = DiscreteCurve(pts_1, isclosed(c0))

# Compute velocities at initial position
compute_velocities!(vel, c_view, CurvatureFlow(SteinerCurvature()))

# Manual Euler step: x_{n+1} = x_n + τ * v_n
for i in eachindex(pts_1)
    pts_1[i] = pts_1[i] + τ * vel[i]
end

c1 = DiscreteCurve(pts_1; closed=isclosed(c0))
L0 = arc_length(c0)
L1 = arc_length(c1)

println("\nAfter one Euler step (τ=$τ):")
println("  Initial arc length:  L₀ = $(L0)")
println("  Final arc length:    L₁ = $(L1)")
println("  Change:              ΔL = $(L1 - L0)")
println("  Relative change:     ΔL/L₀ = $((L1-L0)/L0)")
println("  (Should be negative for curve shortening)")

# Check max velocity decreased
v_max_0 = maximum(norm.(vel))
println("  Max velocity magnitude: $(v_max_0)")
println("  Expected: ≈ $(mean(κs_steiner)) = $(mean(κs_steiner))")

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: Full Evolution Comparison
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 4: FULL EVOLUTION - Track arc length decrease")
println("=" ^ 80)

c_test = ClosedCurve(circle_points)
r = evolve(c_test, CurvatureFlow(SteinerCurvature());
           nsteps=100,
           cfl_safety=0.1,
           save_every=10,
           track=[:L => arc_length, :κmax => c -> maximum(curvatures(c, SteinerCurvature()))])

println("\nEvolution results:")
println("  Steps taken: $(r.nsteps)")
println("  Physical time: $(r.time)")
println("  Initial arc length: $(first(r[:L]))")
println("  Final arc length: $(last(r[:L]))")
println("  Length reduction: $(100 * (1 - last(r[:L])/first(r[:L])))%")

println("\nArc length at saved steps:")
for (i, L) in enumerate(r[:L])
    step = (i-1) * 10
    println("  Step $step: L = $(L)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 5: Different Curvature Models
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 5: CURVATURE MODELS COMPARISON")
println("=" ^ 80)

c_comp = ClosedCurve(circle_points)

# Compare different models
models = [
    ("TurningAngle", TurningAngle()),
    ("SteinerCurvature", SteinerCurvature()),
    ("OsculatingCircle", OsculatingCircle()),
]

for (name, model) in models
    κs = curvatures(c_comp, model)
    println("\n$name:")
    println("  Mean curvature:  $(mean(κs))")
    println("  Std dev:         $(std(κs))")
    println("  Min/Max:         $(minimum(κs)) / $(maximum(κs))")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 6: Velocity Computation Correctness
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 6: VELOCITY COMPUTATION - Manual verification")
println("=" ^ 80)

# Create a simple curve
simple_points = [
    SVector{2,Float64}(0, 0),
    SVector{2,Float64}(1, 0),
    SVector{2,Float64}(1, 1),
    SVector{2,Float64}(0, 1),
]
simple = ClosedCurve(simple_points)

println("\nSimple square:")
for i in 1:nvertices(simple)
    println("  Vertex $i: $(simple_points[i])")
end

# Compute velocities
vel_simple = similar(simple_points)
compute_velocities!(vel_simple, simple, CurvatureFlow(SteinerCurvature()))

println("\nVelocities (curvature flow):")
for i in 1:nvertices(simple)
    κ = curvature(simple, i, SteinerCurvature())
    println("  Vertex $i: κ=$(round(κ; digits=6)), v=$(vel_simple[i])")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 7: CFL Timestep
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 7: CFL TIMESTEP SELECTION")
println("=" ^ 80)

c_cfl = ClosedCurve(circle_points)
τ_auto = cfl_timestep(c_cfl)
edge_lens = edge_lengths(c_cfl)
min_edge = minimum(edge_lens)

println("\nCircle CFL timestep:")
println("  Min edge length: $(min_edge)")
println("  τ_auto: $(τ_auto)")
println("  Expected: ≈ 0.25 × $(min_edge)² = $(0.25 * min_edge^2)")

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("ALL TESTS COMPLETED")
println("=" ^ 80)
