using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DiscreteCurves
using LinearAlgebra
using StaticArrays
using Random

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: Circle with different precisions
# ─────────────────────────────────────────────────────────────────────────────

println("=" ^ 80)
println("Test 1: CIRCLE (should have sum ≈ 2π)")
println("=" ^ 80)

n_points = 64

for FloatType in [Float32, Float64, BigFloat]
    println("\n$(repr(FloatType)):")
    
    # Create a unit circle
    θs = range(0, 2π, length=n_points+1)[1:end-1]
    points = [SVector{2,FloatType}(cos(θ), sin(θ)) for θ in θs]
    
    # Create a closed curve
    c = ClosedCurve(points)
    
    # Compute signed turning angles using the unified API
    angles_signed = turning_angles(c; signed=true)
    sum_signed = sum(angles_signed)
    
    # Expected: 2π for CCW circle
    theoretical_2π = convert(FloatType, 2 * π)
    error_signed = abs(sum_signed - theoretical_2π)
    
    println("  Signed turning angles (sign precision test):")
    println("    Sum:       $sum_signed")
    println("    Expected:  $theoretical_2π")
    println("    Error:     $error_signed")
    println("    Rel.Error: $(error_signed / theoretical_2π)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Square (4 right angles = 2π total)
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 2: SQUARE (4 corners × 90° = 2π)")
println("=" ^ 80)

for FloatType in [Float32, Float64, BigFloat]
    println("\n$(repr(FloatType)):")
    
    # Square corners with extra points on edges
    points_base = [
        SVector{2,FloatType}(0, 0),
        SVector{2,FloatType}(1, 0),
        SVector{2,FloatType}(1, 1),
        SVector{2,FloatType}(0, 1),
    ]
    
    # Add intermediate points for better sampling
    points = SVector{2,FloatType}[]
    for i in eachindex(points_base)
        curr = points_base[i]
        next = points_base[mod1(i+1, length(points_base))]
        push!(points, curr)
        # Add 2 intermediate points
        for t in [1/3, 2/3]
            push!(points, curr + t * (next - curr))
        end
    end
    
    c = ClosedCurve(points)
    
    angles_signed = turning_angles(c; signed=true)
    sum_signed = sum(angles_signed)
    
    theoretical_2π = convert(FloatType, 2 * π)
    theoretical_π_half = convert(FloatType, π / 2)
    
    println("  Turning angles (2D: signed):")
    println("    Sum:       $sum_signed")
    println("    Expected:  $theoretical_2π")
    println("    Error:     $(abs(sum_signed - theoretical_2π))")
    
    # Check that we have 4 significant jumps (90 degrees each) at corners
    corner_indices = [1,2,3,4]  # Approximately where corners are
    println("  Turning angles at approximate corners:")
    for idx in corner_indices
        if idx <= length(angles_signed)
            angle_deg = angles_signed[idx] * 180 / π
            println("    angles[$idx] ≈ $(angle_deg)°")
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Noisy circle (test numerical stability)
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 3: NOISY CIRCLE (numerical stability test)")
println("=" ^ 80)

for FloatType in [Float32, Float64]
    println("\n$(repr(FloatType)):")
    
    local n_points = 32
    θs = range(0, 2π, length=n_points+1)[1:end-1]
    
    # Add small noise to circle
    noise_level = convert(FloatType, 0.01)
    Random.seed!(42)
    points = [
        SVector{2,FloatType}(
            cos(θ) + noise_level * randn(FloatType),
            sin(θ) + noise_level * randn(FloatType)
        ) for θ in θs
    ]
    
    # Normalize back to approximate unit circle
    points = [p / norm(p) for p in points]
    
    c = ClosedCurve(points)
    
    angles_signed = turning_angles(c; signed=true)
    sum_signed = sum(angles_signed)
    
    theoretical_2π = convert(FloatType, 2 * π)
    error = abs(sum_signed - theoretical_2π)
    
    println("  With $(FloatType):")
    println("    Sum:       $sum_signed")
    println("    Expected:  $theoretical_2π")
    println("    Error:     $error")
    println("    Rel.Error: $(error / theoretical_2π)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: Curvature dispatch system
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 4: CURVATURE DISPATCH SYSTEM")
println("=" ^ 80)

n_points = 32
θs = range(0, 2π, length=n_points+1)[1:end-1]
points = [SVector{2,Float64}(cos(θ), sin(θ)) for θ in θs]
c = ClosedCurve(points)

println("\nCircle (unit circle has κ = 1 everywhere):")

try
    κ_turning = curvature(c, TurningAngle())
    println("  TurningAngle mean κ:       $(mean(κ_turning))")
    println("  Expected (≈ 1):            1.0")
catch e
    println("  ERROR: $e")
end

try
    κ_signed = curvature(c, SignedCurvature2D())
    println("  SignedCurvature2D mean κ:  $(mean(κ_signed))")
    println("  Angles sum to:             $(sum(turning_angles(c; signed=true)))")
    println("  Expected (≈ 2π):           $(2π)")
catch e
    println("  ERROR: $e")
end

try
    κ_steiner = curvature(c, SteinerCurvature())
    println("  SteinerCurvature mean κ:   $(mean(κ_steiner))")
catch e
    println("  ERROR: $e")
end

println("\n" * "=" ^ 80)
println("All precision tests completed!")
println("=" ^ 80)
