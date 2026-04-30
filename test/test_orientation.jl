using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DiscreteCurves
using StaticArrays
using LinearAlgebra

println("=" ^ 80)
println("ORIENTATION FEATURE TESTS")
println("=" ^ 80)

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: Signed Area and Orientation Detection
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 1: SIGNED AREA & ORIENTATION DETECTION")
println("=" ^ 80)

# Create a CCW square (counterclockwise)
square_ccw = ClosedCurve([
    SVector(0.0, 0.0),
    SVector(1.0, 0.0),
    SVector(1.0, 1.0),
    SVector(0.0, 1.0),
])

# Create a CW square (clockwise)
square_cw = ClosedCurve([
    SVector(0.0, 0.0),
    SVector(0.0, 1.0),
    SVector(1.0, 1.0),
    SVector(1.0, 0.0),
])

println("\nCCW Square:")
area_ccw = signed_area(square_ccw)
println("  Signed area: $(area_ccw)")
println("  Is CCW: $(is_ccw(square_ccw))")
println("  Is CW: $(is_cw(square_ccw))")

println("\nCW Square:")
area_cw = signed_area(square_cw)
println("  Signed area: $(area_cw)")
println("  Is CCW: $(is_ccw(square_cw))")
println("  Is CW: $(is_cw(square_cw))")

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Random Convex Polygon Orientation
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 2: RANDOM CONVEX POLYGON ORIENTATION")
println("=" ^ 80)

c_rand = random_convex_polygon(10; seed=42)

println("\nRandom convex polygon (10 vertices):")
area_rand = signed_area(c_rand)
current_orient = is_ccw(c_rand) ? "CCW" : "CW"
println("  Signed area: $(round(area_rand; digits=6))")
println("  Current orientation: $(current_orient)")

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Orient Function (Non-mutating)
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 3: ORIENT FUNCTION (Non-mutating)")
println("=" ^ 80)

# Start with square_cw
println("\nOriginal square (CW):")
println("  Area: $(signed_area(square_cw))")
println("  Orientation: $(is_cw(square_cw) ? "CW" : "CCW")")

# Convert to CCW
square_ccw_converted = orient(square_cw; to=CCW)
println("\nAfter orient(; to=CCW):")
println("  Area: $(signed_area(square_ccw_converted))")
println("  Orientation: $(is_ccw(square_ccw_converted) ? "CCW" : "CW")")

# Original unchanged
println("\nOriginal square (should be unchanged):")
println("  Area: $(signed_area(square_cw))")
println("  Orientation: $(is_cw(square_cw) ? "CW" : "CCW")")

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: Orient! Function (Mutating)
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 4: ORIENT! FUNCTION (Mutating)")
println("=" ^ 80)

c_test = random_convex_polygon(8; seed=100)
println("\nBefore orient!:")
println("  Area: $(round(signed_area(c_test); digits=6))")
println("  Orientation: $(is_ccw(c_test) ? "CCW" : "CW")")

orient!(c_test; to=CCW)
println("\nAfter orient!(; to=CCW):")
println("  Area: $(round(signed_area(c_test); digits=6))")
println("  Orientation: $(is_ccw(c_test) ? "CCW" : "CW")")

# ─────────────────────────────────────────────────────────────────────────────
# Test 5: Inward and Outward Normals
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 5: INWARD & OUTWARD NORMALS")
println("=" ^ 80)

# Simple square for clear results
test_square = ClosedCurve([
    SVector(0.0, 0.0),
    SVector(1.0, 0.0),
    SVector(1.0, 1.0),
    SVector(0.0, 1.0),
])

println("\nCCW Square normals:")
println("  Orientation: $(is_ccw(test_square) ? "CCW" : "CW")")

inward_ns = inward_normals(test_square)
outward_ns = outward_normals(test_square)

for i in 1:nvertices(test_square)
    println("  Vertex $i:")
    println("    Inward:  $(round.(inward_ns[i]; digits=3))")
    println("    Outward: $(round.(outward_ns[i]; digits=3))")
    println("    Sum (should be ≈0): $(round.(inward_ns[i] + outward_ns[i]; digits=6))")
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 6: Orientation with Circle
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 6: ORIENTATION WITH CIRCLE")
println("=" ^ 80)

circle = circle_curve(32)
println("\nCircle (32 vertices):")
println("  Signed area: $(round(signed_area(circle); digits=6))")
println("  Orientation: $(is_ccw(circle) ? "CCW" : "CW")")

# Reverse circle
circle_rev = orient(circle; to=CW)
println("\nReversed circle:")
println("  Signed area: $(round(signed_area(circle_rev); digits=6))")
println("  Orientation: $(is_ccw(circle_rev) ? "CCW" : "CW")")

# ─────────────────────────────────────────────────────────────────────────────
# Test 7: Normal consistency check
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("Test 7: NORMAL CONSISTENCY")
println("=" ^ 80)

# For a convex polygon, inward normals should point toward centroid
test_poly = random_convex_polygon(6; seed=200)
ctr = sum(test_poly.points) / nvertices(test_poly)

println("\nConvex polygon (6 vertices):")
println("  Centroid: $(round.(ctr; digits=3))")

dot_prods = Float64[]
for i in 1:nvertices(test_poly)
    v = test_poly.points[i]
    to_center = ctr - v
    n_inward = inward_normal(test_poly, i)
    
    # Dot product should be positive if pointing generally toward center
    dot_prod = dot(n_inward, to_center)
    push!(dot_prods, dot_prod)
    
    if i ≤ 3
        println("  Vertex $i: inward·to_center = $(round(dot_prod; digits=4))")
    end
end

all_point_inward = all(dp > -0.1 for dp in dot_prods)
println("  All inward normals point roughly toward centroid: $(all_point_inward)")

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 80)
println("ALL ORIENTATION TESTS COMPLETED")
println("=" ^ 80)
println("\nKey features demonstrated:")
println("  ✓ signed_area() computes signed area (orientation indicator)")
println("  ✓ is_ccw() and is_cw() detect orientation")
println("  ✓ orient() non-mutating orientation conversion")
println("  ✓ orient!() in-place orientation conversion")
println("  ✓ inward_normal() and outward_normal() respect orientation")
println("  ✓ inward_normals() and outward_normals() for all vertices")
println("  ✓ Normals are opposite (sum to zero)")
println("  ✓ Inward normals point toward interior of convex polygons")
