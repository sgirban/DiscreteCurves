# Curve Orientation System

## Overview

The DiscreteCurves library now includes comprehensive support for curve orientation, distinguishing between **counter-clockwise (CCW)** and **clockwise (CW)** curves. This is essential for mathematically consistent:

- Normal vector computations (inward vs outward)
- Flow simulations (curve shortening respects geometry)
- Signed area calculations
- Topological properties

## Mathematical Foundation

### Signed Area

For a 2D closed curve, the **shoelace formula** computes signed area:

$$A = \frac{1}{2} \sum_{i=0}^{n-1} (x_i y_{i+1} - x_{i+1} y_i)$$

- **Positive area** → CCW orientation
- **Negative area** → CW orientation
- **Zero area** → Degenerate curve

### Normal Vectors

For any 2D curve segment, given tangent direction $\vec{T}$:

- **Left normal** (rotate 90° CCW): $\vec{N}_L = (-T_y, T_x)$
- **Right normal** (rotate 90° CW): $\vec{N}_R = (T_y, -T_x)$

**Inward vs Outward depends on orientation:**

For **CCW curves:**
- Inward = left side = $-\vec{N}_L$
- Outward = right side = $\vec{N}_L$

For **CW curves:**
- Inward = right side = $\vec{N}_R$
- Outward = left side = $-\vec{N}_R$

## API

### Orientation Detection

```julia
using DiscreteCurves

# Create curves with different orientations
c_ccw = random_convex_polygon(8)  # Usually CCW
c_cw = orient(c_ccw; to=CW)       # Convert to CW

# Check orientation
is_ccw(c) :: Bool                 # true if counter-clockwise
is_cw(c) :: Bool                  # true if clockwise

# Get signed area
A = signed_area(c)  # positive=CCW, negative=CW
```

### Orientation Manipulation

```julia
# Non-mutating: returns new curve
c_ccw = orient(c_any; to=CCW)
c_cw = orient(c_any; to=CW)

# In-place: modifies curve
orient!(c; to=CCW)
orient!(c; to=CW)
```

### Normal Vectors

```julia
# Single normals
n_in = inward_normal(c, i)    # At vertex i, inward
n_out = outward_normal(c, i)  # At vertex i, outward

# All normals
ns_in = inward_normals(c)     # Vector of all inward normals
ns_out = outward_normals(c)   # Vector of all outward normals

# Properties
dot(inward_normal(c, i), outward_normal(c, i)) ≈ -1  # Opposite directions
```

## Integration with Existing Code

### Flow Simulations

**Current state:** Flows assume inward is negative curvature direction.

**Updated:**
- `CurvatureFlow()` now respects curve orientation
- Inward-pointing velocities depend on whether curve is CCW or CW
- Curve shortening still works correctly for both orientations

```julia
c_ccw = random_convex_polygon(20)
c_cw = orient(c_ccw; to=CW)

r_ccw = evolve(c_ccw, CurvatureFlow(SteinerCurvature()); nsteps=50)
r_cw = evolve(c_cw, CurvatureFlow(SteinerCurvature()); nsteps=50)

# Both converge to a point, but velocity directions differ based on orientation
```

### Curvature Calculations

**Note:** Unsigned curvature (magnitude) is orientation-independent. Signed curvature (`SignedCurvature2D`) depends on orientation:

```julia
# Unsigned curvature: same regardless of orientation
κ_unsigned = curvatures(c, SteinerCurvature())

# Signed curvature: changes sign with orientation
c_ccw = circle_curve(32)
c_cw = orient(c_ccw; to=CW)

κ_ccw = curvatures(c_ccw, SignedCurvature2D())
κ_cw = curvatures(c_cw, SignedCurvature2D())

# κ_cw ≈ -κ_ccw (opposite signs for CW vs CCW)
```

## Performance Considerations

### 2D vs Higher Dimensions

- **2D:** Full orientation support (signed area, inward/outward normals)
- **3D+:** Limited support (normals require additional frame specification)

For 3D+ curves, we recommend:
1. Work with 2D projections when orientation matters
2. Use explicit normal specification (Frenet frame, chosen orientation)
3. Store orientation as metadata if needed

### Computational Cost

- `signed_area(c)`: **O(n)** (one pass through vertices)
- `is_ccw(c)`, `is_cw(c)`: **O(n)** (calls `signed_area`)
- `inward_normal(c, i)`: **O(1)** (two neighbors + orientation check)
- `orient(c; to)`: **O(n)** (reversal if needed) / **O(1)** (if already correct)
- `orient!(c; to)`: **O(n)** (in-place reversal) / **O(1)** (if already correct)

## Examples

### Example 1: Ensure CCW Orientation

```julia
using DiscreteCurves, GLMakie

# Generate random polygon (orientation unknown)
c = random_convex_polygon(12; seed=42)

# Ensure CCW
orient!(c; to=CCW)

# Now safe to use with flow/geometry that assumes CCW
r = evolve(c, CurvatureFlow(); nsteps=100)
```

### Example 2: Compare CCW vs CW Curves

```julia
c_ccw = circle_curve(64)
c_cw = orient(c_ccw; to=CW)

# Inward normals point in opposite directions
n_in_ccw = inward_normals(c_ccw)
n_in_cw = inward_normals(c_cw)

# For unit circle, these should be opposite
println(maximum(norm.(n_in_ccw .+ n_in_cw)))  # Should be close to 0
```

### Example 3: Geometry Verification

```julia
# Verify convex polygon has consistent inward normals
c = random_convex_polygon(20)
centroid = sum(c.points) / nvertices(c)

all_correct = all(1:nvertices(c)) do i
    to_center = centroid - c.points[i]
    inward_n = inward_normal(c, i)
    dot(inward_n, to_center) > -0.01  # Some tolerance for numerical error
end

println("All inward normals point toward centroid: $all_correct")
```

## Future Enhancements

1. **3D Orientation:** Implement using signed volume and Frenet frames
2. **Orientation Field:** Store orientation per-vertex for mixed-orientation curves
3. **Oriented Flows:** Make flow direction explicitly depend on orientation choice
4. **Oriented Curvature:** Full signed curvature field in all dimensions
5. **Topology Detection:** Use orientation to detect curve self-intersections

## Testing

Run the comprehensive orientation test suite:

```bash
julia test/test_orientation.jl
```

Tests verify:
- ✓ Signed area computation
- ✓ CCW/CW detection
- ✓ Orientation conversion (mutating and non-mutating)
- ✓ Normal vector consistency
- ✓ Normal opposition property
- ✓ Geometric correctness (inward pointing)
