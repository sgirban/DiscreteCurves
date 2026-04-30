# Geometric Properties

Global scalar invariants of discrete curves: enclosed area, centroid, and
isoperimetric ratios.

---

## signed\_area

```julia
signed_area(c) → T
```

Signed area enclosed by a 2D curve using the **shoelace formula** (Green's theorem).

```math
A = \frac{1}{2} \sum_{i=1}^{n} \bigl( x_i y_{i+1} - x_{i+1} y_i \bigr)
```

**Sign convention:**
- ``A > 0`` — curve is oriented **counter-clockwise** (CCW).
- ``A < 0`` — curve is oriented **clockwise** (CW).
- ``A = 0`` — degenerate curve (self-intersecting or collinear points).

For **open** curves, the area is computed as if the endpoints were connected by a
straight line.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{2,T}` | A 2D curve (open or closed). |

### Returns

`T` — signed area. Units are length².

### Examples

```julia
# Unit circle: area ≈ π
c = circle_curve(1024)
signed_area(c)    # → ≈ 3.1416  (CCW, positive)

# CW circle: area ≈ -π
c_cw = orient(circle_curve(1024); to=CW)
signed_area(c_cw)   # → ≈ -3.1416

# Unit square: area = 1
c = square_curve(side=1.0; n=4)
abs(signed_area(c))   # → ≈ 1.0

# Orientation from sign alone
c = random_convex_polygon(32; seed=7)
if signed_area(c) > 0
    println("CCW")
else
    println("CW")
end

# Isoperimetric inequality check
c = fourier_curve(512; seed=3)
A = abs(signed_area(c))
L = arc_length(c)
r = 4π * A / L^2    # ≤ 1.0, equality for circle
println("IQ = $(round(r; digits=4))")

# Track area during curve-shortening flow
c = random_convex_polygon(64)
r = evolve(c, CurvatureFlow();
    nsteps = 500,
    track  = [:A => c -> abs(signed_area(c))])
r[:A]   # → shrinks toward 0
```

### Notes

- **Exact** for polygons (no quadrature error; the shoelace formula is exact).
- Works for **non-convex** and **self-intersecting** curves (though the geometric
  interpretation is then the algebraic area with winding numbers).
- For an unsigned area, use `abs(signed_area(c))`.

### See Also

[`centroid`](@ref), [`isoperimetric_quotient`](@ref), [`orient`](@ref), [`is_ccw`](@ref)

---

## centroid

```julia
centroid(c) → SVector{N,T}
```

Centroid (geometric centre / centre of mass) of a **closed** curve: the
arithmetic mean of all vertex positions.

```math
\bar{\gamma} = \frac{1}{n} \sum_{i=1}^{n} \gamma_i
```

Note: this is the centroid of the **vertex set**, not the centroid of the
enclosed area. These coincide only for very uniform meshes.

Throws `ArgumentError` for open curves.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | **Closed** curve required. |

### Returns

`SVector{N,T}` — the centroid position.

### Examples

```julia
# Unit circle centred at origin: centroid ≈ (0, 0)
c = circle_curve(128)
centroid(c)   # → SVector(≈0.0, ≈0.0)

# Shifted circle: centroid = shift
c = ParametricCurve(t -> SVector(1.0 + cos(t), 2.0 + sin(t)), 0, 2π, 128; closed=true)
centroid(c)   # → SVector(≈1.0, ≈2.0)

# Attraction toward centroid: PointwiseFlow
c = fourier_curve(256; seed=4)
ctr = centroid(c)
flow = PointwiseFlow((c, i) -> 0.1 * (centroid(c) - vertex(c, i)))
r = evolve(c, flow; nsteps=100)

# Monitor centroid displacement during Steiner flow
# (SteinerCurvature flow preserves centroid)
c0 = fourier_curve(256; seed=4)
r  = evolve(c0, CurvatureFlow(SteinerCurvature()); nsteps=500,
    track = [:cx => c -> centroid(c)[1],
             :cy => c -> centroid(c)[2]])
maximum(abs.(r[:cx] .- centroid(c0)[1]))   # → ≈ 0 (centroid preserved)
```

### Notes

For the centroid of the **enclosed area** (the true area centroid), use the
weighted formula:
```math
\bar{x}_A = \frac{1}{6A} \sum_i (x_i + x_{i+1})(x_i y_{i+1} - x_{i+1} y_i)
```
which is not directly provided but can be computed from the vertex data.

### See Also

[`signed_area`](@ref), [`SteinerCurvature`](@ref)

---

## isoperimetric\_quotient

```julia
isoperimetric_quotient(c) → T
```

The **isoperimetric quotient** (IQ) measures how "circle-like" a curve is:

```math
\text{IQ} = \frac{4\pi A}{L^2} \in (0, 1]
```

where ``A = |\text{signed\_area}(c)|`` and ``L = \text{arc\_length}(c)``.

- IQ = 1 if and only if the curve is a perfect circle.
- IQ closer to 0 indicates a long, thin, or highly indented curve.
- IQ is invariant under scaling, rotation, and translation.
- IQ requires the curve to be **closed** (returns 0 for degenerate cases).

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | Any closed 2D curve. |

### Returns

`T ∈ (0, 1]`.

### Examples

```julia
# Perfect circle: IQ = 1
c = circle_curve(1024)
isoperimetric_quotient(c)   # → ≈ 1.0

# Regular polygon: IQ increases with n
for n in [3, 4, 6, 8, 12, 32, 128]
    c  = regular_polygon(n)
    IQ = isoperimetric_quotient(c)
    println("n=$n: IQ = $(round(IQ; digits=4))")
end
# n=3:   0.6046  (equilateral triangle)
# n=4:   0.7854  (square)
# n=6:   0.9069  (hexagon)
# n=∞ → 1.0     (circle)

# Compare shapes
shapes = [
    ("circle",   circle_curve(256)),
    ("star",     star_curve(5)),
    ("fourier",  fourier_curve(256; seed=7)),
    ("square",   square_curve()),
]
for (name, c) in shapes
    println("$name: IQ = $(round(isoperimetric_quotient(c); digits=4))")
end

# Track IQ during curve-shortening flow (should approach 1)
c = random_concave_polygon(64; seed=3)
r = evolve(c, CurvatureFlow();
    nsteps = 1000,
    save_every = 20,
    track = [:IQ => isoperimetric_quotient])
r[:IQ]   # → increases toward 1.0 (approaches circle)
```

### See Also

[`signed_area`](@ref), [`arc_length`](@ref), [`regular_isoperimetric_quotient`](@ref)

---

## regular\_isoperimetric\_quotient

```julia
regular_isoperimetric_quotient(n::Int) → Float64
regular_isoperimetric_quotient(c)      → Float64
```

The **theoretical** isoperimetric quotient of a regular `n`-gon:

```math
\text{IQ}_n = \frac{\pi}{n} \cot\!\left(\frac{\pi}{n}\right)
```

This provides an upper bound on the IQ of any `n`-vertex polygon and a useful
baseline for comparing discrete curves.

When called with a curve `c`, uses `nvertices(c)` as `n`.

### Examples

```julia
regular_isoperimetric_quotient(3)    # → ≈ 0.6046  (triangle)
regular_isoperimetric_quotient(4)    # → ≈ 0.7854  (square)
regular_isoperimetric_quotient(6)    # → ≈ 0.9069  (hexagon)
regular_isoperimetric_quotient(12)   # → ≈ 0.9653
regular_isoperimetric_quotient(100)  # → ≈ 0.9997

# How far is a noisy polygon from the regular ideal?
c   = noisy_polygon(32; amplitude=0.3)
iq  = isoperimetric_quotient(c)
iqr = regular_isoperimetric_quotient(c)   # uses nvertices(c) = 32
efficiency = iq / iqr   # fraction of regular-polygon IQ achieved
println("Shape efficiency: $(round(100*efficiency; digits=1))%")
```

### See Also

[`isoperimetric_quotient`](@ref)