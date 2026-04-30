# Lengths

Functions for computing arc length, individual edge lengths, and cumulative
arc-length parameters along a curve.

---

## arc\_length

```julia
arc_length(c) → T
```

Total arc length of the curve: the sum of all edge lengths.

```math
L = \sum_{i=1}^{\text{nedges}(c)} \| \gamma_{i+1} - \gamma_i \|
```

**Complexity**: O(n), single pass, zero allocations.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | Any discrete curve, open or closed, any dimension. |

### Returns

`T` — total arc length in the same units as the vertex coordinates.

### Examples

```julia
# Unit circle — arc length ≈ 2π
c = circle_curve(1024)
arc_length(c)   # → ≈ 6.2832  (converges to 2π as n → ∞)

# Ellipse with semi-axes 2 and 1
c = ellipse_curve(512; a=2.0, b=1.0)
arc_length(c)   # → ≈ 9.6884  (≈ elliptic integral value)

# Open curve: parabola y = x² from x=-1 to x=1
xs = range(-1.0, 1.0; length=256)
c  = OpenCurve([SVector(x, x^2) for x in xs])
arc_length(c)   # → ≈ 2.958

# 3D helix (one turn of radius 1, pitch 1/(2π))
c = ParametricCurve(t -> SVector(cos(t), sin(t), t/(2π)), 0, 2π, 512)
arc_length(c)   # → ≈ √(1 + (1/2π)²) * 2π ≈ 6.366

# Isoperimetric ratio check (circle maximises area/perimeter²)
c   = random_convex_polygon(64)
L   = arc_length(c)
A   = abs(signed_area(c))
r   = 4π * A / L^2    # ≤ 1, equals 1 only for a circle
println("IQ = $(round(r; digits=4))")

# Track arc length during a flow
r = evolve(c, CurvatureFlow();
    nsteps = 200,
    track  = [:L => arc_length])
r[:L]   # → Vector of arc lengths, one per step
```

### Notes

- For a **closed** curve, the closing edge (vertex `n` → vertex `1`) is included.
- Arc length is a **Riemannian invariant**: it depends only on the curve, not on
  how it is parametrised or oriented.

### See Also

[`edge_lengths`](@ref), [`cumulative_length`](@ref), [`isoperimetric_quotient`](@ref)

---

## edge\_lengths

```julia
edge_lengths(c) → Vector{T}
```

Return a vector of all `nedges(c)` individual edge lengths.

```math
\ell_i = \|\gamma_{i+1} - \gamma_i\|, \quad i = 1,\ldots,\text{nedges}(c)
```

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | Any discrete curve. |

### Returns

`Vector{T}` of length `nedges(c)`. Entry `i` is the length of edge `i`
(from vertex `i` to vertex `i+1`, with wrap for closed curves).

### Examples

```julia
c = circle_curve(8)
ℓ = edge_lengths(c)
# → 8-element Vector{Float64}, all ≈ 2sin(π/8) ≈ 0.765

# Uniformity check: coefficient of variation should be 0 for a regular polygon
cv = std(ℓ) / mean(ℓ)
println("Edge length CV = $cv")   # → 0.0 for regular_polygon, > 0 for noisy curves

# Attach as EdgeData for topology-aware access
c = fourier_curve(256)
ℓ = edge_lengths(c)
ed = attach_edge_data(c, ℓ)
ed[3]    # → third edge length

# CFL timestep is proportional to min edge length squared
τ = 0.25 * minimum(edge_lengths(c))^2

# Dual edge lengths (used in curvature computation)
ℓ = edge_lengths(c)
ℓ_dual = [(ℓ[mod1(i-1,nedges(c))] + ℓ[i])/2 for i in 1:nedges(c)]

# Sum = arc length
sum(edge_lengths(c)) ≈ arc_length(c)   # → true
```

### Notes

Allocates one `Vector{T}` of size `nedges(c)`. If you only need the total arc
length, use [`arc_length`](@ref) instead (zero allocations).

### See Also

[`arc_length`](@ref), [`cumulative_length`](@ref), [`attach_edge_data`](@ref)

---

## cumulative\_length

```julia
cumulative_length(c) → Vector{T}
```

Arc-length parameter `s[i]` at each vertex, starting from zero at the first
vertex. This is the natural arc-length parametrisation of the discrete curve.

```math
s[1] = 0, \quad s[i+1] = s[i] + \|\gamma_{i+1} - \gamma_i\|
```

The returned vector has length `nvertices(c)`, so `s[end]` equals the total
arc length of all edges **except** the closing edge of a closed curve.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | Any discrete curve. |

### Returns

`Vector{T}` of length `nvertices(c)`. Always starts with `0.0`.

### Examples

```julia
c = circle_curve(128)
s = cumulative_length(c)
s[1]     # → 0.0
s[end]   # → arc length up to (but not including) the closing edge
         # ≈ 2π * (127/128)

# Normalised arc-length parameter ∈ [0, 1)
s_norm = s ./ arc_length(c)
s_norm[1]    # → 0.0
s_norm[end]  # → ≈ 127/128

# Attach as vertex data for colour mapping
c  = fourier_curve(256)
s  = cumulative_length(c)
vd = attach_vertex_data(c, s)     # VertexData{Float64}
# or directly in a bundle
b  = bundle(c)
b[:arc_length] = cumulative_length(c)

# Arc-length uniform resampling (find parameter at equally-spaced s values)
s    = cumulative_length(c)
L    = arc_length(c)
s_eq = range(0, L; length=nvertices(c)+1)[1:end-1]
# Interpolate curve at s_eq values using your preferred method

# Use as x-axis in a curvature plot
using Plots
s  = cumulative_length(c)
κs = curvatures(c)
# For closed curves, κs has length nvertices(c); for open, nvertices-2
plot(s[2:end-1], κs)   # interior vertices for open curve
```

### Notes

- For a **closed** curve, the closing edge is **not** included in `s[end]`.
  Use `s[end] + edge_length(c, nedges(c))` if you need the full perimeter.
- Complexity: O(n), single pass.

### See Also

[`arc_length`](@ref), [`edge_lengths`](@ref)