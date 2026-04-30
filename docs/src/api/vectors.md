# Tangents & Normals

Functions for computing tangent and normal vectors at vertices and along edges.
The library distinguishes between **edge-level** quantities (one per edge) and
**vertex-level** quantities (one per vertex), and provides both raw and
orientation-aware variants.

---

## Overview

| Function | Level | Returns | Notes |
|----------|-------|---------|-------|
| `edge_tangent(c, i)` | Edge | Unit tangent of edge `i` | Direction: vertex `i` → `i+1` |
| `edge_tangents(c)` | Edge | All unit edge tangents | Vector of length `nedges(c)` |
| `edge_normal(c, i)` | Edge | Edge tangent rotated 90° CCW | 2D only |
| `edge_normals(c)` | Edge | All edge normals | 2D only |
| `tangent_vector(c, i)` | Vertex | Angle-bisector of adjacent edges | Both endpoints handled for open curves |
| `tangent_vectors(c)` | Vertex | All vertex tangents | Length `nvertices(c)` |
| `normal_vector(c, i)` | Vertex | `normalize(t_out - t_in)` | Points toward concave side |
| `normal_vectors(c)` | Vertex | All vertex normals | Length `nvertices(c)` |
| `inward_normal(c, i)` | Vertex | Orientation-corrected inward normal | Closed curves only |
| `outward_normal(c, i)` | Vertex | `-inward_normal(c, i)` | Closed curves only |
| `inward_normals(c)` | Vertex | All inward normals | Closed curves only |
| `outward_normals(c)` | Vertex | All outward normals | Closed curves only |

---

## edge\_tangent / edge\_tangents

```julia
edge_tangent(c, i)  → SVector{N,T}
edge_tangents(c)    → Vector{SVector{N,T}}
```

Unit tangent of edge `i`: the normalised direction from `vertex(c, i)` toward
`vertex(c, i+1)`.

For closed curves, edge `nedges(c)` wraps from the last vertex back to the
first.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | Any curve, any dimension. |
| `i` | `Int` | 1-based edge index in `1:nedges(c)`. |

### Returns

`SVector{N,T}` — unit vector (length 1, except at degenerate zero-length edges
where the zero vector is returned).

### Examples

```julia
c = circle_curve(4)   # square-like, 4 vertices

t1 = edge_tangent(c, 1)   # tangent of first edge
norm(t1)                   # → 1.0

# Verify: consecutive edge tangents differ by the turning angle
t1 = edge_tangent(c, 1)
t2 = edge_tangent(c, 2)
acos(clamp(dot(t1, t2), -1, 1))   # ≈ π/2 for a square

# All edge tangents
ts = edge_tangents(c)      # Vector of nedges(c) unit vectors

# Batch: check all tangents are unit vectors
all(norm(t) ≈ 1.0 for t in edge_tangents(c))  # → true
```

### See Also

[`tangent_vector`](@ref), [`edge_normal`](@ref)

---

## edge\_normal / edge\_normals

```julia
edge_normal(c, i)  → SVector{2,T}
edge_normals(c)    → Vector{SVector{2,T}}
```

Unit normal of edge `i`: the edge tangent rotated **90° counter-clockwise**
(CCW). Defined for **2D curves only**.

**Orientation note (2D closed curves):**
- **CCW curve** → edge normal points *left* of travel → **inward** for convex shapes.
- **CW curve** → edge normal points *right* of travel → **outward** for convex shapes.

For an orientation-independent inward normal at vertices, use [`inward_normal`](@ref).

### Examples

```julia
c = circle_curve(128)   # CCW unit circle by default

n1 = edge_normal(c, 1)
t1 = edge_tangent(c, 1)
dot(t1, n1)   # → ≈ 0.0  (perpendicular)

# For a CCW circle, edge normal points inward (toward origin)
# Midpoint of edge 1:
p, q  = edge_segment(c, 1)
mid   = (p + q) / 2
# Edge normal points roughly toward origin
dot(n1, -mid/norm(mid))   # → > 0

# All edge normals
ns = edge_normals(c)   # length nedges(c)
```

### See Also

[`edge_tangent`](@ref), [`inward_normal`](@ref), [`normal_vector`](@ref)

---

## tangent\_vector / tangent\_vectors

```julia
tangent_vector(c, i)  → SVector{N,T}
tangent_vectors(c)    → Vector{SVector{N,T}}
```

Unit tangent at vertex `i`: the **angle bisector** of the incoming and outgoing
edge directions.

```math
\hat{T}_i = \operatorname{normalize}\!\bigl(\hat{t}_{i-1} + \hat{t}_i\bigr)
```

where ``\hat{t}_{i-1}`` and ``\hat{t}_i`` are the unit tangents of the incoming
and outgoing edges. For uniform meshes this coincides with the centred difference
``(\gamma_{i+1} - \gamma_{i-1})``; for non-uniform meshes it gives equal angular
weight to both adjacent edges.

**Boundary handling (open curves):**
- Vertex `1`: returns `edge_tangent(c, 1)` (only outgoing edge).
- Vertex `n`: returns `edge_tangent(c, n-1)` (only incoming edge).

Returns the **zero vector** at hairpin vertices (180° turn) where both edge
tangents are anti-parallel.

### Examples

```julia
c = circle_curve(256)

T5 = tangent_vector(c, 5)
norm(T5)   # → 1.0

# Verify: tangent is perpendicular to the inward normal
n5 = inward_normal(c, 5)
dot(T5, n5)   # → ≈ 0.0

# All vertex tangents
Ts = tangent_vectors(c)   # length nvertices(c)

# Open curve: boundary behaviour
c_open = OpenCurve([SVector(x, x^2) for x in range(-1.0, 1.0; length=64)])
tangent_vector(c_open, 1)          # = edge_tangent(c_open, 1)
tangent_vector(c_open, 64)         # = edge_tangent(c_open, 63)
tangent_vector(c_open, 32)         # centred bisector

# Visualise tangents
using GLMakie
curveplot(c) + vectorplot(tangent_vectors(c); anchors=c,
                          arrow_color=:royalblue, normalize=true)
```

### See Also

[`edge_tangent`](@ref), [`normal_vector`](@ref), [`inward_normal`](@ref)

---

## normal\_vector / normal\_vectors

```julia
normal_vector(c, i)  → SVector{N,T}
normal_vectors(c)    → Vector{SVector{N,T}}
```

Unit normal at vertex `i`, defined as:

```math
\hat{N}_i = \operatorname{normalize}\!\bigl(\hat{t}_i - \hat{t}_{i-1}\bigr)
```

This vector **points toward the concave side** of the curve (inward for convex
shapes). Its direction is parallel to (and has the same sign as) the curvature
vector.

For closed curves this is identical to [`inward_normal`](@ref) when the curve is
CCW; for CW curves, `normal_vector` and `inward_normal` point in opposite directions.
Use `inward_normal` when you need orientation-independent results.

**Boundary handling (open curves):**
- Vertex `1`: 90° CCW rotation of `edge_tangent(c, 1)` (2D), or zero vector (N>2).
- Vertex `n`: same for the last edge.

### Examples

```julia
c = circle_curve(256)   # CCW unit circle

n5 = normal_vector(c, 5)
norm(n5)   # → 1.0

# For a CCW convex curve, normal_vector ≈ inward_normal
in5 = inward_normal(c, 5)
dot(n5, in5)   # → ≈ +1.0

# CW circle: normal_vector and inward_normal are anti-parallel
c_cw = orient(circle_curve(256); to=CW)
n5_cw  = normal_vector(c_cw, 5)
in5_cw = inward_normal(c_cw, 5)
dot(n5_cw, in5_cw)   # → ≈ -1.0

# Curvature vector and normal_vector are parallel for convex curves
kv = curvature_vector(c, 5)
dot(normalize(kv), normal_vector(c, 5))   # → ≈ +1.0
```

### See Also

[`inward_normal`](@ref), [`tangent_vector`](@ref), [`curvature_vector`](@ref)

---

## inward\_normal / outward\_normal

```julia
inward_normal(c, i)  → SVector{N,T}
outward_normal(c, i) → SVector{N,T}
```

Orientation-corrected unit normals at vertex `i` for **closed curves**.

`inward_normal` always points **toward the interior** of the enclosed region,
regardless of whether the curve is CCW or CW. `outward_normal` is its negation.

These functions throw `ArgumentError` if called on an open curve.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | **Closed** curve required. |
| `i` | `Int` | Vertex index. |

### Examples

```julia
# Unit circle (CCW by default)
c = circle_curve(256)

n_in  = inward_normal(c, 1)    # → points toward origin ≈ -vertex(c,1)
n_out = outward_normal(c, 1)   # → points away from origin ≈ +vertex(c,1)

# Verify: inward normal points toward interior for CCW circle
dot(n_in, -vertex(c, 1))   # → ≈ 1.0  (origin is inside)

# CW circle: same result (orientation-independent)
c_cw = orient(circle_curve(256); to=CW)
n_in_cw = inward_normal(c_cw, 1)
dot(n_in_cw, -vertex(c_cw, 1))   # → ≈ 1.0  (still points inward)

# Open curve: throws ArgumentError
c_open = OpenCurve([SVector(x, x^2) for x in range(-1.0, 1.0; length=64)])
inward_normal(c_open, 1)   # ERROR: ArgumentError

# Consistent with curvature vector (always co-directional)
kvs = curvature_vectors(c)
ins = inward_normals(c)
all(dot(kvs[i], ins[i]) ≥ -1e-10 for i in eachindex(ins))   # → true
```

### See Also

[`outward_normals`](@ref), [`normal_vector`](@ref), [`orient`](@ref)

---

## inward\_normals / outward\_normals

```julia
inward_normals(c)  → Vector{SVector{N,T}}
outward_normals(c) → Vector{SVector{N,T}}
```

All `nvertices(c)` inward (or outward) unit normals for a closed curve.
Orientation-independent analogues of [`normal_vectors`](@ref).

Both functions throw `ArgumentError` for open curves.

### Examples

```julia
c = fourier_curve(256; seed=2)

ins  = inward_normals(c)    # length nvertices(c)
outs = outward_normals(c)

# Verify: inward + outward = zero
all(norm(ins[i] + outs[i]) < 1e-12 for i in eachindex(ins))   # → true

# Offset curve (shrink inward by 0.1)
offset_pts = [vertex(c, i) + 0.1 * ins[i] for i in 1:nvertices(c)]
c_inner = ClosedCurve(offset_pts)

# Visualise
using GLMakie
curveplot(c; skin=:studio) + vectorplot(inward_normals(c); anchors=c,
                                        arrow_color=:coral, normalize=true)
```

### See Also

[`inward_normal`](@ref), [`outward_normal`](@ref), [`normal_vectors`](@ref)