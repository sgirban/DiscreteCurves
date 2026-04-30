# Constructors

Discrete curves in **DiscreteCurves.jl** are built from an ordered sequence of
points in ``\mathbb{R}^N``. Every constructor returns a
[`DiscreteCurve{N,T}`](@ref) — a type-stable, immutable polyline parameterised by
ambient dimension `N` and scalar type `T <: AbstractFloat`.

---

## DiscreteCurve

```julia
DiscreteCurve(points; closed = false)
DiscreteCurve(points, closed::Bool)
DiscreteCurve(M::AbstractMatrix{T}; closed = false)
DiscreteCurve(f, t0, t1, n; closed = false, T = Float64)
```

The primary concrete curve type. Construct from a list of points, a matrix whose
columns are points, or a parametric function.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `points` | `Vector`, `Vector{SVector}`, `Vector{Tuple}` | Ordered sequence of vertex positions. Tuples and plain `Vector`s are automatically converted to `SVector`. |
| `closed` | `Bool` | If `true`, an implicit edge connects the last vertex back to the first. Default: `false`. |
| `M` | `AbstractMatrix{T}` | Matrix of size `N × n`; column `i` becomes vertex `i`. |
| `f` | `Function` | Parametric map `t → SVector{N,T}` (or Tuple / Vector). |
| `t0`, `t1` | `Real` | Parameter range. |
| `n` | `Int` | Number of sample points (vertices). Must be ≥ 2. |
| `T` | `Type{<:AbstractFloat}` | Scalar precision. Default `Float64`. Use `Float32` for GPU work. |

### Returns

`DiscreteCurve{N,T}` — the assembled polyline.

### Examples

```julia
using DiscreteCurves
using StaticArrays

# ── From a list of SVectors ────────────────────────────────────────────────
pts = [SVector(cos(θ), sin(θ)) for θ in range(0, 2π; length=9)[1:8]]
c   = DiscreteCurve(pts; closed=true)

# ── From plain tuples (auto-converted) ────────────────────────────────────
c = DiscreteCurve([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)]; closed=true)

# ── From a matrix (columns = points) ──────────────────────────────────────
M = rand(2, 64)          # 2D curve with 64 vertices
c = DiscreteCurve(M; closed=false)

# ── 3D curve ──────────────────────────────────────────────────────────────
pts3 = [SVector(cos(t), sin(t), t/10) for t in range(0, 4π; length=128)]
c3   = DiscreteCurve(pts3)          # open 3D helix

# ── From a parametric function (alias for ParametricCurve) ────────────────
c = DiscreteCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)

# ── Float32 precision (useful for GPU pipelines) ───────────────────────────
c32 = DiscreteCurve(Float32[1 0; 0 1; -1 0; 0 -1]'; closed=true)
```

### Notes

- For **closed** curves, `nvertices == n` and `nedges == n`; the closing edge is
  implicit (no duplicate first/last point is stored).
- For **open** curves, `nvertices == n` and `nedges == n - 1`.
- The type parameter `N` is inferred at construction time; it cannot change
  after construction.

### See Also

[`ClosedCurve`](@ref), [`OpenCurve`](@ref), [`ParametricCurve`](@ref), [`@curve`](@ref)

---

## ClosedCurve

```julia
ClosedCurve(pts)
```

Shorthand for `DiscreteCurve(pts; closed=true)`. The last vertex implicitly
connects back to the first — do **not** duplicate the first point.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `pts` | Any point collection | Same as `DiscreteCurve`. Tuples, `SVector`s, and plain `Vector`s all accepted. |

### Returns

`DiscreteCurve{N,T}` with `closed = true`.

### Examples

```julia
# Unit square (4 vertices, 4 edges)
c = ClosedCurve([SVector(1,0), SVector(0,1), SVector(-1,0), SVector(0,-1)])

# Circle from tuples
pts = [(cos(θ), sin(θ)) for θ in range(0, 2π; length=129)[1:128]]
c   = ClosedCurve(pts)

# Verify topology
isclosed(c)    # → true
nvertices(c)   # → 128
nedges(c)      # → 128
```

### See Also

[`OpenCurve`](@ref), [`DiscreteCurve`](@ref)

---

## OpenCurve

```julia
OpenCurve(pts)
```

Shorthand for `DiscreteCurve(pts; closed=false)`. Endpoints are free.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `pts` | Any point collection | Ordered vertex positions. |

### Returns

`DiscreteCurve{N,T}` with `closed = false`.

### Examples

```julia
# Parabolic arc
xs = range(-2.0, 2.0; length=64)
c  = OpenCurve([SVector(x, x^2) for x in xs])

# 3D open helix
c3 = OpenCurve([SVector(cos(t), sin(t), t/(2π)) for t in range(0, 4π; length=256)])

nvertices(c3)   # → 256
nedges(c3)      # → 255  (open: nedges = nvertices - 1)
isopen(c3)      # → true
```

### See Also

[`ClosedCurve`](@ref), [`DiscreteCurve`](@ref)

---

## ParametricCurve

```julia
ParametricCurve(f, t0, t1, n; closed = false, T = Float64)
ParametricCurve(fs::NTuple{N,Function}, t0, t1, n; closed = false, T = Float64)
ParametricCurve(f, t0, t1; step, closed = false, T = Float64)
```

Sample a parametric map ``f: \mathbb{R} \to \mathbb{R}^N`` at `n` equally spaced
parameter values and assemble the result into a `DiscreteCurve`.

For **closed** curves the parameter range `[t0, t1)` is sampled (the endpoint
`t1` is excluded because `f(t0) == f(t1)` by assumption); for **open** curves
the full closed interval `[t0, t1]` is sampled.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `f` | `Function` | Parametric map. May return `SVector{N,T}`, `Tuple`, or `Vector`. |
| `fs` | `NTuple{N,Function}` | Tuple of component functions `(fx, fy, ...)`. Each is called as `fs[d](t)`. |
| `t0` | `Real` | Start of parameter range. |
| `t1` | `Real` | End of parameter range. |
| `n` | `Int` | Number of sample points. Must be ≥ 2. |
| `step` | `Real` | Alternative to `n`: use a fixed parameter step size. |
| `closed` | `Bool` | Topology of the result. Default: `false`. |
| `T` | `Type{<:AbstractFloat}` | Scalar precision. Default: `Float64`. |

### Returns

`DiscreteCurve{N,T}` sampled from `f`.

### Examples

```julia
# ── Basic closed curves ────────────────────────────────────────────────────
circle  = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)
ellipse = ParametricCurve(t -> SVector(2cos(t), sin(t)), 0, 2π, 256; closed=true)

# ── 3D curves ─────────────────────────────────────────────────────────────
helix = ParametricCurve(t -> SVector(cos(t), sin(t), t/(2π)), 0, 4π, 256)
torus_knot = ParametricCurve(
    t -> SVector((2 + cos(3t))*cos(2t),
                 (2 + cos(3t))*sin(2t),
                 sin(3t)),
    0, 2π, 512; closed=true)

# ── Component-function tuple syntax ───────────────────────────────────────
lissajous = ParametricCurve((sin, cos), 0, 2π, 256; closed=true)
#  equivalent to: t -> SVector(sin(t), cos(t))

# ── Using a fixed step size instead of n ──────────────────────────────────
c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π; step=0.05, closed=true)

# ── Float32 precision ──────────────────────────────────────────────────────
c32 = ParametricCurve(t -> SVector(cos(t), sin(t)), 0f0, 2f0*π, 128;
                      closed=true, T=Float32)
eltype(c32)   # → SVector{2,Float32}

# ── Lemniscate (figure eight) ──────────────────────────────────────────────
lemn = ParametricCurve(
    t -> SVector(cos(t)/(1+sin(t)^2), sin(t)*cos(t)/(1+sin(t)^2)),
    0, 2π, 512; closed=true)
```

### Performance tip

Returning a statically-sized `SVector{N,T}` from `f` is the fastest path. The
library infers `N` and `T` from the first evaluation of `f` and pre-allocates
accordingly. Avoid returning different types on different calls.

### See Also

[`@curve`](@ref), [`DiscreteCurve`](@ref), [`circle_curve`](@ref)