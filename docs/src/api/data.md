# Data Attachments

Mechanisms for associating scalar, vector, or arbitrary data with vertices and
edges of a curve. Three storage strategies are provided:

| Type | Storage | Use case |
|------|---------|---------|
| [`VertexData`](@ref) | Materialized `Vector` | Pre-computed fields; fast repeated access |
| [`EdgeData`](@ref) | Materialized `Vector` | Per-edge scalars (lengths, fluxes, …) |
| [`LazyVertexField`](@ref) | Zero storage; evaluated on demand | One-time access or very large curves |
| [`CurveDataBundle`](@ref) | Named dictionary of any of the above | Multi-channel scientific data |

All three are **topology-aware**: indexing wraps for closed curves.

---

## VertexData

```julia
VertexData{D, C}
```

A typed array of values of type `D`, one per vertex of a curve `C`. Indexing
wraps for closed curves (like `vertex(c, i)`).

### Construction

```julia
attach_vertex_data(c, data::Vector{D}) → VertexData{D,C}
attach_vertex_data(c, f)               → VertexData{D,C}
```

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | The curve the data belongs to. |
| `data` | `Vector{D}` | Pre-computed values. Length must equal `nvertices(c)`. |
| `f` | `Function` | Called as `f(vertex(c, i))` for each vertex. Return type inferred from first call. |

### Interface

| Operation | Syntax | Description |
|-----------|--------|-------------|
| Read at vertex `i` | `vd[i]` | Wraps for closed curves |
| Write at vertex `i` | `vd[i] = val` | Wraps for closed curves |
| Length | `length(vd)` | = `nvertices(c)` |
| Element type | `eltype(vd)` | = `D` |
| Iterate | `for x in vd` | |
| Collect | `collect(vd)` | Returns `copy(vd.data)` |
| Broadcast | `vd .* 2` | Works via `broadcastable` |

### Examples

```julia
c = fourier_curve(256; seed=3)

# ── From a pre-computed Vector ─────────────────────────────────────────────
κs = curvatures(c)
vd = attach_vertex_data(c, κs)
vd[1]       # → κs[1]
vd[0]       # → vd[nvertices(c)]  (wraps, closed curve)

# ── From a function: position → scalar ────────────────────────────────────
vd_r = attach_vertex_data(c, p -> norm(p))      # radial distance
vd_x = attach_vertex_data(c, p -> p[1])         # x-coordinate

# ── From a function: position → vector ────────────────────────────────────
vd_v = attach_vertex_data(c, p -> p / norm(p))  # D = SVector{2,Float64}

# ── Iteration ──────────────────────────────────────────────────────────────
for κ in vd
    println(κ)
end

# ── Broadcasting ──────────────────────────────────────────────────────────
vd_norm = vd ./ maximum(vd)      # returns a plain Vector{Float64}

# ── Attach curvature as VertexData ─────────────────────────────────────────
vd = attach_vertex_data(c, curvatures(c))
vertex_data(vd, 5)   # → curvature at vertex 5
```

### See Also

[`EdgeData`](@ref), [`CurveDataBundle`](@ref), [`LazyVertexField`](@ref)

---

## EdgeData

```julia
EdgeData{D, C}
```

A typed array of values of type `D`, one per **edge** of a curve `C`.
Length = `nedges(c)`. Indexing wraps for closed curves.

### Construction

```julia
attach_edge_data(c, data::Vector{D})      → EdgeData{D,C}
attach_edge_data(c, f)                    → EdgeData{D,C}
```

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | The curve. |
| `data` | `Vector{D}` | Pre-computed values. Length must equal `nedges(c)`. |
| `f` | `Function` | Called as `f(src, dst)` for each edge. Return type inferred from first call. |

### Examples

```julia
c = circle_curve(128)

# ── From a pre-computed Vector ─────────────────────────────────────────────
ℓ  = edge_lengths(c)
ed = attach_edge_data(c, ℓ)
ed[1]       # → first edge length
ed[0]       # → wraps for closed curves

# ── From a two-point function ──────────────────────────────────────────────
ed_len  = attach_edge_data(c, (p, q) -> norm(q - p))          # lengths
ed_mid  = attach_edge_data(c, (p, q) -> (p + q) / 2)          # midpoints
ed_dot  = attach_edge_data(c, (p, q) -> dot(p, q))            # dot product

# ── Interface (same as VertexData) ─────────────────────────────────────────
length(ed)       # → 128 = nedges(c)
eltype(ed)       # → Float64
collect(ed)      # → Vector{Float64}

for ℓᵢ in ed
    # iterate
end

# ── Use for edge-level statistics ──────────────────────────────────────────
ℓs = collect(attach_edge_data(c, (p,q) -> norm(q-p)))
println("Min edge: $(minimum(ℓs)), Max: $(maximum(ℓs)), Ratio: $(maximum(ℓs)/minimum(ℓs))")
```

### See Also

[`VertexData`](@ref), [`CurveDataBundle`](@ref)

---

## LazyVertexField

```julia
LazyVertexField{F, C}
lazy_vertex_field(c, f) → LazyVertexField
```

A zero-storage field that evaluates `f(vertex(c, i))` on demand. No values are
stored; each access recomputes from the vertex position.

Use when you need occasional access to a derived quantity and the memory of
a full `VertexData` is undesirable.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | The curve. |
| `f` | `Function` | `p → value` evaluated at each vertex position. |

### Interface

| Operation | Description |
|-----------|-------------|
| `lf[i]` | Evaluate at vertex `i` (wraps for closed curves) |
| `length(lf)` | = `nvertices(c)` |
| `collect(lf)` | Materialise into a `VertexData` |
| `for x in lf` | Iterate (recomputes each step) |

### Examples

```julia
c = fourier_curve(512; seed=7)

# Lazy radial distance
lf = lazy_vertex_field(c, p -> norm(p))
lf[5]          # evaluates norm(vertex(c, 5)) — no storage used
lf[0]          # wraps

# Materialise only when needed
vd = collect(lf)   # → VertexData{Float64}

# Use in a bundle as a lazy channel
b = bundle(c)
b[:radius] = p -> norm(p)     # stored as LazyVertexField automatically
b[:radius][5]                 # evaluates on access
```

### See Also

[`VertexData`](@ref), [`CurveDataBundle`](@ref)

---

## CurveDataBundle

```julia
CurveDataBundle{C}
bundle(c) → CurveDataBundle
```

A named container of multiple data channels (vertex, edge, or lazy) on a single
curve. Channels are keyed by `Symbol`. This is the recommended way to manage
scientific data associated with a curve (temperature, pressure, curvature,
signal, …) in a single object.

### Construction

```julia
b = bundle(c)
```

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | The curve the bundle is attached to. |

### Setting channels

```julia
b[:key] = data::Vector{D}       # auto-wrapped as VertexData or EdgeData by length
b[:key, toEdges=true] = data    # force EdgeData even when length ambiguous
b[:key] = vd::VertexData        # assign directly
b[:key] = ed::EdgeData          # assign directly
b[:key] = f::Function           # stored as LazyVertexField
```

Auto-wrapping rules for `Vector`:
- Length = `nvertices(c)` and `toEdges=false` → `VertexData`
- Length = `nedges(c)` → `EdgeData` (or `VertexData` for closed curves if lengths match)
- Use `toEdges=true` to force `EdgeData` when `nvertices == nedges` (closed curves)

### Reading channels

```julia
b[:key]         # → the full channel (VertexData, EdgeData, or LazyVertexField)
b[:key, i]      # → value at index i
b.key           # → dot-notation shortcut (same as b[:key])
vertex_data(b, i, :key)   # → b[:key, i]
edge_data(b, i, :key)     # → b[:key, i]
```

### Channel inspection

```julia
keys(b)               # → set of channel symbols
haskey(b, :key)       # → Bool
channel_type(b, :key) # → :vertex | :edge | :lazy | :unknown
length(b)             # → number of channels
```

### Examples

```julia
c = fourier_curve(256; seed=3)
b = bundle(c)

# ── Attach various data ────────────────────────────────────────────────────
b[:curvature]   = curvatures(c)                    # VertexData{Float64}
b[:arc_length]  = cumulative_length(c)             # VertexData{Float64}
b[:tangents]    = tangent_vectors(c)               # VertexData{SVector{2,Float64}}
b[:edge_len]    = edge_lengths(c)                  # EdgeData{Float64}
b[:radius]      = p -> norm(p)                     # LazyVertexField

# ── Read channels ──────────────────────────────────────────────────────────
κs = b[:curvature]         # the full VertexData
b[:curvature][5]           # value at vertex 5
b[:curvature, 5]           # same
b.curvature[5]             # dot-notation

# ── Inspect ────────────────────────────────────────────────────────────────
keys(b)                    # → Set([:curvature, :arc_length, :tangents, :edge_len, :radius])
channel_type(b, :curvature)  # → :vertex
channel_type(b, :edge_len)   # → :edge
channel_type(b, :radius)     # → :lazy

# ── Snapshot: all vertex channels at one vertex ────────────────────────────
s = snapshot(b, 5)
# → NamedTuple: (curvature=0.23, arc_length=0.48, radius=0.97, ...)
s.curvature     # field access

# ── Performance tip: extract channel before tight loop ─────────────────────
κ_ch = b[:curvature]        # one Dict lookup
for i in 1:nvertices(c)
    p = vertex(c, i)
    κ = κ_ch[i]             # plain array read — no Dict overhead
    process(p, κ)
end

# ── Iteration over all channels ────────────────────────────────────────────
for (name, channel) in b
    println(":$name  →  $(channel_type(b, name))")
end

# ── Co-iterate vertex channel with geometry ────────────────────────────────
for (p, κ) in zip(vertices(c), b[:curvature])
    println("pos=$p  κ=$κ")
end

# ── Display ────────────────────────────────────────────────────────────────
b
# CurveDataBundle  256v / 256e
#   :arc_length    [vertex, Float64, 256 entries]
#   :curvature     [vertex, Float64, 256 entries]
#   :edge_len      [edge, Float64, 256 entries]
#   :radius        [lazy, on demand, 256 entries]
#   :tangents      [vertex, SVector{2, Float64}, 256 entries]
```

---

## snapshot

```julia
snapshot(b, i) → NamedTuple
```

Collect all **vertex and lazy** channels of bundle `b` into a `NamedTuple` at
vertex `i`. `EdgeData` channels are skipped (they live on edges, not vertices).

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `b` | `CurveDataBundle` | A populated bundle. |
| `i` | `Int` | Vertex index. |

### Returns

`NamedTuple` with one field per vertex/lazy channel, each holding the value at
vertex `i`.

### Examples

```julia
c = circle_curve(64)
b = bundle(c)
b[:curvature]  = curvatures(c)
b[:arc_length] = cumulative_length(c)
b[:radius]     = p -> norm(p)

s = snapshot(b, 10)
# → (curvature = 1.006, arc_length = 0.940, radius = 1.000)
s.curvature     # → 1.006
s.arc_length    # → 0.940

# All vertices as a vector of snapshots
snaps = [snapshot(b, i) for i in 1:nvertices(c)]
```

### See Also

[`CurveDataBundle`](@ref), [`vertex_data`](@ref)

---

## channel\_type

```julia
channel_type(b, key) → Symbol
```

Returns the type tag of channel `key` in bundle `b`:

| Return value | Meaning |
|---|---|
| `:vertex` | `VertexData` — materialised, one value per vertex |
| `:edge` | `EdgeData` — materialised, one value per edge |
| `:lazy` | `LazyVertexField` — computed on access |
| `:unknown` | Some other type was stored directly |

### Examples

```julia
c = circle_curve(64)
b = bundle(c)
b[:kappa]    = curvatures(c)
b[:el]       = edge_lengths(c)
b[:lazy_r]   = p -> norm(p)

channel_type(b, :kappa)    # → :vertex
channel_type(b, :el)       # → :edge
channel_type(b, :lazy_r)   # → :lazy
```

### See Also

[`CurveDataBundle`](@ref), [`snapshot`](@ref)