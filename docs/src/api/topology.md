# Topology & Vertex Access

These functions query the structure of a curve: how many vertices and edges it
has, whether it is open or closed, and how to retrieve individual vertices and
edges. They are the foundation that all other computations build on.

---

## isclosed / isopen

```julia
isclosed(c) → Bool
isopen(c)   → Bool
```

Query the topological type of a curve.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | Any discrete curve. |

### Returns

`Bool`.

### Details

A **closed** curve has an implicit edge from its last vertex back to its first.
This means:
- `nedges(c) == nvertices(c)` for closed curves.
- `nedges(c) == nvertices(c) - 1` for open curves.
- `vertex(c, nvertices(c)+1)` wraps to `vertex(c, 1)` for closed curves.

`isclosed` and `isopen` are strict complements: exactly one is always `true`.

### Examples

```julia
c_open   = OpenCurve([SVector(0.0,0.0), SVector(1.0,0.0), SVector(2.0,1.0)])
c_closed = ClosedCurve([SVector(1.0,0.0), SVector(0.0,1.0), SVector(-1.0,0.0)])

isclosed(c_open)    # → false
isopen(c_open)      # → true
isclosed(c_closed)  # → true
isopen(c_closed)    # → false

# Use in dispatch or conditional logic
if isclosed(c)
    println("Closed curve: $(nedges(c)) edges")
else
    println("Open curve: endpoints at $(vertex(c,1)) and $(vertex(c,nvertices(c)))")
end
```

### See Also

[`nvertices`](@ref), [`nedges`](@ref), [`DiscreteCurve`](@ref)

---

## nvertices

```julia
nvertices(c) → Int
```

Number of distinct stored vertices. For a closed curve with `n` vertices, the
closing edge runs from vertex `n` back to vertex `1`; it is **not** represented
as an extra stored point.

### Examples

```julia
c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)
nvertices(c)   # → 128

# Relation to Base.length
length(c)      # → 128  (DiscreteCurve implements Base.length = nvertices)
```

### See Also

[`nedges`](@ref), [`vertex`](@ref)

---

## nedges

```julia
nedges(c) → Int
```

Number of edges.

| Topology | Value |
|----------|-------|
| Closed   | `nvertices(c)` |
| Open     | `nvertices(c) - 1` |

### Examples

```julia
c_open   = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, π, 64)
c_closed = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 64; closed=true)

nedges(c_open)    # → 63
nedges(c_closed)  # → 64

# Total arc length via edge count
sum(norm(vertex(c,i+1) - vertex(c,i)) for i in 1:nedges(c))  # == arc_length(c)
```

### See Also

[`nvertices`](@ref), [`edge_segment`](@ref), [`edges`](@ref)

---

## vertex

```julia
vertex(c, i) → SVector{N,T}
```

Retrieve the `i`-th vertex position.

- **Open curves**: `i` must be in `1:nvertices(c)`. Out-of-bounds access is
  checked in debug mode.
- **Closed curves**: `i` is taken modulo `nvertices(c)` (1-based), so
  `vertex(c, 0)`, `vertex(c, nvertices(c)+1)`, etc. all wrap correctly.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | The curve. |
| `i` | `Int` | 1-based vertex index. Wraps for closed curves. |

### Returns

`SVector{N,T}` — the position of vertex `i`.

### Examples

```julia
c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 8; closed=true)

vertex(c, 1)    # → SVector(1.0, 0.0)  (first vertex)
vertex(c, 8)    # → last vertex
vertex(c, 9)    # → same as vertex(c, 1) — wraps for closed curves
vertex(c, 0)    # → same as vertex(c, 8) — wraps backwards

# c[i] is identical to vertex(c, i)
c[3]            # → vertex(c, 3)

# Loop over all vertices (direct indexing)
for i in 1:nvertices(c)
    p = vertex(c, i)
    println("v$i = $p")
end
```

### Performance

`vertex` is `@inline` and has zero allocation. For closed curves, the modular
arithmetic is a single `mod1` call (one CPU instruction when `n` is a power of 2).
Prefer `vertex(c, i)` inside tight loops over `collect`-ing all vertices first.

### See Also

[`edge_segment`](@ref), [`vertices`](@ref), [`nvertices`](@ref)

---

## edge\_segment

```julia
edge_segment(c, i) → (SVector{N,T}, SVector{N,T})
```

Return the `i`-th edge as an ordered pair of endpoint positions
`(vertex(c, i), vertex(c, i+1))`.

For closed curves, edge `nedges(c)` wraps from the last vertex back to the
first.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | The curve. |
| `i` | `Int` | 1-based edge index in `1:nedges(c)`. |

### Returns

`Tuple{SVector{N,T}, SVector{N,T}}` — `(src, dst)` endpoint positions.

### Examples

```julia
c = ClosedCurve([SVector(1.0,0.0), SVector(0.0,1.0), SVector(-1.0,0.0)])

p, q = edge_segment(c, 1)   # → (SVector(1,0), SVector(0,1))
len  = norm(q - p)           # edge length

p, q = edge_segment(c, 3)   # last edge: wraps from vertex 3 → vertex 1
norm(q - p)

# Total arc length
L = sum(norm(q-p) for (p,q) in (edge_segment(c,i) for i in 1:nedges(c)))
# Equivalent to arc_length(c)
```

### See Also

[`vertex`](@ref), [`edges`](@ref), [`edge_lengths`](@ref)

---

## vertices

```julia
vertices(c) → VertexIterator
```

Lazy iterator over all vertices of `c`. Yields `SVector{N,T}` values in order
from vertex 1 to vertex `nvertices(c)`. Zero allocations — nothing is collected
until you explicitly call `collect`.

### Examples

```julia
c = circle_curve(128)

# Sum all vertex positions (centroid × n)
total = sum(p for p in vertices(c))

# Find all vertices in the upper half-plane
upper = filter(p -> p[2] > 0, vertices(c))

# Co-iterate with curvature data
κs = curvatures(c)
for (p, κ) in zip(vertices(c), κs)
    println("position=$(round.(p; digits=3))  κ=$κ")
end

# Collect into a plain Vector{SVector{2,Float64}} if needed
pts = collect(vertices(c))

# Check a property
any(p -> norm(p) > 1.5, vertices(c))
all(p -> p[1]^2 + p[2]^2 < 1.01, vertices(c))
```

### See Also

[`vertex`](@ref), [`edges`](@ref), [`@where`](@ref)

---

## edges

```julia
edges(c) → EdgeIterator
```

Lazy iterator over all edges of `c`. Each element is a **named tuple**
`(; src::SVector{N,T}, dst::SVector{N,T})`, representing one directed edge.

Yields `nedges(c)` elements. The closing edge of a closed curve (from vertex
`n` to vertex `1`) is included automatically.

### Examples

```julia
c = circle_curve(128)

# Arc length — zero allocations
L = sum(norm(dst - src) for (; src, dst) in edges(c))

# Edge midpoints
mids = [SVector((src+dst)/2) for (; src, dst) in edges(c)]

# Mean edge length
mean_ℓ = sum(norm(dst - src) for (; src, dst) in edges(c)) / nedges(c)

# Collect all edges
edge_list = collect(edges(c))
edge_list[1].src   # → first vertex
edge_list[1].dst   # → second vertex
```

### See Also

[`edge_segment`](@ref), [`edge_lengths`](@ref), [`vertices`](@ref)

---

## edge\_windows

```julia
edge_windows(c; k = 3) → WindowIterator
```

Lazy iterator yielding `k`-point sliding windows centred at each interior vertex.
Each element is an `NTuple{k, SVector{N,T}}`.

- **Closed curves**: yields `nvertices(c)` windows (wraps at both ends).
- **Open curves**: yields `nvertices(c) - 2⌊k/2⌋` windows (boundary vertices
  are excluded because a centred window of half-width `k÷2` would reach
  out-of-bounds).

`k` must be **odd** (centred windows only) and ≥ 3.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | The curve. |
| `k` | `Int` | Window size. Must be odd and ≥ 3. Default: `3`. |

### Returns

`WindowIterator{k,C}` — lazy, zero-allocation.

### Examples

```julia
c = fourier_curve(128)

# 3-point windows (default) — basis of all curvature computation
for (prev, curr, next) in edge_windows(c)
    # prev = vertex(c, i-1), curr = vertex(c, i), next = vertex(c, i+1)
    t1 = curr - prev
    t2 = next - curr
    println("turning angle ≈ ", acos(clamp(dot(t1,t2)/(norm(t1)*norm(t2)), -1, 1)))
end

# 5-point windows for wider stencils
for (p1, p2, p3, p4, p5) in edge_windows(c; k=5)
    # p3 is the centre vertex
end

# Count windows (open vs closed)
length(edge_windows(c))         # nvertices (closed)
length(edge_windows(c; k=5))    # nvertices - 4 (closed, half=2)
```

### See Also

[`curvature`](@ref), [`turning_angle`](@ref), [`vertices`](@ref)