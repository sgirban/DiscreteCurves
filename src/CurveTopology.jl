module CurveTopology

using StaticArrays
using ..AbstractTypes
using ..CurveTypes

export isclosed, isopen, nvertices, nedges, vertex, edge_segment

"""
    isclosed(c)  →  Bool
 
Returns `true` if the curve's last vertex connects back to its first.
"""
@inline isclosed(c::AbstractDiscreteCurve) = is_closed_topology(c)
 
"""
    isopen(c)  →  Bool
"""
@inline isopen(c::AbstractDiscreteCurve)   = is_open_topology(c)
 
"""
    nvertices(c)  →  Int
 
Number of distinct vertices.  For a closed polygon with `n` vertices, the
closing edge runs from vertex `n` to vertex `1`; it is NOT stored as an
extra point.
"""
@inline nvertices(c::AbstractDiscreteCurve) = length(c.points)
 
"""
    nedges(c)  →  Int
 
Number of edges:
  - Open curve:   n - 1
  - Closed curve: n        (the closing edge vertex[n]→vertex[1] is implicit)
"""
@inline function nedges(c::AbstractDiscreteCurve)
    n = nvertices(c)
    isclosed(c) ? n : n - 1
end

"""
    vertex(c, i)
 
Access the `i`-th vertex.  For closed curves, `i` is taken modulo `n`,
so `vertex(c, n+1) === vertex(c, 1)`.
 
```julia
c = ClosedCurve([SVector(1.0,0.0), SVector(0.0,1.0), SVector(-1.0,0.0)])
vertex(c, 4)  # same as vertex(c, 1)
```
"""
@inline vertex(c::AbstractDiscreteCurve, i::Int) =
    _vertex(topology_trait(c), c, i)
 
@inline function _vertex(::OpenTopology, c, i)
    @boundscheck checkbounds(c.points, i)
    @inbounds c.points[i]
end
 
@inline function _vertex(::ClosedTopology, c, i)
    # mod1: maps any integer to 1..n (1-based modular arithmetic)
    # This is a single CPU instruction when n is a power of 2.
    @inbounds c.points[mod1(i, nvertices(c))]
end
# ── Edge access ────────────────────────────────────────────────────────────────
#
"""
    edge_segment(c, i)  →  (p_i, p_{i+1})  as a 2-Tuple of SVectors
 
Returns the `i`-th edge as an ordered pair of endpoint positions.
For closed curves, edge `n` wraps from vertex `n` back to vertex `1`.
 
```julia
p, q = edge_segment(c, 3)
len  = norm(q - p)
```
"""
@inline function edge_segment(c::AbstractDiscreteCurve, i::Int)
    (vertex(c, i), vertex(c, i + 1))
end

# ── Base interface integration ─────────────────────────────────────────────────
#

Base.length(c::AbstractDiscreteCurve) = nvertices(c)
Base.size(c::AbstractDiscreteCurve) = (nvertices(c),)
Base.eltype(::Type{<:AbstractDiscreteCurve{N,T}}) where {N,T} = SVector{N,T}

@inline Base.getindex(c::AbstractDiscreteCurve, i::Int) = vertex(c, i)

Base.iterate(c::AbstractDiscreteCurve,i::Int=1) =
    i > nvertices(c) ? nothing : (vertex(c, i), i + 1)
Base.firstindex(c::AbstractDiscreteCurve) = 1
Base.lastindex(c::AbstractDiscreteCurve) = nvertices(c)


# ── Pretty printing ────────────────────────────────────────────────────────────
#
function Base.show(io::IO, c::AbstractDiscreteCurve{N,T}) where {N,T}
    topo = isclosed(c) ? "closed" : "open"
    print(io, "DiscreteCurve{$N,$T}($topo, $(nvertices(c)) vertices)")
end

function Base.show(io::IO, ::MIME"text/plain", c::AbstractDiscreteCurve{N,T}) where {N,T}
    topo   = isclosed(c) ? "closed" : "open"
    print(io, "DiscreteCurve{$N, $T}  —  $topo  —  $(nvertices(c)) vertices, $(nedges(c)) edges\n")
    n = min(3, nvertices(c))
    for i in 1:n
        print(io, "  [$i] $(c.points[i])\n")
    end
    nvertices(c) > n && print(io, "  … $(nvertices(c)-n) more vertices\n")
end

end