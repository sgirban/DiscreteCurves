module Iterators
    using StaticArrays
    using ..AbstractTypes
    using ..CurveTypes
    using ..CurveTopology

    export vertices, edges, edge_windows
    struct VertexIterator{C<:AbstractDiscreteCurve}
        curve::C
    end

    """
    vertices(c)
 
Lazy iterator over all vertices of `c`.
 
```julia
for p in vertices(c)
    println(norm(p))
end
```
"""
@inline vertices(c::AbstractDiscreteCurve) = VertexIterator(c)

Base.length(it::VertexIterator) = nvertices(it.curve)
Base.eltype(it::VertexIterator{C}) where {C} = eltype(it.curve)

Base.iterate(it::VertexIterator, state=1) = state > nvertices(it.curve) ? nothing : (vertex(it.curve, state), state + 1)

# ── EdgeIterator ───────────────────────────────────────────────────────────────
#
struct EdgeIterator{C<:AbstractDiscreteCurve}
    curve::C
end

"""
    edges(c)
 
Lazy iterator over all edges `(src, dst)` of `c`.
 
```julia
# Sum of edge lengths (arc length) – zero allocations
arc_len = sum(norm(dst - src) for (; src, dst) in edges(c))
 
# Collect into a vector when needed
edge_list = collect(edges(c))
```
"""
@inline edges(c::AbstractDiscreteCurve) = EdgeIterator(c)

Base.length(it::EdgeIterator) = nedges(it.curve)
Base.eltype(::EdgeIterator{C}) where {C<:AbstractDiscreteCurve{N,T}} where {N,T} = 
    @NamedTuple{src::Svector{N,T}, dst::SVector{N,T}}

function Base.iterate(it::EdgeIterator, state=1) 
    state > nedges(it.curve) && return nothing
    p = vertex(it.curve, state)
    q = vertex(it.curve, state + 1)
    (; src=p, dst=q), state + 1   
end

# ── WindowIterator ─────────────────────────────────────────────────────────────
#
struct WindowIterator{k,C<:AbstractDiscreteCurve}
    curve::C
    function WindowIterator{k}(c::C) where {k,C}
        k ≥ 2 || throw(ArgumentError("Window size k must be ≥ 2, got $k"))
        isodd(k) || throw(ArgumentError("Winow must be cenetred, thus k must be odd, got $k"))
        new{k,C}(c)
    end
        
end
"""
    edge_windows(c; k=3)
 
Lazy iterator yielding `k`-point sliding windows `(p_{i-1}, p_i, p_{i+1})`.
For closed curves the window wraps.  For open curves, interior vertices only
(endpoints are excluded because the centred window would reach out of bounds).
 
```julia
for (prev, curr, next) in edge_windows(c; k=3)
    θ = turning_angle_at(prev, curr, next)
end
```
"""
@inline function edge_windows(c::AbstractDiscreteCurve; k::Int=3)
    WindowIterator{k}(c)
end

@inline _half(k) = k ÷ 2
 
function Base.length(it::WindowIterator{k}) where k
    c = it.curve
    isclosed(c) ? nvertices(c) : nvertices(c) - 2*_half(k)
end
Base.eltype(::WindowIterator{k, C}) where {k, C<:AbstractDiscreteCurve{N,T}} where {N,T} =
    NTuple{k, SVector{N,T}}
 
function Base.iterate(it::WindowIterator{k}, state=nothing) where k
    c    = it.curve
    half = _half(k)
    n    = nvertices(c)
 
    i_start = isclosed(c) ? 1       : half+1
    i_end   = isclosed(c) ? n       : n-half
    i       = state === nothing ? i_start : state
 
    i > i_end && return nothing
 
    win = ntuple(j -> vertex(c, i - half + j - 1), Val(k))
    win, i+1
end
end