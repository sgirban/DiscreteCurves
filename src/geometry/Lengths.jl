module Lengths
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators

    export arc_length, edge_lengths, cumulative_length
@inline _edge_length(src, dst) = norm(dst - src)

"""
    arc_length(c)  →  T
 
Total arc length of the curve: sum of all edge lengths.
 
Complexity: O(n).
 
```julia
L = arc_length(c)
```
"""
function arc_length(c::AbstractDiscreteCurve{N,T}) where {N,T}
    s = zero(T)
    for (; src, dst) in edges(c)
        s += _edge_length(src, dst)
    end
    s
end

"""
    edge_lengths(c)  →  Vector{T}
 
Returns a vector of the `nedges(c)` edge lengths.
 
```julia
ℓ = edge_lengths(c)   # ℓ[i] = ||vertex(c,i+1) - vertex(c,i)||
```
"""
function edge_lengths(c::AbstractDiscreteCurve{N,T}) where {N,T}
    n = nedges(c)
    l = Vector{T}(undef, n)
    @inbounds for (i, (; src, dst)) in enumerate(edges(c))
        l[i] = _edge_length(src, dst)
    end
    l
end

"""
    cumulative_length(c)  →  Vector{T}
 
Arc-length parameter `s[i]` at each vertex.
`s[1] == 0`, `s[i+1] = s[i] + ||p_{i+1} - p_i||`.
 
Useful for arc-length parametrisation and reparametrisation.
 
```julia
s = cumulative_length(c)
# s is length nvertices(c): s[1]=0, s[end]=arc_length(c)
```
"""
function cumulative_length(c::AbstractDiscreteCurve{N,T}) where {N,T}
    n  = nvertices(c)
    s  = Vector{T}(undef, n)
    s[1] = zero(T)
    @inbounds for (i, (; src, dst)) in enumerate(edges(c))
        s[i+1] = s[i] + _edge_length(src, dst)
    end
    s
end
end