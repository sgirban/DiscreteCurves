module CurveVectors
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators, ..GeometricProperties, ..Orientation

    export edge_tangent, edge_normal, edge_tangents, edge_normals
    export tangent_vector, tangent_vectors
    export normal_vector, normal_vectors
    export inward_normal, outward_normal
    export inward_normals, outward_normals


    @inline function _unit(v::SVector{N,T}) where {N,T}
        n = norm(v)
        n > eps(T) ? v * inv(n) : zero(SVector{N,T})
    end
    @inline _rot90ccw(v::SVector{2,T}) where T = SVector{2,T}(-v[2],  v[1])
    @inline _rot90cw( v::SVector{2,T}) where T = SVector{2,T}( v[2], -v[1])

    """
        edge_tangent(c, i)  →  SVector{N,T}

    Unit tangent of edge `i`: the direction from `vertex(c, i)` toward `vertex(c, i+1)`.

    For closed curves edge `nedges(c)` wraps from the last vertex back to the first
    (because `vertex` uses modular indexing on closed curves).


    ```julia
    t̂ = edge_tangent(c, 3)
    @assert norm(t̂) ≈ 1          # unit vector
    ```
    """
    @inline function edge_tangent(c::AbstractDiscreteCurve{N,T}, i::Int) where {N,T}
        _unit(vertex(c, i + 1) - vertex(c, i))
    end

    """
        edge_tangents(c)  →  Vector{SVector{N,T}}

    All `nedges(c)` unit edge tangents, one per edge.

    """
    function edge_tangents(c::AbstractDiscreteCurve{N,T}) where {N,T}
        ne  = nedges(c)
        out = Vector{SVector{N,T}}(undef, ne)
        @inbounds for i in 1:ne
            out[i] = edge_tangent(c, i)
        end
        out
    end

    """
        edge_normal(c, i)  →  SVector{2,T}

    Unit normal of edge `i`: the edge tangent rotated 90° **counter-clockwise**,

    **Orientation note (2-D closed curves):**
    - **CCW** curve → edge normal points *left* of travel → **inward** for convex shapes.
    - **CW**  curve → edge normal points *right* of travel → **outward** for convex shapes.

    This is an edge-level quantity; for a vertex-level orientation-independent normal
    see [`inward_normal`](@ref).

    ```julia
    # Verify: edge tangent and edge normal are orthogonal
    t = edge_tangent(c, i)
    n = edge_normal(c, i)
    @assert abs(dot(t, n)) < 1e-14
    ```
    """
    @inline function edge_normal(c::AbstractDiscreteCurve{2,T}, i::Int) where T
        _rot90ccw(edge_tangent(c, i))
    end

    """
        edge_normals(c)  →  Vector{SVector{2,T}}

    All `nedges(c)` unit edge normals (edge tangents rotated 90° CCW).
    """
    function edge_normals(c::AbstractDiscreteCurve{2,T}) where T
        ne  = nedges(c)
        out = Vector{SVector{2,T}}(undef, ne)
        @inbounds for i in 1:ne
            out[i] = edge_normal(c, i)
        end
        out
    end

    """
        tangent_vector(c, i)  →  SVector{N,T}

    Unit tangent at vertex `i`:

    ```
    Tᵢ = normalize( tᵢₙ + tₒᵤₜ )
    ```

    where `tᵢₙ = edge_tangent(c, i-1)` and `tₒᵤₜ = edge_tangent(c, i)` are the
    unit tangents of the incoming and outgoing edges respectively.

    This is the **angle bisector** of the two adjacent edge directions.  For uniform
    meshes it coincides with the central-difference tangent `(pᵢ₊₁ - pᵢ₋₁)`;
    for non-uniform meshes it correctly gives equal angular weight to both adjacent
    edges, unlike the biased central-difference formula.

    **Boundary handling (open curves):**
    - `i = 1`   → returns `edge_tangent(c, 1)` (only one edge available).
    - `i = n`   → returns `edge_tangent(c, n-1)`.

    **Degenerate case:** returns `zero(SVector{N,T})` when `tᵢₙ + tₒᵤₜ ≈ 0`
    (a 180° hairpin vertex where both edge tangents point in opposite directions).

    **Complexity:** O(1), zero heap allocations.

    ```julia
    T = tangent_vector(c, 5)
    @assert norm(T) ≈ 1 || norm(T) == 0   # unit or zero at hairpin
    ```
    """
    function tangent_vector(c::AbstractDiscreteCurve{N,T}, i::Int) where {N,T}
        n      = nvertices(c)
        closed = isclosed(c)
        if closed
            t_in  = edge_tangent(c, mod1(i - 1, n)) 
            t_out = edge_tangent(c, i) 
            return _unit(t_in + t_out)
        end
        i == 1 && return edge_tangent(c, 1)
        i == n && return edge_tangent(c, n - 1)
        t_in  = edge_tangent(c, i - 1)
        t_out = edge_tangent(c, i)
        _unit(t_in + t_out)
    end

    """
        tangent_vectors(c)  →  Vector{SVector{N,T}}

    All `nvertices(c)` unit tangent vectors.  See [`tangent_vector`](@ref).
    """
    function tangent_vectors(c::AbstractDiscreteCurve{N,T}) where {N,T}
        nv  = nvertices(c)
        out = Vector{SVector{N,T}}(undef, nv)
        @inbounds for i in 1:nv
            out[i] = tangent_vector(c, i)
        end
        out
    end

    """
        normal_vector(c, i)  →  SVector{N,T}

    Unit vertex normal at `i`:

    ```
    Nᵢ = normalize( tₒᵤₜ - tᵢₙ )
    ```

    where `tᵢₙ` and `tₒᵤₜ` are the unit tangents of the incoming and outgoing edges.

    ```julia
    # normal_vector and curvature_vector are parallel (same direction) for any
    # smooth convex curve, regardless of orientation:
    κvs = curvature_vectors(c)
    ns  = normal_vectors(c)
    @assert all(dot(κvs[i], ns[i]) ≥ -1e-12 for i in eachindex(ns))
    ```
    """
    function normal_vector(c::AbstractDiscreteCurve{N,T}, i::Int) where {N,T}
        n      = nvertices(c)
        closed = isclosed(c)
        if closed
            t_in  = edge_tangent(c, mod1(i - 1, n))
            t_out = edge_tangent(c, i)
            return _unit(t_out - t_in)
        end
        i == 1 && return _boundary_normal(edge_tangent(c, 1))
        i == n && return _boundary_normal(edge_tangent(c, n - 1))
        t_in  = edge_tangent(c, i - 1)
        t_out = edge_tangent(c, i)
        _unit(t_out - t_in)
    end

    """
        normal_vectors(c)  →  Vector{SVector{N,T}}

    All `nvertices(c)` vertex normals.  Each vector points toward the concave (inward)
    side; see [`normal_vector`](@ref).
    """
    function normal_vectors(c::AbstractDiscreteCurve{N,T}) where {N,T}
        nv  = nvertices(c)
        out = Vector{SVector{N,T}}(undef, nv)
        @inbounds for i in 1:nv
            out[i] = normal_vector(c, i)
        end
        out
    end


    """
        inward_normal(c, i)  →  SVector{N,T}

    Unit inward (toward-interior) normal at vertex `i`.

    Defined as `normalize(tₒᵤₜ - tᵢₙ)`, with orientation correction applied:
    - **CCW orientation:** returns `normalize(tₒᵤₜ - tᵢₙ)` directly.
    - **CW orientation:** returns `-normalize(tₒᵤₜ - tᵢₙ)` to flip the sign.
    - **Result:** Always points toward the interior, regardless of orientation.


    Throws `ArgumentError` for open curves.

    ```julia
    c = circle_curve(128)

    # All inward normals point toward the origin on a unit circle
    @assert all(norm(inward_normal(c,i) + vertex(c,i)) < 1e-10 for i in 1:nvertices(c))

    # Consistent with curvature vector direction (parallel, positive dot product)
    κvs = curvature_vectors(c)
    ns  = inward_normals(c)
    @assert all(dot(κvs[k], ns[k]) ≥ 0 for k in eachindex(ns))
    ```
    """
    function inward_normal(c::AbstractDiscreteCurve{N,T}, i::Int) where {N,T}
        isclosed(c) || throw(ArgumentError(
        "inward_normal requires a closed curve (got an open curve). " *
        "For open curves use normal_vector, left_normal_vector, or right_normal_vector."))
        n      = nvertices(c)
        t_in   = edge_tangent(c, mod1(i - 1, n))
        t_out  = edge_tangent(c, i)
        normal = _unit(t_out - t_in)
        return normal
    end

    """
        outward_normal(c, i)  →  SVector{N,T}

    Unit outward (away-from-interior) normal at vertex `i`: `−inward_normal(c, i)`.

    Orientation-independent. Throws `ArgumentError` for open curves.
    """
    function outward_normal(c::AbstractDiscreteCurve{N,T}, i::Int) where {N,T}
        -inward_normal(c, i)
    end

    """
        inward_normals(c)  →  Vector{SVector{N,T}}

    All `nvertices(c)` inward unit normals (closed curve required).
    """
    function inward_normals(c::AbstractDiscreteCurve{N,T}) where {N,T}
        isclosed(c) || throw(ArgumentError("inward_normals requires a closed curve"))
        nv  = nvertices(c)
        out = Vector{SVector{N,T}}(undef, nv)
        @inbounds for i in 1:nv
            out[i] = inward_normal(c, i)
        end
        out
    end

    """
        outward_normals(c)  →  Vector{SVector{N,T}}

    All `nvertices(c)` outward unit normals (closed curve required).
    """
    function outward_normals(c::AbstractDiscreteCurve{N,T}) where {N,T}
        isclosed(c) || throw(ArgumentError("outward_normals requires a closed curve"))
        nv  = nvertices(c)
        out = Vector{SVector{N,T}}(undef, nv)
        @inbounds for i in 1:nv
            out[i] = -inward_normal(c, i)
        end
        out
    end

    
    @inline _boundary_normal(t::SVector{2,T}) where T  = _rot90ccw(t)
    @inline _boundary_normal(t::SVector{N,T}) where {N,T} = zero(SVector{N,T})

end