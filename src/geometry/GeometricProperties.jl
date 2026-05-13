module GeometricProperties
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators, ..Lengths
    export is_convex_vertex, is_reflex_vertex, vertex_topology
    export VertexTopology, ConvexVertex, ReflexVertex
    export signed_area, centroid, isoperimetric_quotient, regular_isoperimetric_quotient, bounding_box

    """
        signed_area(c::AbstractDiscreteCurve{2,T})  →  T

    Signed area enclosed by a 2D curve using the shoelace formula.
    
    For closed curves: area of the enclosed region.
    For open curves: area enclosed by the curve and the line connecting endpoints.
    
    - Positive: counter-clockwise orientation
    - Negative: clockwise orientation
    - Zero: degenerate curve
    
    ```julia
    c = random_convex_polygon(8)
    A = signed_area(c)
    A > 0 ? println("CCW") : println("CW")
    ```
    """
    function signed_area(c::AbstractDiscreteCurve{2,T}) where T<:AbstractFloat
        n = nvertices(c)
        area = zero(T)
        
        if isclosed(c)
            @inbounds for i in 1:n
                p = c.points[i]
                q = c.points[mod1(i+1, n)]
                area += p[1] * q[2] - q[1] * p[2]
            end
        else
            @inbounds for i in 1:(n-1)
                p = c.points[i]
                q = c.points[i+1]
                area += p[1] * q[2] - q[1] * p[2]
            end
            p = c.points[n]
            q = c.points[1]
            area += p[1] * q[2] - q[1] * p[2]
        end
        
        area / 2
    end


    """
    centroid(c) → SVector{N,T}

    Vertex-averaged centroid. Works for both open and closed curves.
    """
    function centroid(c::AbstractDiscreteCurve{N,T}) where {N,T}
        sum(c.points) / nvertices(c)
    end

        """
        bounding_box(c) → (lo::SVector{N,T}, hi::SVector{N,T})

    Axis-aligned bounding box of all vertices.
    """
    function bounding_box(c::AbstractDiscreteCurve{N,T}) where {N,T}
        lo = c.points[1]; hi = c.points[1]
        @inbounds for i in 2:nvertices(c)
            p  = c.points[i]
            lo = min.(lo, p)
            hi = max.(hi, p)
        end
        lo, hi
    end
    function isoperimetric_quotient(c::AbstractDiscreteCurve{N,T}) where {N,T}
        A = abs(signed_area(c))
        L = arc_length(c)
        A > 0 && L > 0 ? (4 * π * A) / (L^2) : zero(T)
    end

    function regular_isoperimetric_quotient(n::Int)
       return π /n * cot(π/n)
    end
    
    regular_isoperimetric_quotient(c::AbstractDiscreteCurve{N,T}) where {N,T} = regular_isoperimetric_quotient(nvertices(c))

# ─────────────────────────────────────────────────────────────────────────────
    #  Vertex Topology Classification
    # ─────────────────────────────────────────────────────────────────────────────

    """
        VertexTopology

    Trait indicating whether a vertex is locally convex or reflex (concave).
    """
    abstract type VertexTopology end
    struct ConvexVertex <: VertexTopology end
    struct ReflexVertex <: VertexTopology end

    """
        vertex_topology(c::AbstractDiscreteCurve{2}, i::Int) → VertexTopology

    Classify vertex `i` as locally convex or reflex.

    **Definition (for 2D closed curves):**
    - **Convex vertex**: The angle at `i` turns in the direction consistent with the 
      curve's orientation. For CCW curves, this is a left turn; for CW curves, a right turn.
    - **Reflex (concave) vertex**: The angle turns opposite to the curve's orientation.

    **Implementation:** Checks the sign of the 2D cross product `(t_out - t_in) × t_out`:
    - Positive (CCW curve) or negative (CW curve) → ConvexVertex
    - Opposite sign → ReflexVertex

    # Example
    ```julia
    c = circle_curve(100)  # All vertices convex
    @assert all(vertex_topology(c, i) isa ConvexVertex for i in 1:nvertices(c))

    c_star = star_polygon(5)  # Mix of convex and reflex
    @assert any(vertex_topology(c_star, i) isa ReflexVertex for i in 1:nvertices(c_star))
    ```
    """
    function vertex_topology(c::AbstractDiscreteCurve{2,T}, i::Int)::VertexTopology where T
        isclosed(c) || throw(ArgumentError("vertex_topology requires a closed curve"))
        
        n = nvertices(c)
        
        t_in  = edge_tangent(c, mod1(i - 1, n))
        t_out = edge_tangent(c, i)
        
        # 2D cross product: t_in × t_out
        # Positive = left turn (CCW-consistent), Negative = right turn (CW-consistent)
        cross_prod = t_in[1] * t_out[2] - t_in[2] * t_out[1]
        
        curve_is_ccw = is_ccw(c)
        # For CCW curve: positive cross product (left turn) → convex
        # For CW curve: negative cross product (right turn) → convex
        is_convex = curve_is_ccw ? (cross_prod > zero(T)) : (cross_prod < zero(T))
        
        return is_convex ? ConvexVertex() : ReflexVertex()
    end

    """
        is_convex_vertex(c::AbstractDiscreteCurve{2}, i::Int) → Bool

    Check if vertex `i` is locally convex. Equivalent to `vertex_topology(c, i) isa ConvexVertex`.
    """
    function is_convex_vertex(c::AbstractDiscreteCurve{2}, i::Int)::Bool
        vertex_topology(c, i) isa ConvexVertex
    end

    """
        is_reflex_vertex(c::AbstractDiscreteCurve{2}, i::Int) → Bool

    Check if vertex `i` is locally reflex (concave). Equivalent to `vertex_topology(c, i) isa ReflexVertex`.
    """
    function is_reflex_vertex(c::AbstractDiscreteCurve{2}, i::Int)::Bool
        vertex_topology(c, i) isa ReflexVertex
    end
end
