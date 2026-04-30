module GeometricProperties
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators, ..Lengths

    export signed_area, centroid, isoperimetric_quotient, regular_isoperimetric_quotient

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
        centroid(c::AbstractDiscreteCurve{N,T})  →  SVector{N,T}

    Compute the centroid (center of mass) of a closed curve.
    
    For a discrete curve, centroid is the average of all vertex positions.
    """
    function centroid(c::AbstractDiscreteCurve{N,T}) where {N, T}
        !isclosed(c) && error("centroid requires a closed curve")
        sum(c.points) / nvertices(c)
    end

    """
        isoperimetric_quotient(c::AbstractDiscreteCurve{N,T})  →  T

    Compute the isoperimetric quotient of a closed curve.
    
    ```julia
    c = random_convex_polygon(8)
    q = isoperimetric_quotient(c)
    println("Isoperimetric quotient: $q")
    ```
    """
    function isoperimetric_quotient(c::AbstractDiscreteCurve{N,T}) where {N,T}
        A = abs(signed_area(c))
        L = arc_length(c)
        A > 0 && L > 0 ? (4 * π * A) / (L^2) : zero(T)
    end

    function regular_isoperimetric_quotient(n::Int)
       return π /n * cot(π/n)
    end
    function regular_isoperimetric_quotient(c::AbstractDiscreteCurve{N,T}) where {N,T}
       n = nvertices(c)
       return regular_isoperimetric_quotient(n)
    end

end
