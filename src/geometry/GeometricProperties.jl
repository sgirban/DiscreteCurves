module GeometricProperties
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators

    export signed_area, centroid

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
            # Closed curve: standard shoelace formula
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

end
