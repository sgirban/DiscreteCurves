module Orientation
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology
    using ..GeometricProperties

    export Orientation, CCW, CW, CounterClockwise, Clockwise
    export orient, orient!, is_ccw, is_cw

    """
        Orientation

    Enum for curve orientation:
    - `CCW`: Counter-clockwise (left-turn positive)
    - `CW`: Clockwise (right-turn positive)

    Aliases:
    - `CounterClockwise` = `CCW`
    - `Clockwise` = `CW`

    Mathematically, for a 2D closed curve:
    - CCW: signed area > 0
    - CW: signed area < 0
    """
    @enum Orientation CCW CW

    const CounterClockwise = CCW
    const Clockwise = CW

    """
        is_ccw(c::AbstractDiscreteCurve{2})  →  Bool

    Check if a 2D curve has counter-clockwise orientation.
    
    Works for both closed and open curves. For open curves, 
    checks the orientation of the polygon formed by closing with a line from last to first point.
    """
    function is_ccw(c::AbstractDiscreteCurve{2})
        signed_area(c) > 0
    end

    """
        is_cw(c::AbstractDiscreteCurve{2})  →  Bool

    Check if a 2D curve has clockwise orientation.
    
    Works for both closed and open curves. For open curves, 
    checks the orientation of the polygon formed by closing with a line from last to first point.
    """
    function is_cw(c::AbstractDiscreteCurve{2})
        signed_area(c) < 0
    end

    """
        orient(c::AbstractDiscreteCurve{2}; to::Orientation=CCW)  →  DiscreteCurve

    Return a new curve with the specified orientation (CCW or CW).
    
    For 2D closed curves only. Non-mutating.
    
    ```julia
    c = random_convex_polygon(10)
    c_ccw = orient(c; to=CCW)
    c_cw = orient(c; to=CW)
    ```
    """
    function orient(c::AbstractDiscreteCurve{2}; to::Orientation=CCW)
        !isclosed(c) && error("orient requires a closed curve")
        
        current_orientation = is_ccw(c) ? CCW : CW
        if current_orientation === to
            return DiscreteCurve(copy(c.points); closed=true)
        end
        reversed_pts = reverse(c.points)
        DiscreteCurve(reversed_pts; closed=true)
    end

    """
        orient!(c::DiscreteCurve{2}; to::Orientation=CCW)  →  DiscreteCurve

    In-place orientation adjustment. Modifies the curve in place.
    
    ```julia
    c = random_convex_polygon(10)
    orient!(c; to=CCW)
    ```
    """
    function orient!(c::DiscreteCurve{2}; to::Orientation=CCW)
        !isclosed(c) && error("orient! requires a closed curve")
        
        current_orientation = is_ccw(c) ? CCW : CW
        
        if current_orientation !== to
            reverse!(c.points)
        end
        
        c
    end

end
