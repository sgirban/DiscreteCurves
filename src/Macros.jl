module Macros
    using StaticArrays
    using ..AbstractTypes
    using ..CurveTypes
    using ..CurveTopology
    using ..Iterators
    export @curve,@where,@any_vertex,@all_vertices,@map_verts,@functional
    """
    @curve expr var ∈ t0..t1  [n=N] [closed[=true]]
 
Construct a `DiscreteCurve` from a parametric expression.
 
The expression can be a tuple `(fx, fy, ...)`, a vector `[fx, fy, ...]`,
or an `SVector(fx, fy, ...)`.
 
```julia
# 2D unit circle
c = @curve (cos(t), sin(t))  t ∈ 0..2π  n=128  closed
 
# 3D helix
c = @curve (cos(t), sin(t), t/(2π))  t ∈ 0..4π  n=256
 
# Lissajous
c = @curve (sin(3t), sin(2t))  t ∈ 0..2π  n=512  closed
 
# Trefoil knot
c = @curve (sin(t)+2sin(2t), cos(t)-2cos(2t), -sin(3t))  t ∈ 0..2π  n=512  closed
 
# Custom variable name
c = @curve (r*cos(θ), r*sin(θ)) θ ∈ 0..2π  n=64  closed
```
"""
macro curve(expr,range_expr,opts...)
    #  var ∈ t0..t1
    (range_expr.head == :call && 
        (range_expr.args[1] == :∈ || range_expr.args[1] == :in)) || 
        error("Expected range expression of the form `var ∈ t0..t1`,got $range_expr")
    var = range_expr.args[2]
    bounds = range_expr.args[3]
    (bounds.head == :call && bounds.args[1] == :..) || 
        error("Expected range expression of the form `var ∈ t0..t1`, got $range_expr")
    t0 = bounds.args[2]
    t1 = bounds.args[3]
    
    n_val = 128
    closed_val = false
    for opt in opts
        if opt isa Expr && opt.head == :(=) && opt.args[1] == :n
            n_val = opt.args[2]
        elseif opt == :closed
            closed_val = true
        elseif opt isa Expr && opt.head == :(=) && opt.args[1] == :closed
            closed_val = opt.args[2]
        end
    end
    sv_expr = _expr_to_svector(expr,var)
    quote
        ParametricCurve($var -> $sv_expr, $t0, $t1, $n_val; closed=$closed_val)
    end
end
function _expr_to_svector(expr, var)
    if expr.head == :tuple
        # (cos(t), sin(t)) → SVector(cos(t), sin(t))
        N = length(expr.args)
        :( StaticArrays.SVector{$N}($(expr.args...)) )
    elseif expr.head == :vect
        # [cos(t), sin(t)] → SVector(cos(t), sin(t))
        N = length(expr.args)
        :( StaticArrays.SVector{$N}($(expr.args...)) )
    elseif expr.head == :call && expr.args[1] == :SVector
        # Already SVector(...)  — pass through unchanged
        expr
    else
        # Single expression — assume 1D
        :( StaticArrays.SVector{1}($expr) )
    end
end

"""
    @where c p  expr
 
Return all vertices `p` of curve `c` satisfying `expr`.
 
```julia
# All vertices with radius > 0.5
pts = @where c p  norm(p) > 0.5
 
# Vertices in the upper half-plane
pts = @where c p  p[2] > 0
 
# Build a new open curve from the filtered vertices
c2 = DiscreteCurve(@where c p  p[1] > 0)
```
"""
macro where(curve_expr, var, pred_expr)
    quote
        filter($var -> $pred_expr, vertices($(esc(curve_expr))))
    end
end
 

"""
    @any_vertex c p  expr  →  Bool
 
Returns `true` if `expr` holds for at least one vertex `p` of `c`.
 
```julia
has_large = @any_vertex c p  norm(p) > 2.0
```
"""
macro any_vertex(curve_expr, var, pred_expr)
    quote
        any($var -> $pred_expr, vertices($(esc(curve_expr))))
    end
end

"""
    @all_vertices c p  expr  →  Bool
 
Returns `true` if `expr` holds for every vertex `p` of `c`.
 
```julia
all_unit = @all_vertices c p  abs(norm(p) - 1.0) < 1e-6
```
"""
macro all_vertices(curve_expr, var, pred_expr)
    quote
        all($var -> $pred_expr, vertices($(esc(curve_expr))))
    end
end

macro map_verts(curve_expr, var, expr)
    quote
        let _c = $(esc(curve_expr))
            DiscreteCurve(
                [let $var = p; $expr; end for p in vertices(_c)];
                closed = isclosed(_c))
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# @functional — Define and evaluate a curve functional
# ─────────────────────────────────────────────────────────────────────────────

"""
    @functional name c = expr
 
Define a named curve functional (a function curve → value).
 
```julia
@functional bending_energy c = sum(κ^2 for κ in curvatures(c)) / arc_length(c)
@functional centroid c = sum(vertices(c)) / nvertices(c)
 
bending_energy(c)   # evaluate
```
"""
macro functional(name, def)
    # def should be an assignment expression: c = expr
    if def.head != :(=)
        error("@functional: expected `name var = expr`, got $(def.head)")
    end
    arg = def.args[1]
    expr = def.args[2]
    quote
        $(esc(name))($(esc(arg))::AbstractDiscreteCurve) = $(esc(expr))
    end
end
 
end