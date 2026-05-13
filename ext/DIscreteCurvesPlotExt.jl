module DiscreteCurvesPlotExt

using DiscreteCurves
using DiscreteCurves.AbstractTypes: AbstractDiscreteCurve
using DiscreteCurves.CurveTypes: DiscreteCurve
using DiscreteCurves.CurveTopology: isclosed, nvertices, nedges
using DiscreteCurves.Curvature: TurningAngle, SteinerCurvature
using DiscreteCurves.Flows: FlowResult
using LinearAlgebra
using Plots

function _curve_xy(c::AbstractDiscreteCurve{2})
    n  = nvertices(c)
    xs = [c.points[i][1] for i in 1:n]
    ys = [c.points[i][2] for i in 1:n]
    if isclosed(c)
        push!(xs, xs[1]); push!(ys, ys[1])
    end
    xs, ys
end

function _curve_xyz(c::AbstractDiscreteCurve{3})
    n  = nvertices(c)
    xs = [c.points[i][1] for i in 1:n]
    ys = [c.points[i][2] for i in 1:n]
    zs = [c.points[i][3] for i in 1:n]
    if isclosed(c)
        push!(xs, xs[1]); push!(ys, ys[1]); push!(zs, zs[1])
    end
    xs, ys, zs
end

function _curvature_colors(c::AbstractDiscreteCurve{N,T}, cmap) where {N,T}
    κs = DiscreteCurves.curvatures(c, TurningAngle())
    n  = nvertices(c)
    
    κ_full = if isclosed(c)
        κs
    else
        [κs[1]; κs; κs[end]]
    end
    lo, hi = extrema(κ_full)
    rng    = hi - lo
    normed = rng ≈ 0 ? fill(0.5, n) : clamp.((κ_full .- lo) ./ rng, 0, 1)
    cmap, normed
end


@recipe function f(c::AbstractDiscreteCurve{2,T};
                   show_vertices = false,
                   color_by      = :none,
                   curvature_cmap = :plasma) where T

    aspect_ratio --> :equal
    label        --> ""
    linewidth    --> 2.0
    framestyle   --> :box

    xs, ys = _curve_xy(c)
    n      = nvertices(c)

    if color_by === :curvature
        cmap_fn = cgrad(curvature_cmap)
        κs = DiscreteCurves.curvatures(c, TurningAngle())
        κ_full = if isclosed(c)
            κs
        else
            [κs[1]; κs; κs[end]]
        end
        lo, hi = extrema(κ_full)
        rng = hi - lo

        for i in 1:length(xs)-1
            t  = rng ≈ 0 ? 0.5 : clamp((κ_full[i] - lo) / rng, 0, 1)
            col = cmap_fn[t]
            @series begin
                label --> ""
                linecolor --> col
                [xs[i], xs[i+1]], [ys[i], ys[i+1]]
            end
        end
        return 
    end

    if show_vertices
        @series begin
            seriestype  := :scatter
            markershape --> :circle
            markersize  --> 5
            label       --> ""
            xs[1:n], ys[1:n]
        end
    end

    xs, ys
end


@recipe function f(c::AbstractDiscreteCurve{3,T};
                   show_vertices = false) where T
    aspect_ratio --> :equal
    label        --> ""
    linewidth    --> 2.0

    xs, ys, zs = _curve_xyz(c)
    n          = nvertices(c)

    if show_vertices
        @series begin
            seriestype  := :scatter3d
            markershape --> :circle
            markersize  --> 5
            label       --> ""
            xs[1:n], ys[1:n], zs[1:n]
        end
    end

    seriestype := :path3d
    xs, ys, zs
end

@recipe function f(curves::AbstractVector{<:AbstractDiscreteCurve{2}};
                   color_by  = :index,
                   palette_  = :tab10,
                   alpha_min = 0.25)

    aspect_ratio --> :equal
    framestyle   --> :box
    legend       --> false

    nc   = length(curves)
    cpal = cgrad(palette_, nc; categorical = true)

    for (k, c) in enumerate(curves)
        xs, ys = _curve_xy(c)
        α = if color_by === :index
            alpha_min + (1 - alpha_min) * (k - 1) / max(nc - 1, 1)
        else
            1.0
        end
        @series begin
            label     --> ""
            linecolor --> cpal[k]
            alpha     --> α
            linewidth --> 1.5
            xs, ys
        end
    end
end

@recipe function f(curves::AbstractVector{<:AbstractDiscreteCurve{3}};
                   color_by  = :index,
                   palette_  = :tab10,
                   alpha_min = 0.25)

    framestyle --> :box
    legend     --> false
    nc   = length(curves)
    cpal = cgrad(palette_, nc; categorical = true)

    for (k, c) in enumerate(curves)
        xs, ys, zs = _curve_xyz(c)
        α = alpha_min + (1 - alpha_min) * (k - 1) / max(nc - 1, 1)
        @series begin
            seriestype := :path3d
            label      --> ""
            linecolor  --> cpal[k]
            alpha      --> α
            xs, ys, zs
        end
    end
end

@recipe function f(r::FlowResult{N};
                   show_tracking = false,
                   track_key     = :L,
                   palette_      = :viridis,
                   alpha_min     = 0.15,
                   show_final    = true) where N

    if isempty(r.curves)
       
        @series begin; r.curve; end
        return
    end

    if show_tracking && haskey(r, track_key)
        layout --> (1, 2)
    end

    nc   = length(r.curves)
    cpal = cgrad(palette_, nc)

    for (k, c) in enumerate(r.curves)
        t_norm = (k - 1) / max(nc - 1, 1)
        α      = alpha_min + (1 - alpha_min) * t_norm
        @series begin
            subplot   --> 1
            label     --> ""
            linecolor --> cpal[t_norm]
            alpha     --> α
            linewidth --> (k == nc && show_final) ? 2.5 : 1.2
            c
        end
    end

    if show_tracking && haskey(r, track_key)
        vals = r[track_key]
        steps = range(0, r.nsteps; length = length(vals))
        @series begin
            subplot   --> 2
            label     --> String(track_key)
            linewidth --> 2
            legend    --> :best
            xlabel    --> "step"
            ylabel    --> String(track_key)
            steps, vals
        end
    end
end


"""
    plot_tracking(r::FlowResult, keys...; kwargs...)

Plot tracked quantities from a `FlowResult` using Plots.jl.

```julia
using Plots
r = evolve(c, CurvatureFlow(); nsteps=500, track=[:L=>arc_length, :A=>signed_area])
plot_tracking(r, :L, :A)
```
"""
function DiscreteCurves.Graphics.plot_tracking_plots(
        r::FlowResult, keys...; kwargs...)
    isempty(keys) && (keys = filter(!=(:max_v), collect(keys(r.history))))
    isempty(keys) && error("No tracked quantities found in FlowResult")

    nk    = length(keys)
    steps = range(0, r.nsteps; length = length(r.max_v))

    plts = map(keys) do k
        haskey(r, k) || error("Key :$k not found in FlowResult")
        vals = r[k]
        sv   = range(0, r.nsteps; length = length(vals))
        Plots.plot(sv, vals;
            xlabel    = "step",
            ylabel    = String(k),
            title     = String(k),
            linewidth = 2,
            legend    = false,
            kwargs...)
    end
    nk == 1 ? plts[1] : Plots.plot(plts...; layout = (1, nk))
end

end