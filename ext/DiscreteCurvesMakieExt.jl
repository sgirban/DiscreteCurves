module DiscreteCurvesMakieExt

__precompile__(false)

using DiscreteCurves
using DiscreteCurves.AbstractTypes: AbstractDiscreteCurve
using DiscreteCurves.Graphics: CurveSkin, CurveFigure, VectorLayer, SKIN_NAMES, _parse_skin, BUILTIN_VERTEX_METRICS, SKIN_STUDIO, SKIN_DESIGNER, SKIN_GEOMETRIC, SKIN_SCIENTIFIC, SKIN_NEON, SKIN_PAPER, SKIN_PASTEL, SKIN_DARK, SKIN_WATERCOLOR, _require_makie
using Makie: Makie
using Makie: latexstring   # LaTeXString constructor for runtime-built LaTeX labels

using LinearAlgebra

using DiscreteCurves: AbstractTypes, CurveTypes, nvertices, nedges, isclosed, vertex, vertices

function _ensure_backend()
	if !isdefined(Makie, :current_backend) || Makie.current_backend() === nothing
		for backend in (:GLMakie, :CairoMakie, :WGLMakie)
			try
				eval(:(using $backend))
				return
			catch
				continue
			end
		end
	end
end

const SKINS = Dict{CurveSkin, NamedTuple}(
	# :studio
	SKIN_STUDIO => (
		background     = Makie.RGBf(1.0, 1.0, 1.0),
		curve_color    = Makie.RGBf(0.200, 0.490, 0.850),   # cerulean blue
		curve_width    = 3.0,
		colormap       = :blues,
		glow           = false,
		vertex_visible = true,
		vertex_color   = Makie.RGBf(0.095, 0.275, 0.553),   # deep navy, same hue
		vertex_size    = 15.0,
		axis_visible   = false,
		grid_visible   = false,
		text_color     = Makie.RGBf(0.20, 0.20, 0.25),
		palette        = [
		Makie.RGBf(0.200, 0.490, 0.850),   # cerulean
		Makie.RGBf(0.898, 0.337, 0.322),   # coral red
		Makie.RGBf(0.196, 0.706, 0.549),   # mint green
		Makie.RGBf(0.945, 0.612, 0.220),   # amber
		Makie.RGBf(0.647, 0.357, 0.820),   # soft violet
		Makie.RGBf(0.243, 0.659, 0.859)   # sky blue
		],
	),

	# :designer
	SKIN_DESIGNER => (
		background     = Makie.RGBf(0.027, 0.027, 0.067),
		curve_color    = :gradient,
		curve_width    = 3.0,
		colormap       = :plasma,
		glow           = true,
		glow_layers    = 4,
		glow_intensity = 0.13,
		vertex_visible = true,
		vertex_color   = :white,
		vertex_size    = 6.0,
		axis_visible   = false,
		grid_visible   = false,
		text_color     = Makie.RGBf(0.85, 0.85, 0.95),
		palette        = [:plasma, :viridis, :inferno, :magma],
	),

	# :geometric
	SKIN_GEOMETRIC => (
		background     = Makie.RGBf(0.98, 0.98, 0.98),
		curve_color    = Makie.RGBf(0.08, 0.08, 0.08),
		curve_width    = 2.0,
		colormap       = :grays,
		glow           = false,
		vertex_visible = true,
		vertex_color   = Makie.RGBf(0.08, 0.08, 0.08),
		vertex_size    = 6.5,
		axis_visible   = true,
		grid_visible   = true,
		text_color     = Makie.RGBf(0.1, 0.1, 0.1),
		palette        = [Makie.RGBf(0.08, 0.08, 0.08), Makie.RGBf(0.4, 0.4, 0.4), Makie.RGBf(0.65, 0.65, 0.65)],
	),

	# :scientific
	SKIN_SCIENTIFIC => (
		background     = :white,
		curve_color    = :palette,
		curve_width    = 1.8,
		colormap       = :viridis,
		glow           = false,
		vertex_visible = false,
		vertex_size    = 5.0,
		axis_visible   = true,
		grid_visible   = true,
		text_color     = :black,
		palette        = [
		Makie.RGBf(0.000, 0.447, 0.698),
		Makie.RGBf(0.835, 0.369, 0.000),
		Makie.RGBf(0.000, 0.620, 0.451),
		Makie.RGBf(0.941, 0.894, 0.259),
		Makie.RGBf(0.337, 0.706, 0.914),
		Makie.RGBf(0.800, 0.475, 0.655),
		Makie.RGBf(0.902, 0.624, 0.000)
		],
	),

	# :neon
	SKIN_NEON => (
		background     = Makie.RGBf(0.0, 0.0, 0.0),
		curve_color    = :palette,
		curve_width    = 2.0,
		colormap       = :rainbow,
		glow           = true,
		glow_layers    = 6,
		glow_intensity = 0.10,
		vertex_visible = true,
		vertex_size    = 5.5,
		axis_visible   = false,
		grid_visible   = false,
		text_color     = Makie.RGBf(0.9, 0.9, 0.9),
		palette        = [
		Makie.RGBf(0.0, 1.0, 1.0), Makie.RGBf(1.0, 0.0, 1.0),
		Makie.RGBf(1.0, 0.9, 0.0), Makie.RGBf(0.0, 1.0, 0.4),
		Makie.RGBf(1.0, 0.3, 0.0), Makie.RGBf(0.4, 0.4, 1.0)
		],
	),

	# :paper
	SKIN_PAPER => (
		background     = Makie.RGBf(0.961, 0.941, 0.902),
		curve_color    = Makie.RGBf(0.169, 0.110, 0.055),
		curve_width    = 1.9,
		colormap       = :copper,
		glow           = false,
		vertex_visible = true,
		vertex_color   = Makie.RGBf(0.169, 0.110, 0.055),
		vertex_size    = 5.5,
		axis_visible   = false,
		grid_visible   = false,
		text_color     = Makie.RGBf(0.169, 0.110, 0.055),
		palette        = [Makie.RGBf(0.169, 0.110, 0.055), Makie.RGBf(0.4, 0.2, 0.05), Makie.RGBf(0.6, 0.35, 0.1)],
	),

	# :pastel
	SKIN_PASTEL => (
		background     = :white,
		curve_color    = :gradient,
		curve_width    = 2.2,
		colormap       = :rainbow,
		glow           = false,
		vertex_visible = true,
		vertex_size    = 6.0,
		axis_visible   = false,
		grid_visible   = false,
		text_color     = Makie.RGBf(0.4, 0.4, 0.4),
		palette        = [Makie.RGBf(0.686, 0.769, 0.929), Makie.RGBf(0.961, 0.702, 0.749),
		Makie.RGBf(0.729, 0.902, 0.745), Makie.RGBf(0.988, 0.886, 0.694)],
	),

	# :dark
	SKIN_DARK => (
		background     = Makie.RGBf(0.118, 0.118, 0.180),
		curve_color    = :gradient,
		curve_width    = 2.3,
		colormap       = :cool,
		glow           = true,
		glow_layers    = 3,
		glow_intensity = 0.10,
		vertex_visible = true,
		vertex_color   = Makie.RGBf(0.537, 0.706, 0.980),
		vertex_size    = 6.0,
		axis_visible   = true,
		grid_visible   = true,
		text_color     = Makie.RGBf(0.804, 0.839, 0.957),
		palette        = [
		Makie.RGBf(0.537, 0.706, 0.980), Makie.RGBf(0.953, 0.545, 0.659),
		Makie.RGBf(0.647, 0.902, 0.678), Makie.RGBf(0.988, 0.843, 0.498)
		],
	),

	# :watercolor
	SKIN_WATERCOLOR => (
		background     = :white,
		curve_color    = :gradient,
		curve_width    = 3.0,
		curve_alpha    = 0.6,
		colormap       = :turbo,
		glow           = false,
		vertex_visible = true,
		vertex_size    = 8.0,
		vertex_alpha   = 0.5,
		axis_visible   = false,
		grid_visible   = false,
		text_color     = Makie.RGBf(0.3, 0.3, 0.3),
		palette        = [Makie.RGBAf(0.2, 0.5, 0.9, 0.6), Makie.RGBAf(0.9, 0.2, 0.5, 0.6),
		Makie.RGBAf(0.2, 0.8, 0.4, 0.6)],
	),
)

function _resolve(skin::CurveSkin, user_kw)
	sk = get(SKINS, skin, NamedTuple())

	function opt(key, default)
		haskey(user_kw, key) && return user_kw[key]
		hasfield(typeof(sk), key) && return getfield(sk, key)
		return default
	end

	(
		skin                = skin,
		background          = opt(:background, :white),
		curve_color         = opt(:curve_color, Makie.RGBf(0.2, 0.49, 0.85)),
		curve_width         = opt(:curve_width, 2.0),
		curve_alpha         = opt(:curve_alpha, 1.0),
		curve_style         = opt(:curve_style, :solid),
		colormap            = opt(:colormap, :plasma),
		clims               = opt(:clims, :auto),
		glow                = opt(:glow, false),
		glow_layers         = opt(:glow_layers, 4),
		glow_intensity      = opt(:glow_intensity, 0.12),
		vertex_visible      = opt(:vertex_visible, false),
		vertex_color        = opt(:vertex_color, nothing),
		vertex_size         = opt(:vertex_size, 14.0),
		vertex_marker       = opt(:vertex_marker, :circle),
		vertex_alpha        = opt(:vertex_alpha, 1.0),
		vertex_labels       = opt(:vertex_labels, false),
		label_position      = opt(:label_position, :N),
		label_spacing       = opt(:label_spacing, 0.25),   # fraction of bounding-box extent (increased significantly for visibility)
		label_fontsize      = opt(:label_fontsize, 16),
		label_color         = opt(:label_color, :auto),
		edge_labels         = opt(:edge_labels, false),
		edge_label_position = opt(:edge_label_position, :above),
		edge_label_spacing  = opt(:edge_label_spacing, 0.20),
		tooltip             = opt(:tooltip, false),
		vectors             = opt(:vectors, nothing),
		vector_scale        = opt(:vector_scale, 0.15),
		vector_width        = opt(:vector_width, 1.5),
		vector_color        = opt(:vector_color, :auto),
		vector_colormap     = opt(:vector_colormap, :viridis),
		axis_visible        = opt(:axis_visible, true),
		grid_visible        = opt(:grid_visible, true),
		aspect              = opt(:aspect, :equal),
		title               = opt(:title, ""),
		xlabel              = opt(:xlabel, ""),
		ylabel              = opt(:ylabel, ""),
		legend              = opt(:legend, false),
		curve_labels        = opt(:curve_labels, nothing),
		size                = opt(:size, (900, 700)),
		text_color          = opt(:text_color, :black),
		palette             = opt(:palette, _default_palette()),
		grid_color          = opt(:grid_color, Makie.RGBAf(0, 0, 0, 0.12)),
		controls            = opt(:controls, false),
		mouse_mode          = opt(:mouse_mode, :default),
		curve_name          = opt(:curve_name, "\\gamma"),
	)
end

_default_palette() = [
	Makie.RGBf(0.2, 0.49, 0.85), Makie.RGBf(0.9, 0.34, 0.32),
	Makie.RGBf(0.2, 0.71, 0.55), Makie.RGBf(0.94, 0.61, 0.22),
	Makie.RGBf(0.65, 0.36, 0.82), Makie.RGBf(0.24, 0.66, 0.86),
]

function _eval_vertex_colors(c, spec, cmap_sym, clims)
	n  = nvertices(c)
	cm = Makie.to_colormap(cmap_sym)

	if spec === nothing
		return fill(Makie.RGBf(0.2, 0.49, 0.85), n), nothing
	elseif (spec isa Symbol && spec ∈ BUILTIN_VERTEX_METRICS) || spec == :gradient
		vals = _get_metric(c, spec)
		return _scalars_to_colors(vals, cm, clims), vals
	elseif spec isa Pair
		metric, cmap2 = spec
		cm2 = Makie.to_colormap(cmap2)
		vals = _get_metric(c, metric)
		return _scalars_to_colors(vals, cm2, clims), vals
	elseif spec isa Function
		vals = _eval_func_colors(c, spec)
		return _scalars_to_colors(vals, cm, clims), vals
	elseif spec isa AbstractVector{<:Real}
		length(spec) == n || throw(DimensionMismatch("color vector length $(length(spec)) ≠ nvertices $n"))
		vals = collect(Float64, spec)
		return _scalars_to_colors(vals, cm, clims), vals
	elseif spec isa AbstractVector
		length(spec) == n || throw(DimensionMismatch("color vector length $(length(spec)) ≠ nvertices $n"))
		return Makie.to_color.(spec), nothing
	else
		col = Makie.to_color(spec)
		return fill(col, n), nothing
	end
end

function _get_metric(c, spec)
	n = nvertices(c)
	if spec == :curvature || spec == :gradient
		κ = curvatures(c, TurningAngle())
		return isclosed(c) ? Float64.(κ) : Float64.([κ[1]; κ; κ[end]])
	elseif spec == :arc_length
		return Float64.(cumulative_length(c))
	elseif spec == :arc_length_norm
		s = cumulative_length(c);
		return Float64.(s ./ s[end])
	elseif spec == :turning_angle
		θ = turning_angles(c; signed = false)
		return isclosed(c) ? Float64.(θ) : Float64.([θ[1]; θ; θ[end]])
	elseif spec == :vertex_index
		return Float64.(1:n)
	elseif spec == :radius
		ctr = centroid(c)
		return Float64.([norm(vertex(c, i)-ctr) for i in 1:n])
	elseif spec == :normal_x
		return Float64.([normals(c)[i][1] for i in 1:n])
	elseif spec == :normal_y
		return Float64.([normals(c)[i][2] for i in 1:n])
	else
		error("Unknown metric :$spec")
	end
end

function _eval_func_colors(c, f)
	n = nvertices(c)
	vals = Vector{Float64}(undef, n)
	# Try (c, i) signature first, fall back to (p)
	try
		vals[1] = Float64(f(c, 1))
		for i in 2:n
			;
			vals[i] = Float64(f(c, i));
		end
	catch
		for i in 1:n
			;
			vals[i] = Float64(f(vertex(c, i)));
		end
	end
	vals
end

function _scalars_to_colors(vals::Vector{Float64}, cm, clims)
	lo, hi = clims == :auto ? extrema(vals) : Float64.(clims)
	rng    = hi - lo
	normed = rng ≈ 0.0 ? fill(0.5, length(vals)) : clamp.((vals .- lo) ./ rng, 0.0, 1.0)
	cm_vec = Makie.to_colormap(cm)          # Vector{RGBAf} — stable across all versions
	[_sample_colormap(cm_vec, t) for t in normed]
end

@inline function _sample_colormap(cm::AbstractVector, t::Real)
	n  = length(cm)
	n == 1 && return Makie.RGBAf(cm[1])
	s  = clamp(Float64(t), 0.0, 1.0) * (n - 1)
	i0 = clamp(floor(Int, s) + 1, 1, n)
	i1 = min(i0 + 1, n)
	α  = s - floor(s)
	c0 = Makie.RGBAf(cm[i0])
	c1 = Makie.RGBAf(cm[i1])
	Makie.RGBAf(c0.r + α*(c1.r - c0.r),
	            c0.g + α*(c1.g - c0.g),
	            c0.b + α*(c1.b - c0.b),
	            c0.alpha + α*(c1.alpha - c0.alpha))
end
function _make_figure_axis(opt, is3d::Bool)
	fig = Makie.Figure(; size = opt.size, backgroundcolor = opt.background)
	ax  = if is3d
		Makie.Axis3(fig[1, 1]; backgroundcolor = opt.background, title = opt.title,
		titlecolor = opt.text_color)
	else
		_tc = opt.text_color
		Makie.Axis(fig[1, 1];
		backgroundcolor    = opt.background,
		aspect             = opt.aspect == :equal ? Makie.DataAspect() : Makie.AxisAspect(1),
		title              = opt.title, titlecolor      = _tc,
		xlabel             = opt.xlabel, ylabel          = opt.ylabel,
		xlabelcolor        = _tc, ylabelcolor     = _tc,
		xticklabelcolor    = _tc, yticklabelcolor = _tc,
		xtickcolor         = _tc, ytickcolor      = _tc,
		# ── Generous margins so vertex labels are never clipped ──────────────
		# Default is 0.05 (5%). 0.15 (15%) gives labels room on all sides.
		# This affects the autolimits calculation, not the clip region itself.
		xautolimitmargin   = (0.15f0, 0.15f0),
		yautolimitmargin   = (0.15f0, 0.15f0),
		)
	end

	if !opt.axis_visible
		Makie.hidedecorations!(ax)
		Makie.hidespines!(ax)
	end
	if !opt.grid_visible && !is3d
		ax.xgridvisible = false
		ax.ygridvisible = false
	end

	if !is3d
		_apply_mouse_mode!(ax, opt.mouse_mode)
	end

	if opt.controls
		_add_controls!(fig, ax, is3d, opt)
	end

	fig, ax
end


function _apply_mouse_mode!(ax::Makie.Axis, mode::Symbol)
	mode == :default && return 

	if mode == :pan
		
		try
			Makie.deregister_interaction!(ax, :scrollzoom)
		catch; end
		try
			Makie.register_interaction!(ax, :scrollpan) do event::Makie.ScrollEvent, ax
				dx = event.x * 0.05
				dy = event.y * 0.05
				xmin, xmax = ax.xaxis.attributes.limits[]
				ymin, ymax = ax.yaxis.attributes.limits[]
				rng_x = xmax - xmin
				rng_y = ymax - ymin
				Makie.xlims!(ax, xmin - dx*rng_x, xmax - dx*rng_x)
				Makie.ylims!(ax, ymin - dy*rng_y, ymax - dy*rng_y)
				Makie.consume(event)
			end
		catch; end

	elseif mode == :zoom
		try
			Makie.deregister_interaction!(ax, :rectanglezoom)
		catch; end

	else
		@warn "Unknown mouse_mode :$mode. Use :default, :pan, or :zoom."
	end
end


function _add_controls!(fig, ax, is3d::Bool, opt)
	isdefined(Makie, :Button) || return

	bg  = Makie.to_color(opt.background)
	tc  = Makie.to_color(opt.text_color)
	btn_bg = Makie.RGBf(
		clamp(bg.r + 0.08f0, 0f0, 1f0),
		clamp(bg.g + 0.08f0, 0f0, 1f0),
		clamp(bg.b + 0.08f0, 0f0, 1f0),
	)
	btn_bg_hover = Makie.RGBf(
		clamp(bg.r + 0.18f0, 0f0, 1f0),
		clamp(bg.g + 0.18f0, 0f0, 1f0),
		clamp(bg.b + 0.18f0, 0f0, 1f0),
	)

	toolbar = fig[2, 1] = Makie.GridLayout(tellwidth = false)

	reset_btn = Makie.Button(toolbar[1, 1];
		label        = "Reset",
		buttoncolor  = btn_bg,
		buttoncolor_hover   = btn_bg_hover,
		buttoncolor_active  = btn_bg_hover,
		labelcolor   = tc,
		height       = 28,
		width        = 80,
	)
	Makie.on(reset_btn.clicks) do _
		Makie.reset_limits!(ax)
	end

	if !is3d
		pan_btn = Makie.Button(toolbar[1, 2];
			label       = "Pan",
			buttoncolor = btn_bg,
			buttoncolor_hover  = btn_bg_hover,
			buttoncolor_active = btn_bg_hover,
			labelcolor  = tc,
			height = 28, width = 80,
		)
		Makie.on(pan_btn.clicks) do _
			_apply_mouse_mode!(ax, :pan)
		end

		zoom_btn = Makie.Button(toolbar[1, 3];
			label       = "Zoom",
			buttoncolor = btn_bg,
			buttoncolor_hover  = btn_bg_hover,
			buttoncolor_active = btn_bg_hover,
			labelcolor  = tc,
			height = 28, width = 80,
		)
		Makie.on(zoom_btn.clicks) do _
			_apply_mouse_mode!(ax, :zoom)
		end

		fit_btn = Makie.Button(toolbar[1, 4];
			label       = "Fit",
			buttoncolor = btn_bg,
			buttoncolor_hover  = btn_bg_hover,
			buttoncolor_active = btn_bg_hover,
			labelcolor  = tc,
			height = 28, width = 80,
		)
		Makie.on(fit_btn.clicks) do _
			Makie.autolimits!(ax)
		end
	else
		rot_btn = Makie.Button(toolbar[1, 2];
			label       = "Rotate: drag",
			buttoncolor = btn_bg,
			labelcolor  = tc,
			height = 28, width = 120,
		)
		Makie.on(rot_btn.clicks) do _; end   # informational only
	end

	Makie.rowsize!(fig.layout, 2, Makie.Auto())
	Makie.rowgap!(fig.layout, 4)
end

function _draw_curve!(ax, c, vertex_colors, opt, is3d::Bool)
	n     = nvertices(c)
	w     = Float32(isa(opt.curve_width, Function) ? opt.curve_width(c) : opt.curve_width)
	alpha = Float32(opt.curve_alpha)

	xs = [c.points[i][1] for i in 1:n]
	ys = [c.points[i][2] for i in 1:n]
	cs = copy(vertex_colors)
	zs = is3d ? [c.points[i][3] for i in 1:n] : nothing
	if isclosed(c)
		push!(xs, xs[1]);
		push!(ys, ys[1]);
		push!(cs, cs[1])
		is3d && push!(zs, zs[1])
	end

	draw_cs = alpha < 1.0f0 ? [(Makie.to_color(col), alpha) for col in cs] : cs
	pts     = is3d ? Makie.Point3f.(xs, ys, zs) : Makie.Point2f.(xs, ys)

	if opt.glow
		for k in opt.glow_layers:-1:1
			gα = Float32(opt.glow_intensity / k)
			gcols = [(Makie.to_color(col), gα) for col in cs]
			Makie.lines!(ax, pts; color = gcols, linewidth = w + Float32(k*3.5),
				linestyle = opt.curve_style)
		end
	end

	Makie.lines!(ax, pts; color = draw_cs, linewidth = w, linestyle = opt.curve_style)
end

function _draw_vertices!(ax, c, vertex_colors, opt, is3d::Bool; tooltip_fn = nothing)
	opt.vertex_visible || return
	n      = nvertices(c)
	vs     = Float32(opt.vertex_size)
	valpha = Float32(opt.vertex_alpha)

	vcols = if opt.vertex_color === nothing || opt.vertex_color == :auto
		vertex_colors[1:n]
	else
		vc, _ = _eval_vertex_colors(c, opt.vertex_color, opt.colormap, opt.clims)
		vc
	end

	final_cols = valpha < 1.0f0 ? [(Makie.to_color(col), valpha) for col in vcols] : vcols

	pts = if is3d
		Makie.Point3f.([c.points[i][1] for i in 1:n],
			[c.points[i][2] for i in 1:n],
			[c.points[i][3] for i in 1:n])
	else
		Makie.Point2f.([c.points[i][1] for i in 1:n],
			[c.points[i][2] for i in 1:n])
	end

	if opt.glow
		for k in 2:-1:1
			gα = Float32(opt.glow_intensity * 2.0 / k)
			Makie.scatter!(ax, pts; color = [(Makie.to_color(col), gα) for col in vcols],
				markersize = vs+Float32(k*4), marker = opt.vertex_marker)
		end
	end

	scatter_plt = Makie.scatter!(ax, pts; color = final_cols, markersize = vs,
		marker = opt.vertex_marker)

	if tooltip_fn !== nothing
		scatter_plt.inspector_label = (self, idx, pos) -> tooltip_fn(c, idx)
		Makie.DataInspector(ax)
	end
end

function _label_offset(position::Symbol, spacing::Float64, n_dir::Tuple)
	dx, dy = n_dir
	if position == :N
		;
		return 0.0, spacing
	elseif position == :S
		;
		return 0.0, -spacing
	elseif position == :E
		;
		return spacing, 0.0
	elseif position == :W
		;
		return -spacing, 0.0
	else
		return dx*spacing, dy*spacing
	end
end

function _draw_labels!(ax, c, opt, is3d::Bool)
	opt.vertex_labels == false && return
	n    = nvertices(c)
	lcol = opt.label_color == :auto ? opt.text_color : opt.label_color
	lcol = Makie.to_color(lcol)  # Convert to actual Makie color (handles both colors and symbols)
	pos  = opt.label_position
	cname = String(opt.curve_name)    # e.g. "\\gamma" or "p"

	# ── Offset: fraction of bounding-box extent ────────────────────────────────
	lo_bb, hi_bb = DiscreteCurves.bounding_box(c)
	extent = max(Float64(hi_bb[1] - lo_bb[1]), Float64(hi_bb[2] - lo_bb[2]), 0.01)
	sp     = Float64(opt.label_spacing) * extent

	al = _pos_to_align(pos)

	# ── Precompute normals for :auto mode; ignored for :N/:S/:E/:W ────────────
	ns = DiscreteCurves.normals(c)

	# ── Build and draw labels one by one ──────────────────────────────────────
	# We use per-element calls rather than the vectorised form because:
	#   - The vectorised Makie.text!(ax, strings; position=points, ...) form
	#     changed between Makie versions and can silently render all strings
	#     as a single multi-line block at one position.
	#   - The single-element form  text!(ax, pt; text="...", ...)
	#     is stable and unambiguous in every Makie v0.19–0.21+ release.

	for i in 1:n
		p  = c.points[i]
		nv = ns[i]
		dx, dy = _label_offset(pos, sp, (Float64(nv[1]), Float64(nv[2])))
		pt  = is3d ? Makie.Point3f(p[1]+dx, p[2]+dy, p[3]) :
		             Makie.Point2f(p[1]+dx, p[2]+dy)

		label = if opt.vertex_labels == true || opt.vertex_labels == :index
			string(i)
		elseif opt.vertex_labels == :coords
			if is3d
				"($(round(Float64(p[1]);digits=2)), $(round(Float64(p[2]);digits=2)), $(round(Float64(p[3]);digits=2)))"
			else
				"($(round(Float64(p[1]);digits=2)), $(round(Float64(p[2]);digits=2)))"
			end
		elseif opt.vertex_labels == :latex || opt.vertex_labels == :gamma
			# Produce LaTeX "γ_i" with the user-chosen curve name.
			# We build a proper LaTeXString so subscript renders correctly.
			Makie.latexstring("$(cname)_{$(i)}")
		elseif opt.vertex_labels isa Function
			string(opt.vertex_labels(i))
		else
			string(i)
		end

		Makie.text!(ax, label;
			position = pt,
			color    = lcol,
			fontsize = Float32(opt.label_fontsize),
			align    = al,
		)
	end
end

# Derive Makie text alignment from cardinal direction symbol.
function _pos_to_align(pos::Symbol)
	pos == :N && return (:center, :bottom)
	pos == :S && return (:center, :top)
	pos == :E && return (:left,   :center)
	pos == :W && return (:right,  :center)
	return (:center, :bottom)
end

function _draw_edge_labels!(ax, c, opt, is3d::Bool)
	opt.edge_labels == false && return
	n  = nvertices(c)
	ne = nedges(c)
	lcol = opt.label_color == :auto ? opt.text_color : opt.label_color
	lcol = Makie.to_color(lcol)  # Convert to actual Makie color

	lo_bb, hi_bb = DiscreteCurves.bounding_box(c)
	extent = max(Float64(hi_bb[1] - lo_bb[1]), Float64(hi_bb[2] - lo_bb[2]), 0.01)
	sp  = Float64(opt.edge_label_spacing) * extent
	dir = opt.edge_label_position == :above ? 1.0 : -1.0
	ns  = DiscreteCurves.normals(c)

	for i in 1:ne
		p1  = c.points[i]
		p2  = c.points[mod1(i+1, n)]
		mid = (p1 + p2) / 2
		nv  = (ns[i] + ns[mod1(i+1,n)]) / 2
		nn  = norm(nv)
		nhat = nn > 1e-8 ? nv/nn : typeof(nv)(0,1)
		tx  = Float64(mid[1]) + Float64(nhat[1])*sp*dir
		ty  = Float64(mid[2]) + Float64(nhat[2])*sp*dir

		label = if opt.edge_labels == true
			string(i)
		elseif opt.edge_labels isa Function
			string(opt.edge_labels(c, i))
		elseif opt.edge_labels isa AbstractVector && i ≤ length(opt.edge_labels)
			string(opt.edge_labels[i])
		else
			string(i)
		end

		pt = is3d ? Makie.Point3f(tx, ty, Float64(mid[3])) : Makie.Point2f(tx, ty)
		Makie.text!(ax, label;
			position = pt,
			color    = lcol,
			fontsize = Float32(opt.label_fontsize),
			align    = (:center, :center),
		)
	end
end

function _draw_vectors_overlay!(ax, c, vertex_colors, opt, is3d::Bool)
	opt.vectors === nothing && return
	n = nvertices(c)
	L̄ = arc_length(c) / n
	sc = Float64(opt.vector_scale) * L̄

	vecs = if opt.vectors == :normals
		;
		normals(c)
	elseif opt.vectors == :tangents
		;
		unit_tangents(c; at = :vertex)
	elseif opt.vectors == :curvature_vectors
		;
		curvature_vectors(c)
	elseif opt.vectors isa AbstractVector
		;
		opt.vectors
	else
		error("vectors must be :normals, :tangents, :curvature_vectors or Vector")
	end

	arrow_cols = if opt.vector_color == :auto || opt.vector_color == :matching
		vertex_colors[1:n]
	elseif opt.vector_color == :magnitude
		mags = Float64[norm(v) for v in vecs]
		cm2  = Makie.to_colormap(opt.vector_colormap)
		_scalars_to_colors(mags, cm2, :auto)
	else
		fill(Makie.to_color(opt.vector_color), n)
	end

	w = Float32(opt.vector_width)
	# Proportional arrowhead: tipwidth and tiplength scale with shaft width.
	# arrows2d!/arrows3d! are the Makie v0.21+ replacements for the deprecated arrows!.
	tip_w = w * 1.8f0   # arrowhead base width
	tip_l = w * 2.8f0   # arrowhead length

	if is3d
		ps = Makie.Point3f.([c.points[i][1] for i in 1:n], [c.points[i][2] for i in 1:n], [c.points[i][3] for i in 1:n])
		ds = Makie.Point3f.([vecs[i][1]*sc for i in 1:n], [vecs[i][2]*sc for i in 1:n], [vecs[i][3]*sc for i in 1:n])
		Makie.arrows3d!(ax, ps, ds; color = arrow_cols, linewidth = w,
		                tipwidth = tip_w, tiplength = tip_l)
	else
		ps = Makie.Point2f.([c.points[i][1] for i in 1:n], [c.points[i][2] for i in 1:n])
		ds = Makie.Point2f.([vecs[i][1]*sc for i in 1:n], [vecs[i][2]*sc for i in 1:n])
		Makie.arrows2d!(ax, ps, ds; color = arrow_cols, linewidth = w,
		                tipwidth = tip_w, tiplength = tip_l)
	end
end

# ─────────────────────────────────────────────────────────────────────────────
#  Tooltip builder
# ─────────────────────────────────────────────────────────────────────────────

function _build_tooltip_fn(c, tooltip_spec)
	tooltip_spec == false && return nothing
	tooltip_spec == true && return _build_tooltip_fn(c, [:curvature, :arc_length])

	if tooltip_spec isa AbstractVector{Symbol}
		κs = :curvature ∈ tooltip_spec ? curvatures(c, TurningAngle()) : nothing
		ss = :arc_length ∈ tooltip_spec ? cumulative_length(c) : nothing
		θs = :turning_angle ∈ tooltip_spec ? turning_angles(c; signed = false) : nothing

		return function (c2, i)
			lines = ["vertex $i  pos = ($(round(c2.points[i][1];digits=4)), $(round(c2.points[i][2];digits=4)))"]
			κs !== nothing && i ≤ length(κs) && push!(lines, "κ = $(round(κs[i];sigdigits=5))")
			ss !== nothing && push!(lines, "s = $(round(ss[i];sigdigits=5))")
			θs !== nothing && i ≤ length(θs) && push!(lines, "θ = $(round(θs[i];sigdigits=5))")
			join(lines, "\n")
		end
	elseif tooltip_spec isa Function
		return (c2, i) -> tooltip_spec(c2, i)
	end
	return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Draw one complete curve
# ─────────────────────────────────────────────────────────────────────────────

function _plot_one_curve!(ax, c, opt, cspec, is3d::Bool, palette_idx::Int)
	# Resolve color spec (palette cycling for multi-curve)
	actual_spec = if cspec == :palette
		pal = opt.palette
		pal[mod1(palette_idx, length(pal))]
	else
		cspec
	end

	vertex_colors, scalars = _eval_vertex_colors(c, actual_spec, opt.colormap, opt.clims)

	_draw_curve!(ax, c, vertex_colors, opt, is3d)

	tooltip_fn = _build_tooltip_fn(c, opt.tooltip)
	_draw_vertices!(ax, c, vertex_colors, opt, is3d; tooltip_fn)

	_draw_labels!(ax, c, opt, is3d)
	_draw_edge_labels!(ax, c, opt, is3d)
	_draw_vectors_overlay!(ax, c, vertex_colors, opt, is3d)

	vertex_colors
end
# ─────────────────────────────────────────────────────────────────────────────
#  curveplot
# ─────────────────────────────────────────────────────────────────────────────

function DiscreteCurves.Graphics.curveplot(cs::AbstractVector{<:AbstractDiscreteCurve};
	skin = :studio, kwargs...)
	_ensure_backend()
	skin_val = _parse_skin(skin)
	opt      = _resolve(skin_val, kwargs)
	is3d     = any(c -> length(c.points[1]) == 3, cs)
	fig, ax  = _make_figure_axis(opt, is3d)

	for (k, c) in enumerate(cs)
		_plot_one_curve!(ax, c, opt, opt.curve_color, is3d, k)
	end

	if opt.legend && opt.curve_labels !== nothing
		Makie.Legend(fig[1, 2], ax, opt.curve_labels)
	end

	CurveFigure(fig, ax, skin_val, opt)
end

function DiscreteCurves.Graphics.curveplot(c::AbstractDiscreteCurve; kwargs...)
	DiscreteCurves.Graphics.curveplot([c]; kwargs...)
end

function DiscreteCurves.Graphics.curveplot(c1::AbstractDiscreteCurve,
	rest::AbstractDiscreteCurve...; kwargs...)
	DiscreteCurves.Graphics.curveplot([c1, rest...]; kwargs...)
end

function DiscreteCurves.Graphics.curveplot!(ax::Union{Makie.Axis, Makie.Axis3},
	c::AbstractDiscreteCurve; skin = :studio, kwargs...)
	opt  = _resolve(_parse_skin(skin), kwargs)
	is3d = length(c.points[1]) == 3
	_plot_one_curve!(ax, c, opt, opt.curve_color, is3d, 1)
	ax
end

function DiscreteCurves.Graphics.curveplot!(ax::Union{Makie.Axis, Makie.Axis3},
	cs::AbstractVector{<:AbstractDiscreteCurve};
	skin = :studio, kwargs...)
	opt  = _resolve(_parse_skin(skin), kwargs)
	is3d = any(c -> length(c.points[1]) == 3, cs)
	for (k, c) in enumerate(cs)
		_plot_one_curve!(ax, c, opt, opt.curve_color, is3d, k)
	end
	ax
end

# Makie's own plot() — so `plot(c)` works
Makie.plot(c::AbstractDiscreteCurve; kwargs...) =
	DiscreteCurves.Graphics.curveplot(c; kwargs...)

# ─────────────────────────────────────────────────────────────────────────────
#  CurveFigure display + Base.:+
# ─────────────────────────────────────────────────────────────────────────────

# Transparent display — REPL/Jupyter see the Makie Figure
Base.show(io::IO, m::MIME, cf::CurveFigure) = show(io, m, cf._fig)
Base.showable(m::MIME, cf::CurveFigure)     = showable(m, cf._fig)
Base.display(cf::CurveFigure)               = (_ensure_backend(); display(cf._fig))

# ── Base.:+ — overlay VectorLayer onto CurveFigure ───────────────────────────
function Base.:+(cf::CurveFigure, vl::VectorLayer)
	opt  = cf._opt
	ax   = cf._ax
	is3d = ax isa Makie.Axis3
	_render_vector_layer!(ax, vl, opt, is3d)
	cf
end

# ─────────────────────────────────────────────────────────────────────────────
# save_figure
# ─────────────────────────────────────────────────────────────────────────────

function DiscreteCurves.Graphics.save_figure(filename::AbstractString,
	cf::CurveFigure; kwargs...)
	_ensure_backend()
	Makie.save(filename, cf._fig; kwargs...)
	@info "Saved figure → $filename"
	filename
end

# ─────────────────────────────────────────────────────────────────────────────
# vectorplot / vectorplot!
# ─────────────────────────────────────────────────────────────────────────────

"""
	vectorplot(vecs; anchors, skin, kwargs...)   → CurveFigure (standalone)
	vectorplot(vecs; kwargs...)                  → VectorLayer (for composition with +)
 
When called without `anchors` or called for composition, returns a VectorLayer.
When called standalone (not composing), creates its own figure.
"""
function DiscreteCurves.Graphics.vectorplot(vecs; skin = :auto, kwargs...)

	if skin == :auto
		return VectorLayer(vecs, get(kwargs, :anchors, nothing), kwargs)
	end
	_ensure_backend()
	skin_val = _parse_skin(skin == :auto ? :studio : skin)
	opt      = _resolve(skin_val, kwargs)
	n_vecs   = vecs isa AbstractVector ? length(vecs) : 1
	is3d     = _vecs_is3d(vecs)
	fig, ax  = _make_figure_axis(opt, is3d)
	anchors  = get(kwargs, :anchors, nothing)
	_render_vector_layer!(ax, VectorLayer(vecs, anchors, kwargs), opt, is3d)
	CurveFigure(fig, ax, skin_val, opt)
end


DiscreteCurves.Graphics.vectorplot(vecs, ::Nothing; kwargs...) =
	VectorLayer(vecs, get(kwargs, :anchors, nothing), kwargs)

function DiscreteCurves.Graphics.vectorplot!(ax::Union{Makie.Axis, Makie.Axis3},
	vecs; skin = :studio, kwargs...)
	opt = _resolve(_parse_skin(skin), kwargs)
	is3d = ax isa Makie.Axis3
	anchors = get(kwargs, :anchors, nothing)
	_render_vector_layer!(ax, VectorLayer(vecs, anchors, kwargs), opt, is3d)
	ax
end


function _render_vector_layer!(ax, vl::VectorLayer, opt, is3d::Bool)
	vecs    = _normalise_vecs(vl.vectors)
	anchors = _resolve_anchors(vl.anchors, vecs)
	kw      = vl.kwargs

	N = length(vecs[1])

	arrow_color     = get(kw, :arrow_color, Makie.RGBf(0.200, 0.490, 0.850))
	arrow_width     = Float32(get(kw, :arrow_width, 5.0))
	arrowhead_scale = Float32(get(kw, :arrowhead_scale, 6.0))
	lengthscale     = Float32(get(kw, :lengthscale, 1.2))
	anchor_visible  = Bool(get(kw, :anchor_visible, true))
	anchor_size     = Float32(get(kw, :anchor_size, arrow_width * 3.0))
	anchor_marker   = get(kw, :anchor_marker, :circle)
	anchor_color    = get(kw, :anchor_color, nothing)
	normalize_vecs  = Bool(get(kw, :normalize, false))
	scale           = Float64(get(kw, :scale, 1.0))
	cmap_sym        = get(kw, :colormap, :viridis)
	labels          = get(kw, :labels, false)
	label_pos       = get(kw, :label_position, :N)
	label_sp        = Float64(get(kw, :label_spacing, 0.12))
	label_fs        = get(kw, :label_fontsize, 11)
	label_color     = get(kw, :label_color, opt.text_color)
	label_color     = Makie.to_color(label_color)  # Convert to actual Makie color

	n = length(vecs)
	mags = [norm(v) for v in vecs]

	# Resolve colors
	draw_cols = if arrow_color == :magnitude
		cm = Makie.to_colormap(cmap_sym)
		_scalars_to_colors(Float64.(mags), cm, :auto)
	elseif arrow_color isa AbstractVector
		length(arrow_color) == n ? Makie.to_color.(arrow_color) :
		fill(Makie.to_color(arrow_color[1]), n)
	else
		fill(Makie.to_color(arrow_color), n)
	end


	directions = if normalize_vecs
		[v / (norm(v) + 1e-14) * scale for v in vecs]
	else
		[v * scale for v in vecs]
	end


	# Proportional arrowhead: tipwidth and tiplength scale with shaft width.
	tip_w = arrow_width * arrowhead_scale * 0.4f0
	tip_l = arrow_width * arrowhead_scale * 0.6f0

	if is3d || N == 3
		ps = Makie.Point3f.([anchors[i][1] for i in 1:n],
			[anchors[i][2] for i in 1:n],
			[anchors[i][3] for i in 1:n])
		ds = Makie.Point3f.([directions[i][1] for i in 1:n],
			[directions[i][2] for i in 1:n],
			[directions[i][3] for i in 1:n])
		Makie.arrows3d!(ax, ps, ds; color = draw_cols, shaftwidth = arrow_width,
		                tipwidth = tip_w, tiplength = tip_l, lengthscale = lengthscale)
	else
		ps = Makie.Point2f.([anchors[i][1] for i in 1:n],
			[anchors[i][2] for i in 1:n])
		ds = Makie.Point2f.([directions[i][1] for i in 1:n],
			[directions[i][2] for i in 1:n])
		Makie.arrows2d!(ax, ps, ds; color = draw_cols, shaftwidth = arrow_width,
		                tipwidth = tip_w, tiplength = tip_l, lengthscale = lengthscale)
	end

	if anchor_visible && anchor_marker !== :none && anchor_marker !== nothing
		acol = anchor_color === nothing ? draw_cols : fill(Makie.to_color(anchor_color), n)
		if N == 3 || is3d
			aps = Makie.Point3f.([anchors[i][1] for i in 1:n],
				[anchors[i][2] for i in 1:n],
				[anchors[i][3] for i in 1:n])
			Makie.scatter!(ax, aps; color = acol, markersize = anchor_size, marker = anchor_marker)
		else
			aps = Makie.Point2f.([anchors[i][1] for i in 1:n],
				[anchors[i][2] for i in 1:n])
			Makie.scatter!(ax, aps; color = acol, markersize = anchor_size, marker = anchor_marker)
		end
	end

	if labels != false
		label_strs = if labels == true || labels == :index
			[string(i) for i in 1:n]
		elseif labels == :latex
			[Makie.L"\mathbf{v}_{%$i}" for i in 1:n]
		elseif labels isa Function
			[string(labels(i)) for i in 1:n]
		elseif labels isa AbstractVector
			string.(labels)
		else
			[string(i) for i in 1:n]
		end

		for i in 1:n
			tip = anchors[i] + directions[i]
			d = normalize_vecs ? vecs[i]/(norm(vecs[i])+1e-14) : vecs[i]
			dn = norm(d)
			perp = dn > 1e-8 ? d/dn : typeof(d)(0, 1)
			sp = label_sp
			offset = if label_pos == :N
				;
				(0.0, sp)
			elseif label_pos == :S
				;
				(0.0, -sp)
			elseif label_pos == :E
				;
				(sp, 0.0)
			elseif label_pos == :W
				;
				(-sp, 0.0)
			else
				(perp[1]*sp, perp[2]*sp)
			end
			if N == 2
				Makie.text!(ax, label_strs[i];
					position = Makie.Point2f(tip[1]+offset[1], tip[2]+offset[2]), color = label_color, fontsize = label_fs,
					align = (:center, :bottom))
			else
				Makie.text!(ax, label_strs[i];
					position = Makie.Point3f(tip[1]+offset[1], tip[2]+offset[2], tip[3]), color = label_color, fontsize = label_fs,
					align = (:center, :bottom))
			end
		end
	end
end


function _resolve_anchors(anchors, vecs)
	n = length(vecs)
	N = length(vecs[1])
	T = eltype(vecs[1])
	zero_v = zero(typeof(vecs[1]))

	anchors === nothing && return fill(zero_v, n)

	if anchors isa AbstractDiscreteCurve
		return [anchors.points[mod1(i, nvertices(anchors))] for i in 1:n]
	elseif anchors isa AbstractVector && length(anchors) == n
		return collect(anchors)
	elseif anchors isa AbstractVector && length(anchors) != n
		throw(DimensionMismatch("anchors length $(length(anchors)) ≠ n=$n"))
	elseif anchors isa Function
		return [anchors(vecs[i]) for i in 1:n]
	else

		a = typeof(vecs[1])(anchors)
		return fill(a, n)
	end
end

_normalise_vecs(v::AbstractVector) = v
_normalise_vecs(v) = [v]  # single vector

_vecs_is3d(v::AbstractVector) = !isempty(v) && length(v[1]) == 3
_vecs_is3d(v) = length(v) == 3
function _subsample(v::AbstractVector, n::Int)
	length(v) ≤ n && return v
	idxs = round.(Int, range(1, length(v); length = n))
	v[idxs]
end
end