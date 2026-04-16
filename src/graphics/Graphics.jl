module Graphics
using ..AbstractTypes, ..CurveTypes
export CurveSkin, CurveFigure, VectorLayer
export curveplot, curveplot!, plotcurve
export vectorplot, vectorplot!
export save_figure
export @curveplot, @vectorplot
@enum CurveSkin begin
	SKIN_STUDIO     = 1
	SKIN_DESIGNER   = 2
	SKIN_GEOMETRIC  = 3
	SKIN_SCIENTIFIC = 4
	SKIN_NEON       = 5
	SKIN_PAPER      = 6
	SKIN_PASTEL     = 7
	SKIN_DARK       = 8
	SKIN_WATERCOLOR = 9
end

const SKIN_NAMES = Dict(
	:studio => SKIN_STUDIO,
	:designer => SKIN_DESIGNER,
	:geometric => SKIN_GEOMETRIC,
	:scientific => SKIN_SCIENTIFIC,
	:neon => SKIN_NEON,
	:paper => SKIN_PAPER,
	:pastel => SKIN_PASTEL,
	:dark => SKIN_DARK,
	:watercolor => SKIN_WATERCOLOR,
)
function _parse_skin(skin)
	skin isa CurveSkin && return skin
	haskey(SKIN_NAMES, skin) && return SKIN_NAMES[skin]
	throw(ArgumentError("Invalid skin: $skin. Valid options are: $(keys(SKIN_NAMES))"))
end

const BUILTIN_VERTEX_METRICS = Set([
	:curvature,
	:arc_length,
	:arc_length_norm,
	:turning_angle,
	:normal_x,
	:normal_y,
	:vertex_index,
	:radius,
	:gradient])
# ─────────────────────────────────────────────────────────────────────────────
#  CurveFigure
# ─────────────────────────────────────────────────────────────────────────────

"""
	CurveFigure
 
Returned by `curveplot`. Wraps a Makie Figure transparently.
 
Supports:
- Direct display in REPL / Jupyter
- `save_figure("output.pdf", fig)` — any Makie-supported format
- `curveplot(c) + vectorplot(vecs)` — overlay vector plots
 
```julia
cf = curveplot(c; skin=:studio)
cf + vectorplot(normals(c); anchors=vertices(c))
save_figure("my_curve.pdf", cf)
```
"""
mutable struct CurveFigure
	_fig::Any
	_ax::Any
	_skin::CurveSkin
	_opt::Any
end

Base.display(cf::CurveFigure) = _require_makie("display(CurveFigure)")
Base.:+(cf::CurveFigure, layer::Any) = _require_makie("CurveFigure +")
function Base.show(io::IO, cf::CurveFigure)
	print(io, "CurveFigure(skin=$(cf._skin), $(cf._ax === nothing ? "no axis" : "ready"))")
end

# ─────────────────────────────────────────────────────────────────────────────
#  VectorLayer
# ─────────────────────────────────────────────────────────────────────────────

"""
	VectorLayer
 
Returned by `vectorplot(...)`. Holds vector data for deferred rendering.
Combine with a `CurveFigure` via `+`:
```julia
curveplot(c) + vectorplot(curvature_vectors(c); anchors=c, color=:curvature)
```
"""
struct VectorLayer
	vectors :: Any
	anchors :: Any
	kwargs  :: Any
end

# ─────────────────────────────────────────────────────────────────────────────
#  curveplot / curveplot!
# ─────────────────────────────────────────────────────────────────────────────

"""
	curveplot(c; skin=:studio, kwargs...)         → CurveFigure
	curveplot([c1,c2,...]; kwargs...)             → CurveFigure
	curveplot(c1, c2, ...; kwargs...)             → CurveFigure
 
Create a beautiful figure from one or more `DiscreteCurve`s.
 
The returned `CurveFigure` supports:
- `+` to overlay vector plots: `curveplot(c) + vectorplot(vecs)`
- `save_figure(filename, fig)` to export to PDF, PNG, SVG, GIF, MP4, …
 
# Quick start
```julia
using GLMakie
c = fourier_curve(256; seed=42)
 
curveplot(c)                                         # studio (default)
curveplot(c; skin=:designer)                         # dramatic dark
curveplot(c; skin=:scientific, vertex_color=:curvature)
curveplot([c1, c2, c3]; skin=:neon)                  # multi-curve
 
# Composition
curveplot(c) + vectorplot(normals(c); anchors=c, vector_scale=0.15)
```
 
# Full option reference
 
## Skin
- `skin` — `:studio` (default), `:designer`, `:neon`, `:paper`, `:geometric`, `:scientific`, `:dark`, `:pastel`, `:watercolor`
 
## Curve line
- `curve_color` — color spec (see below). Skin default if `nothing`.
- `curve_width` — `Real` or `Function(c)→Real`
- `curve_alpha` — opacity [0,1]
- `curve_style` — `:solid` | `:dash` | `:dot`
 
## Glow
- `glow` — `Bool`. Skin default if `nothing`.
- `glow_layers`, `glow_intensity`
 
## Vertices
- `vertex_visible` — `Bool`. Skin default.
- `vertex_color` — color spec. `nothing` = matches curve.
- `vertex_size` — `Real` or `Function(c,i)→Real`
- `vertex_marker` — `:circle` | `:diamond` | `:star5` | `:rect` | `:cross`
- `vertex_alpha`
 
## Vertex labels
- `vertex_labels` — `false` | `true` | `:index` | `:coords` | `:latex` | `Function(i)→String`
- `label_position` — `:N` | `:S` | `:E` | `:W` | `:auto`
- `label_spacing` — distance from vertex (fraction of mean edge length)
- `label_fontsize`
- `label_color` — `:auto` (high-contrast from bg)
 
## Edge labels
- `edge_labels` — `false` | `true` | `Function(c,i)→String` | `Vector{String}`
- `edge_label_position` — `:above` | `:below` | `:left` | `:right`
- `edge_label_spacing`
 
## Hover tooltips (interactive only)
- `tooltip` — `false` | `[:curvature, :arc_length, ...]` | `Function(c,i)→String`
 
## Vector overlay
- `vectors` — `:normals` | `:tangents` | `:curvature_vectors` | `Vector{SVector}`
- `vector_scale`, `vector_width`, `vector_color`, `vector_colormap`
 
## Colors
- `colormap` — any Makie colormap. Skin default if `nothing`.
- `clims` — `(lo,hi)` or `:auto`
 
## Axis / layout
- `axis_visible`, `grid_visible`, `aspect`, `title`, `xlabel`, `ylabel`
- `legend`, `curve_labels` — `Vector{String}` to name each curve
- `size` — figure size in pixels
 
## Color specification
Any of these work for `vertex_color`, `curve_color`, `edge_color`:
- Named color: `:red`, `:royalblue`, `:coral`
- Built-in metric: `:curvature`, `:arc_length`, `:turning_angle`, `:radius`, `:vertex_index`
- `:gradient` — smooth gradient along arc length via colormap
- `:interpolate` — edges interpolate between vertex colors
- `(:curvature, :plasma)` — metric + explicit colormap
- `Function(c, i) → Real` — scalar at vertex i
- `Function(p::SVector) → Real` — scalar from position
- `Vector{Real}` — explicit scalars per vertex
- `Vector{<:Colorant}` — explicit color per vertex
"""
function curveplot(cs; kwargs...)
	_require_makie("curveplot")
end

function curveplot(cs...; kwargs...)
	_require_makie("curveplot")
end

"""
	curveplot!(ax, c; kwargs...)
 
Add curve(s) to an existing Makie `Axis` or `Axis3`.
"""
function curveplot!(args...; kwargs...)
	_require_makie("curveplot!")
end
# ─────────────────────────────────────────────────────────────────────────────
#  vectorplot / vectorplot!
# ─────────────────────────────────────────────────────────────────────────────

"""
	vectorplot(vecs; anchors=nothing, kwargs...)          → VectorLayer or CurveFigure
	vectorplot(vecs, skin=:studio; kwargs...)
 
Draw free vectors as beautiful proportionally-scaled quiver arrows.
 
`vecs` can be:
- A single `SVector{N,T}` — one arrow
- `Vector{SVector{N,T}}` — one arrow per vector
- `Vector{<:AbstractVector}` — auto-converted
 
`anchors` specifies where each arrow starts:
- `nothing` / not given → all anchors at origin `[0,0]`
- `SVector{N,T}` → same anchor for all
- `Vector{SVector}` → one anchor per arrow
- `AbstractDiscreteCurve` → anchors at each vertex (auto-matched by length)
- `Function(v) → SVector` → anchor computed from each vector
 
# Visual options
- `arrow_color` — color spec. `:magnitude` = color by arrow length.
- `arrow_width` — linewidth AND arrowhead size scale proportionally.
- `arrowhead_scale` — arrowhead size = `arrow_width × arrowhead_scale`. Default 3.5.
- `anchor_visible` — show a marker at each anchor. Default `true`.
- `anchor_size` — size of anchor marker. Default proportional to `arrow_width`.
- `anchor_marker` — marker shape. Default `:circle`. Use `:none` to hide.
- `anchor_color` — color of anchor marker. Default matches arrow.
- `normalize` — if `true`, normalise all arrows to unit length (colour by magnitude).
- `scale` — global scale factor on all arrow lengths.
- `colormap`
 
# Labels
- `labels` — `false` | `true` | `:index` | `:latex` | `Function(i)→String`
- `label_position` — `:N` | `:S` | `:E` | `:W` (at arrow tip)
- `label_spacing`
- `label_fontsize`
 
# Examples
```julia
using GLMakie
 
c   = fourier_curve(128; seed=1)
ns  = normals(c)       # Vector{SVector{2,Float64}}
κvs = curvature_vectors(c)
 
# Standalone
vectorplot(ns; anchors=c)
 
# Overlay on curve
curveplot(c) + vectorplot(ns; anchors=c, arrow_color=:royalblue, arrow_width=2.0)
curveplot(c) + vectorplot(κvs; anchors=c, arrow_color=:magnitude, normalize=true)
 
# With labels
curveplot(c) + vectorplot(ns; anchors=c, labels=:latex,
	label_position=:N, label_fontsize=10)
 
# Custom anchors (single anchor = all start from same point)
vectorplot(ns; anchors=SVector(0.0, 0.0))
 
# Free vectors (basis)
e1 = SVector(1.0, 0.0)
e2 = SVector(0.0, 1.0)
vectorplot([e1, e2]; labels=["e₁", "e₂"], anchor_color=:black, arrow_width=2.5)
```
"""
function vectorplot(vecs; kwargs...)
	_require_makie("vectorplot")
end

function vectorplot!(args...; kwargs...)
	_require_makie("vectorplot!")
end
# ─────────────────────────────────────────────────────────────────────────────
#  save_figure
# ─────────────────────────────────────────────────────────────────────────────

"""
	save_figure(filename, fig; kwargs...)
 
Export a `CurveFigure` to file. Format is inferred from the extension.
 
Supported formats (all via Makie/CairoMakie):
- **Vector:** `.pdf`, `.svg`, `.eps`
- **Raster:** `.png`, `.jpeg`, `.tiff`
- **Animation:** `.gif`, `.mp4`, `.webm`
 
Extra kwargs are forwarded to `Makie.save`.
 
```julia
save_figure("curve.pdf", fig)                    # vector, publication-ready
save_figure("curve.png", fig; px_per_unit=2)     # 2× resolution PNG
```
 
Note: PDF/SVG/EPS require CairoMakie to be loaded.
"""
function save_figure(filename::AbstractString, fig::CurveFigure; kwargs...)
	_require_makie("save_figure")
end

# ─────────────────────────────────────────────────────────────────────────────
#  Macros
# ─────────────────────────────────────────────────────────────────────────────

"""
	@curveplot curve [key=val ...]
 
```julia
@curveplot c
@curveplot c skin=:neon vertex_color=:curvature
@curveplot [c1, c2] skin=:paper legend=true
```
"""
macro curveplot(curve_expr, kwargs...)
	kw_exprs = [Expr(:kw, kw.args[1], esc(kw.args[2])) for kw in kwargs
														   if kw isa Expr && kw.head == :(=)]
	quote
		curveplot($(esc(curve_expr)); $(kw_exprs...))
	end
end

"""
	@vectorplot vecs [key=val ...]
 
```julia
@vectorplot normals(c) anchors=c arrow_color=:royalblue
```
"""
macro vectorplot(vecs_expr, kwargs...)
	kw_exprs = [Expr(:kw, kw.args[1], esc(kw.args[2])) for kw in kwargs
														   if kw isa Expr && kw.head == :(=)]
	quote
		vectorplot($(esc(vecs_expr)); $(kw_exprs...))
	end
end

function _require_makie(fn::String)
	@warn """
	DiscreteCurves: `$fn` requires a Makie backend.

	For interactive plots:
		pkg> add GLMakie
		julia> using GLMakie

	For PDF/PNG/SVG file output:
		pkg> add CairoMakie
		julia> using CairoMakie

	For Jupyter/browser:
		pkg> add WGLMakie
		julia> using WGLMakie

	Minimal fallback (using Plots.jl):
		If you only have Plots.jl, `using Plots; plot(c)` will work
		with basic rendering (no skins, no glow, no composition).
	"""
	nothing
end
end
