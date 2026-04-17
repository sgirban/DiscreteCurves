# Getting Started

## Installation

```julia
using Pkg
Pkg.add("DiscreteCurves")
```

For visualisation, add one of the Makie backends:

```julia
Pkg.add("GLMakie")    # interactive desktop window
Pkg.add("CairoMakie") # PDF/SVG/PNG file output
```

## Your first curve

```julia
using DiscreteCurves
using GLMakie

# ── Construct ─────────────────────────────────────────────────────────────────

# Parametric: give a function t → (x, y)
c = ParametricCurve(t -> (cos(t), sin(t)), 0, 2π, 128; closed=true)

# Or use the macro
c = @curve (cos(t), sin(t))  t ∈ 0..2π  n=128  closed

# Or use a generator
c = circle_curve(128)
c = fourier_curve(256; modes=6, roughness=0.5, seed=42)
c = random_convex_polygon(16; irregularity=0.4)

# ── Inspect ───────────────────────────────────────────────────────────────────

nvertices(c)      # → 128
nedges(c)         # → 128  (closed: nedges = nvertices)
isclosed(c)       # → true
arc_length(c)     # → ≈ 2π
vertex(c, 5)      # → SVector{2,Float64} — the 5th vertex
vertex(c, 129)    # → same as vertex(c, 1) for closed curves (modular)

# ── Geometry ──────────────────────────────────────────────────────────────────

κs = curvatures(c)                        # Vector of curvatures (κᴬ by default)
κs = curvatures(c, OsculatingCircle())    # 1/R model
ns = normals(c)                           # unit outward normals
ts = unit_tangents(c; at=:vertex)         # unit vertex tangents
θs = tangent_angles(c)                    # atan(t_y, t_x) for each edge

L  = arc_length(c)
A  = enclosed_area(c)
s  = cumulative_length(c)                 # arc-length parameter at each vertex

# ── Transform ─────────────────────────────────────────────────────────────────

c2 = c + SVector(1.0, 0.0)               # translate
c2 = c * 2.0                             # scale
c2 = rotate(c, π/4)                      # rotate 45° CCW
c2 = resample(c, 512)                    # arc-length uniform resampling

# ── Flow ──────────────────────────────────────────────────────────────────────

r = evolve(c, CurvatureFlow();
    nsteps     = 300,
    save_every = 10,           # save a snapshot every 10 steps
    track      = [:L => arc_length, :A => enclosed_area])

r.curve        # final curve
r.time         # total physical time elapsed
r[:L]          # Vector of arc lengths (one per tracked step)
r[0]           # initial curve (σ=0)
r[50]          # curve after 50 steps

# ── Plot ──────────────────────────────────────────────────────────────────────

curveplot(c)                                       # studio skin (default)
curveplot(c; skin=:designer)
curveplot(c; vertex_color=:curvature, colormap=:plasma)
curveplot(c; vertex_labels=:latex, curve_name="\\gamma")
curveplot(c) + vectorplot(normals(c); anchors=c)

plot_flow_evolution(r; nframes=9)
flow_slider(r)                                     # interactive
animate_flow(r; filename="flow.gif", fps=30)
```

## Key concepts

### σ vs τ

- **τ** (tau): the *physical timestep* — how far forward in time each Euler step moves.
  Too large → numerical blow-up. Auto-CFL selects `τ = 0.25 × min_edge²`.
- **σ** (sigma): the *step index* — an integer counting how many Euler steps ran.
  `r[σ]` returns the curve after σ steps. Physical time = Σ τₖ ≈ σ·τ.

### Curvature models

Four models from Crane & Wardetzky (2017):

| Symbol | Julia type | Formula | Flow property |
|--------|-----------|---------|--------------|
| κᴬ | `TurningAngle()` | θᵢ/ℓ̄ᵢ | Preserves total curvature |
| κᴮ | `SteinerCurvature()` | 2sin(θ/2)/ℓ̄ | Preserves centre of mass |
| κᶜ | `TangentCurvature()` | 2tan(θ/2)/ℓ̄ | Kirchhoff analogy |
| κᴰ | `OsculatingCircle()` | 2sin(θ)/w = 1/R | Circles are fixed points |

### Colour specifications

Any plotting kwarg accepting a colour also accepts:

- `:curvature`, `:arc_length`, `:radius`, `:vertex_index` — built-in metrics
- `:gradient` — smooth gradient along arc length
- `(:curvature, :viridis)` — metric + colormap
- `Function(c, i) → Real` — custom per-vertex scalar
- `Vector{Real}` — explicit scalars per vertex
