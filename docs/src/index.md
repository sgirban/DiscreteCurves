# DiscreteCurves.jl

**DiscreteCurves.jl** is a high-performance Julia library for discrete curve geometry.

It is designed for **geometers, applied mathematicians, and computer scientists**
who need correct, fast, and beautiful discrete differential geometry on polylines.

## Philosophy

> *A discrete curve is not merely an approximation of a smooth one —  
> it is a differential-geometric object in its own right.*  
> — Crane & Wardetzky, *A Glimpse into Discrete Differential Geometry*, 2017

This library implements the **mimetic viewpoint**: discrete objects preserve key
structural properties of their smooth counterparts exactly, independent of mesh size.

## Features

- **Four discrete curvature models** from Crane & Wardetzky (κᴬ, κᴮ, κᶜ, κᴰ) + user-extensible protocol
- **Generalised curve flows** — any velocity law, composable, trackable, animated
- **40+ curve generators** — parametric zoo, random polygons, Fourier/Perlin noise
- **Beautiful visualisation** — Makie extension with 8 skins, glow, labels, vector overlays
- **Type-stable, SIMD-friendly kernels** — all hot paths are `@inbounds`, pre-allocated, single-pass
- **Broadcasting and arithmetic operators** — `c + v`, `M * c`, `c .* 2`
- **Named data attachments** — `VertexData`, `EdgeData`, `CurveDataBundle`
- **Arc-length reparametrisation, Hausdorff distance, Frenet frames, writhe, ...**

## Quick Start

```julia
pkg> add DiscreteCurves

using DiscreteCurves
using GLMakie

# Generate a curve
c = fourier_curve(256; modes=8, roughness=0.4)

# Geometric quantities
L  = arc_length(c)
κs = curvatures(c)                  # TurningAngle model (default)
κs = curvatures(c, OsculatingCircle())  # 1/R circumradius model
ns = normals(c)

# Evolve under curve-shortening flow
r = evolve(c, CurvatureFlow(); nsteps=300, save_every=10,
           track=[:L => arc_length])

# Visualise
curveplot(c; skin=:studio)
curveplot(c; vertex_labels=:latex, vertex_visible=true)
curveplot(c) + vectorplot(normals(c); anchors=c)

# Animate the flow
animate_flow(r; filename="flow.gif", fps=30)
flow_slider(r)     # interactive scrubber (GLMakie only)
```

## Installation

```julia
] add DiscreteCurves
```

Optional extensions (loaded automatically when you `using` them):

| Backend | Purpose |
|---------|---------|
| `GLMakie` | Interactive desktop plots |
| `CairoMakie` | PDF/PNG/SVG export |
| `WGLMakie` | Jupyter/browser plots |
| `Plots` | Minimal fallback (basic curves only) |
| `CUDA` | GPU acceleration |

## Contents

```@contents
Pages = ["getting_started.md", "guide/constructors.md", "guide/geometry.md",
         "guide/flows.md", "guide/plotting.md", "mathematics.md", "performance.md"]
Depth = 2
```
