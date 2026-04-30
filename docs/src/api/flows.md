# Flows

The flow system evolves a discrete curve in time by applying a velocity field to
its vertices. **DiscreteCurves.jl** provides a composable set of flow laws,
explicit Euler time integration, CFL timestep control, convergence tracking, and
snapshot saving.

---

## Architecture

Every flow is a **law** — a struct that implements
`compute_velocities!(vel, c, law)`. Laws can be composed with `ScaledFlow`,
`SumFlow`, `ConstrainedFlow`, and `TangentialCorrection`. The top-level entry
point is [`evolve`](@ref), which returns a [`FlowResult`](@ref).

```
your curve c₀
     │
     ▼
evolve(c₀, law; kwargs...)
     │  ┌──────────────────────────────┐
     │  │ for σ in 1:nsteps            │
     │  │   compute_velocities!(v,c,law) │
     │  │   c ← c + τ·v               │
     │  └──────────────────────────────┘
     ▼
FlowResult  (final curve + history + snapshots)
```

---

## Flow Laws

### CurvatureFlow

```julia
CurvatureFlow(model = SteinerCurvature(); tangential = false)
```

Curvature-driven flow. Each vertex moves with velocity proportional to its
curvature in the normal direction:

```math
\dot{\gamma}_i = \kappa_i \hat{N}_i
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model` | Curvature model | `SteinerCurvature()` | Which discrete curvature formula to use. |
| `tangential` | `Bool` | `false` | Add arc-length-equalising tangential component (DeTurck trick). |

**Tangential reparametrisation** (`tangential=true`):

Adds a component
``v_i^T = \frac{\ell_i - \ell_{i-1}}{\ell_i + \ell_{i-1}} \hat{T}_i``
that keeps vertices equally spaced along arc length. This does **not** change
the geometric curve shape, only the parametrisation. Use for long flows to
prevent mesh degeneration.

```julia
CurvatureFlow()                                          # κ^B, normal only
CurvatureFlow(TurningAngle())                            # κ^A
CurvatureFlow(OsculatingCircle())                        # κ^D = 1/R
CurvatureFlow(SteinerCurvature(); tangential=true)       # κ^B + DeTurck
CurvatureFlow(TurningAngle(); tangential=false)          # explicit false
```

```julia
c = fourier_curve(256; seed=3)

# Minimal usage
r = evolve(c, CurvatureFlow())
r.curve   # → shrunken, rounder curve

# With tracking and snapshots
r = evolve(c, CurvatureFlow(SteinerCurvature(); tangential=true);
    nsteps     = 1000,
    save_every = 50,
    track      = [:L => arc_length, :A => c -> abs(signed_area(c)),
                  :IQ => isoperimetric_quotient])

r[:IQ]    # → approaches 1.0 (converges to circle)
r[:L]     # → shrinks toward 0
r[0]      # initial curve
r[500]    # curve after 500 steps
```

---

### PointwiseFlow

```julia
PointwiseFlow(f)
```

Flow where `f(c, i) → SVector{N,T}` gives the velocity at vertex `i`. The
function receives the **full curve** `c`, so it can reference any global
quantity (centroid, curvature at a different vertex, etc.).

| Parameter | Type | Description |
|-----------|------|-------------|
| `f` | `Function` | `(c, i) → SVector{N,T}`. |

```julia
# Attraction toward a fixed point
target = SVector(0.5, 0.0)
r = evolve(c, PointwiseFlow((c, i) -> target - vertex(c, i)))

# Attraction toward centroid (shrink without favouring any direction)
r = evolve(c, PointwiseFlow((c, i) -> centroid(c) - vertex(c, i)))

# Repulsion from origin
r = evolve(c, PointwiseFlow((c, i) -> vertex(c, i) / norm(vertex(c, i))^2))

# Uniform outward expansion
r = evolve(c, PointwiseFlow((c, i) -> inward_normal(c, i) * (-0.1)))

# Combine with curvature via SumFlow
r = evolve(c, SumFlow(
    CurvatureFlow(),
    ScaledFlow(PointwiseFlow((c, i) -> centroid(c) - vertex(c, i)), 0.05)
))
```

---

### BulkFlow

```julia
BulkFlow(f)
```

Flow where `f(vel, c) → nothing` fills the entire velocity buffer at once.
Preferred when all velocities need to be computed in a single coordinated pass
(e.g., non-local interactions, smoothing kernels).

| Parameter | Type | Description |
|-----------|------|-------------|
| `f` | `Function` | `(vel::Vector{SVector{N,T}}, c) → nothing`. Must fill all entries of `vel`. |

```julia
# Laplacian smoothing
laplacian = BulkFlow() do vel, c
    n = nvertices(c)
    for i in 1:n
        prev = vertex(c, mod1(i-1, n))
        next = vertex(c, mod1(i+1, n))
        vel[i] = (prev + next) / 2 - vertex(c, i)
    end
end
r = evolve(c, laplacian; nsteps=50)

# Curvature vectors via BulkFlow
r = evolve(c, BulkFlow() do vel, c
    kvs = curvature_vectors(c)
    for i in eachindex(kvs)
        vel[i] = kvs[i]
    end
end)
```

---

### GradientFlow

```julia
GradientFlow(f; h = 1e-6)
```

Gradient-descent flow that minimises a scalar functional `f(c) → T`.
Velocity: ``v_i = -\nabla_{\gamma_i} f``, computed via central finite differences.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `f` | `Function` | required | `c → Real`. The functional to minimise. |
| `h` | `Real` | `1e-6` | Finite-difference step. |

**Cost**: O(2·N·n) functional evaluations per step. Suitable for `n ≤ 1000`
with inexpensive `f`. For arc-length minimisation, prefer `CurvatureFlow(SteinerCurvature())`
which is exact and O(n).

```julia
# Minimise bending energy ∫κ² ds
bending_energy(c) = sum(curvatures(c).^2) * arc_length(c) / nvertices(c)
r = evolve(c, GradientFlow(bending_energy); nsteps=200,
    track = [:E => bending_energy])

# Minimise arc length (equivalent to CurvatureFlow for closed curves)
r = evolve(c, GradientFlow(arc_length); nsteps=300)

# Minimise enclosed area (collapses the curve to a point)
r = evolve(c, GradientFlow(c -> abs(signed_area(c))); nsteps=100)

# Custom functional: distance to target centroid
target_ctr = SVector(1.0, 0.0)
r = evolve(c, GradientFlow(c -> norm(centroid(c) - target_ctr); h=1e-5))
```

---

## Flow Compositors

Compositors wrap existing laws to modify or combine their velocity fields.

### ScaledFlow

```julia
ScaledFlow(law, α)
```

Multiply all velocities by scalar `α`. Use to control speed or reverse direction.

```julia
ScaledFlow(CurvatureFlow(), 0.5)    # half-speed curve shortening
ScaledFlow(CurvatureFlow(), -1.0)  # curve expansion (reverse)
ScaledFlow(CurvatureFlow(), 2.0)   # double speed (may be less stable)
```

---

### SumFlow

```julia
SumFlow(law1, law2)
```

Additive superposition: ``v_i = v_{1,i} + v_{2,i}``. Nest multiple `SumFlow`
for more than two laws.

```julia
# Curve shortening + centroid attraction
r = evolve(c, SumFlow(
    CurvatureFlow(),
    ScaledFlow(PointwiseFlow((c,i) -> centroid(c) - vertex(c,i)), 0.1)))

# Three-way sum
r = evolve(c, SumFlow(
    SumFlow(CurvatureFlow(), laplacian_flow),
    ScaledFlow(some_other_law, 0.05)))
```

---

### ConstrainedFlow

```julia
ConstrainedFlow(law, g!)
```

Apply `law`, then project velocities via `g!(vel, c) → nothing`. Use for
hard constraints (fixed endpoints, zero mean drift, prescribed velocity, …).

```julia
# Fix both endpoints of an open curve
pinned = ConstrainedFlow(CurvatureFlow(), (vel, c) -> begin
    vel[1]   = zero(eltype(vel))
    vel[end] = zero(eltype(vel))
end)
r = evolve(open_curve, pinned; nsteps=200)

# Zero-mean constraint (remove rigid-body drift)
zero_drift = ConstrainedFlow(CurvatureFlow(), (vel, c) -> begin
    drift = sum(vel) / length(vel)
    for i in eachindex(vel); vel[i] -= drift; end
end)

# Preserve area (project out normal component of mean velocity)
```

---

### TangentialCorrection

```julia
TangentialCorrection(law)
```

Wrap **any** law with an arc-length-equalising tangential component. Equivalent
to `CurvatureFlow(model; tangential=true)` for curvature laws, but works with
**any** velocity law.

The tangential velocity:
``v_i^T = \tfrac{1}{2}(\ell_i - \ell_{i-1})\hat{T}_i``
moves vertices toward the midpoint of their neighbours without changing the
geometric shape of the curve.

```julia
TangentialCorrection(CurvatureFlow(TurningAngle()))
TangentialCorrection(GradientFlow(arc_length))
TangentialCorrection(BulkFlow(laplacian_fn))
```

---

## Custom Laws

Define your own flow by extending `compute_velocities!`:

```julia
struct MyFlow
    speed :: Float64
end

function DiscreteCurves.compute_velocities!(vel, c, law::MyFlow)
    n = nvertices(c)
    for i in 1:n
        # example: move each vertex toward the origin
        vel[i] = -law.speed * vertex(c, i)
    end
end

r = evolve(c, MyFlow(0.1); nsteps=100)
```

The custom law integrates automatically with `evolve`, `cfl_timestep`,
`ConstrainedFlow`, `ScaledFlow`, `SumFlow`, and `TangentialCorrection`.

---

## cfl\_timestep

```julia
cfl_timestep(c; safety = 0.25) → T
```

CFL-stable timestep for curvature flows:

```math
\tau = \text{safety} \times h_{\min}^2
```

where ``h_{\min}`` is the minimum edge length.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `c` | required | The initial curve. |
| `safety` | `0.25` | Safety factor. `0.5` is the theoretical stability boundary; keep below. |

```julia
c = fourier_curve(256)
τ = cfl_timestep(c)                    # 0.25 × min_edge²
τ = cfl_timestep(c; safety=0.1)        # more conservative
τ = cfl_timestep(c; safety=0.5)        # aggressive (use with care)

# For very fine meshes:
c_fine = fourier_curve(2048)
τ_fine = cfl_timestep(c_fine)   # ≈ 16× smaller than 512-vertex curve
```

---

## evolve

```julia
evolve(c₀, law; kwargs...) → FlowResult
```

Evolve curve `c₀` under velocity `law` using explicit Euler time stepping.
This is the main entry point for all curve flows.

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `c₀` | `AbstractDiscreteCurve{N,T}` | required | Initial curve. |
| `law` | Any flow law | required | Velocity law to apply. |

### Keyword Arguments

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `τ` | `Real` or `nothing` | auto CFL | Fixed timestep. `nothing` = auto-select via CFL. |
| `cfl_safety` | `Real` | `0.25` | CFL safety factor (used when `τ=nothing`). |
| `nsteps` | `Int` | `200` | Maximum number of Euler steps σ. |
| `tol` | `Real` | `eps^(1/3)` | Stop early when `max‖vᵢ‖ < tol`. |
| `resample_every` | `Int` | `0` | Resample to uniform arc-length every `k` steps. `0` = never. |
| `adapt_τ_every` | `Int` | `50` | Reduce τ if min edge shrinks. `0` = never. |
| `save_every` | `Int` | `0` | Save intermediate curve every `k` steps. `0` = never. |
| `track` | `Vector` | `[]` | Named functionals: `[:name => function, ...]`. |
| `track_every` | `Int` | `1` | Evaluate tracked functionals every `k` steps. |
| `callback` | `Function` or `nothing` | `nothing` | `f(σ, c, max_v) → Bool`. Return `true` to stop early. |

### Returns

[`FlowResult{N,T}`](@ref).

### Examples

```julia
c = fourier_curve(256; seed=3)

# ── Minimal ────────────────────────────────────────────────────────────────
r = evolve(c, CurvatureFlow())
r.curve         # final curve after 200 steps

# ── With snapshots and tracking ───────────────────────────────────────────
r = evolve(c, CurvatureFlow();
    nsteps     = 500,
    save_every = 10,
    track      = [
        :L   => arc_length,
        :A   => c -> abs(signed_area(c)),
        :IQ  => isoperimetric_quotient,
        :κmax => c -> maximum(curvatures(c)),
    ])

r[0]       # initial curve (requires save_every > 0)
r[10]      # curve after 10 steps
r[100]     # curve after 100 steps
r[:L]      # Vector of arc lengths at each tracked step
r[:IQ]     # Vector of isoperimetric quotients

# ── Fixed timestep ────────────────────────────────────────────────────────
r = evolve(c, CurvatureFlow(); τ=1e-4, nsteps=1000)

# ── Early stopping via callback ───────────────────────────────────────────
r = evolve(c, CurvatureFlow();
    nsteps   = 10000,
    callback = (σ, c, max_v) -> begin
        isoperimetric_quotient(c) > 0.99
    end)
println("Converged in $(r.nsteps) steps")

# ── Gradient flow ─────────────────────────────────────────────────────────
bending(c) = sum(curvatures(c).^2) / nvertices(c)
r = evolve(c, GradientFlow(bending); nsteps=300,
    track = [:E => bending])

# ── Composed law: shortening + attraction to target ───────────────────────
target = SVector(0.5, 0.0)
r = evolve(c, SumFlow(
    CurvatureFlow(),
    ScaledFlow(PointwiseFlow((c,i) -> target - vertex(c,i)), 0.05));
    nsteps = 400)

# ── Constrained: fixed endpoints (open curve) ─────────────────────────────
c_open = spiral_curve(2, 256)
r = evolve(c_open, ConstrainedFlow(CurvatureFlow(), (vel, c) -> begin
    vel[1] = vel[end] = zero(SVector{2,Float64})
end); nsteps=200)
```

---

## evolve\_step!

```julia
evolve_step!(pts, vel, c_view, law, τ) → max_v::T
```

Single explicit Euler step. Lower-level than `evolve`; use for maximum control.

| Argument | Type | Description |
|----------|------|-------------|
| `pts` | `Vector{SVector{N,T}}` | Current point buffer (modified in-place). |
| `vel` | `Vector{SVector{N,T}}` | Velocity buffer (overwritten). |
| `c_view` | `DiscreteCurve{N,T}` | View into `pts` (shares data). |
| `law` | Any flow law | Velocity law. |
| `τ` | `T` | Timestep. |

Returns `max_v::T` — maximum vertex speed (for convergence monitoring).

```julia
c   = fourier_curve(256; seed=5)
N, T = 2, Float64
pts = copy(c.points)
vel = similar(pts)
c_view = DiscreteCurve{N,T}(pts, isclosed(c))
τ = cfl_timestep(c)
law = CurvatureFlow()

for σ in 1:10_000
    max_v = evolve_step!(pts, vel, c_view, law, τ)
    σ % 100 == 0 && println("σ=$σ  max_v=$max_v")
    max_v < 1e-9 && (println("Converged at σ=$σ"); break)
end

final = DiscreteCurve(copy(pts); closed=isclosed(c))
```

---

## FlowResult

```julia
FlowResult{N, T}
```

Returned by [`evolve`](@ref). Carries the final curve, time history of tracked
functionals, and optional intermediate curve snapshots.

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `curve` | `AbstractDiscreteCurve{N,T}` | Final curve after all steps. |
| `nsteps` | `Int` | Actual steps taken (≤ requested `nsteps`). |
| `time` | `T` | Total physical time: ``\sum_\sigma \tau_\sigma``. |
| `max_v` | `Vector{T}` | Maximum vertex speed at each step (auto-tracked). |
| `curves` | `Vector{...}` | Saved intermediate curves (empty if `save_every=0`). |
| `times` | `Vector{T}` | Physical times of each saved snapshot. |

### Indexing

| Syntax | Returns |
|--------|---------|
| `r[:sym]` | `Vector` of tracked functional `sym` values |
| `r[σ]` | `DiscreteCurve` at step `σ` (requires `save_every > 0`) |
| `r[σ => :sym]` | Tracked value of `:sym` at step `σ` |
| `r.curve` | Final curve |
| `r.time` | Total elapsed physical time |
| `r.nsteps` | Number of steps completed |

### Examples

```julia
c = fourier_curve(256; seed=3)
r = evolve(c, CurvatureFlow();
    nsteps    = 300,
    save_every = 10,
    track     = [:L => arc_length, :A => c -> abs(signed_area(c))])

# Access final state
r.curve          # final DiscreteCurve
r.nsteps         # → 300 (or less if converged early)
r.time           # → total physical flow time

# Access snapshots
r[0]             # initial curve
r[10]            # curve at step 10
r[300]           # same as r.curve

# Access tracked history
r[:L]            # → Vector{Float64} of arc lengths
r[:A]            # → Vector{Float64} of areas
r.max_v          # → Vector{Float64} of max velocities

# Convergence diagnostics
using Plots
plot(r[:L]; xlabel="step", ylabel="arc length", label="L")
plot!(r[:A]; label="area")
plot(r.max_v; yscale=:log10, xlabel="step", ylabel="max velocity")

# Check convergence
r.max_v[end] < 1e-6   # → true if flow converged
```

---

## curve\_shortening\_flow

```julia
curve_shortening_flow(c, model = SteinerCurvature(); tangential = false, kwargs...)
```

Convenience wrapper around `evolve(c, CurvatureFlow(model; tangential); kwargs...)`.

```julia
r = curve_shortening_flow(c)                              # κ^B, normal only
r = curve_shortening_flow(c, TurningAngle())              # κ^A
r = curve_shortening_flow(c; tangential=true, nsteps=2000) # with DeTurck
r = curve_shortening_flow(c; save_every=20, track=[:L=>arc_length])
```