# Curvature

Discrete curvature is one of the central quantities in differential geometry.
**DiscreteCurves.jl** implements four classical models from Crane & Wardetzky (2017),
a signed 2D variant, a length-weighted variant, and a fully user-extensible
`CustomCurvature` protocol — all sharing the same iterator infrastructure.

---

## Curvature Models

Each model is a lightweight sentinel struct passed as the second argument to
[`curvature`](@ref), [`curvatures`](@ref), [`curvature_vector`](@ref), and
[`curvature_vectors`](@ref).

| Type | Symbol | Formula | Highlight |
|------|--------|---------|-----------|
| `TurningAngle()` | κᴬ | ``\theta_i / \bar\ell_i`` | Exact Gauss-Bonnet: ``\sum \kappa_i^A \bar\ell_i = 2\pi`` |
| `SteinerCurvature()` | κᴮ | ``2\sin(\theta_i/2) / \bar\ell_i`` | ``-\nabla_{\gamma_i} L = \kappa_i^B N_i`` (gradient of arc length) |
| `TangentCurvature()` | κᶜ | ``2\tan(\theta_i/2) / \bar\ell_i`` | Kirchhoff elastic analogy |
| `OsculatingCircle()` | κᴰ | ``2\sin(\theta_i) / w_i = 1/R_i`` | Circles are exact fixed points |
| `LengthWeightedCurvature()` | — | ``|\theta_i|`` | Integrated curvature, no length normalisation |
| `SignedCurvature2D()` | — | ``\theta_i / \bar\ell_i`` (signed) | 2D only; positive = left turn (CCW) |
| `CustomCurvature(f)` | — | ``f(p, q, r)`` | User-defined 3-point kernel |

Here ``\theta_i`` is the turning angle at vertex `i`, ``\bar\ell_i = (\ell_{i-1} + \ell_i)/2``
is the dual edge length, and ``w_i = \|\gamma_{i+1} - \gamma_{i-1}\|``.

**No free lunch:** no single discrete curvature preserves all properties of
smooth curvature simultaneously. Choose the model that preserves the invariant
you care about.

---

## TurningAngle

```julia
TurningAngle()
```

**κᴬ** — turning angle divided by dual edge length.

- Default model throughout the library.
- Satisfies the discrete Gauss-Bonnet theorem: ``\sum_i \kappa_i^A \bar\ell_i = 2\pi``
  for any simple closed curve traversed once.
- Does **not** equal ``-\nabla_{\gamma_i} L`` (use `SteinerCurvature` for that).

---

## SteinerCurvature

```julia
SteinerCurvature()
```

**κᴮ** — Steiner / DEC curvature.

- ``\kappa_i^B = 2\sin(\theta_i/2) / \bar\ell_i``.
- The associated curvature vector ``\kappa_i^B N_i`` equals the exact negative
  gradient of arc length with respect to ``\gamma_i``. This makes it the
  theoretically correct velocity for **curve-shortening flow**.
- The default model in [`CurvatureFlow`](@ref).
- Converges to smooth curvature at **O(h²)** for uniform meshes.

---

## TangentCurvature

```julia
TangentCurvature()
```

**κᶜ** — tangent curvature.

- ``\kappa_i^C = 2\tan(\theta_i/2) / \bar\ell_i``.
- Arises from the **Kirchhoff elastic rod** analogy.
- **Warning**: diverges when ``|\theta_i| \to \pi`` (180° hairpin). Avoid for
  curves with near-cusps.

---

## OsculatingCircle

```julia
OsculatingCircle()
```

**κᴰ** — inverse circumradius of the osculating circle.

- ``\kappa_i^D = 2\sin(\theta_i) / \|\gamma_{i+1} - \gamma_{i-1}\| = 1/R_i``.
- Circles are **exact fixed points** of the osculating curvature: a discrete
  circle has constant ``\kappa^D`` equal to its smooth curvature.
- Best choice when you need accurate point-wise curvature values, not gradient flow.

---

## LengthWeightedCurvature

```julia
LengthWeightedCurvature()
```

Raw turning angle ``|\theta_i|``, with no division by dual length. This is the
**integrated curvature** per vertex — more of a combinatorial quantity than a
pointwise one. Useful when you need the raw angular deficit (e.g. for discrete
Gauss maps).

---

## SignedCurvature2D

```julia
SignedCurvature2D()
```

Signed version of `TurningAngle`, valid only for 2D curves (`N = 2`).
Sign convention: **positive = left turn = CCW direction**.

- For a CCW circle: all values are positive.
- For a CW circle: all values are negative.
- Throws an error if called on a 3D or higher-dimensional curve.

---

## CustomCurvature

```julia
CustomCurvature(f)
```

Wrap a user-supplied 3-point kernel `f(prev, curr, next) → T` into the standard
curvature dispatch system. The library handles topology (open/closed boundary
wrapping), iteration, and pre-allocation. You only write the formula.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `f` | `Function` | A callable with signature `f(prev::SVector, curr::SVector, next::SVector) → Real`. |

### Examples

```julia
# Menger curvature: 1/R of the circumcircle through (p, q, r)
menger = CustomCurvature() do p, q, r
    a = norm(p - q);  b = norm(q - r);  c = norm(r - p)
    area = abs((q-p)[1]*(r-p)[2] - (r-p)[1]*(q-p)[2]) / 2
    area > 1e-14 ? (a*b*c) / (4*area) |> inv : 0.0
end

c  = circle_curve(128)
κs = curvatures(c, menger)   # all ≈ 1.0 for a unit circle

# Discrete mean curvature (cotangent weights, 2D)
cot_model = CustomCurvature() do p, q, r
    t1 = q - p;  t2 = r - q
    n1 = norm(t1);  n2 = norm(t2)
    (n1 < 1e-14 || n2 < 1e-14) && return 0.0
    cosθ = dot(t1,t2)/(n1*n2)
    sinθ = sqrt(max(0.0, 1 - cosθ^2))
    sinθ < 1e-14 ? 0.0 : cosθ/sinθ / ((n1+n2)/2)
end
```

---

## turning\_angle

```julia
turning_angle(c, i; signed = false) → T
```

Turning angle at vertex `i`: the angle between the incoming and outgoing edge
directions.

```math
\theta_i = \angle(\gamma_i - \gamma_{i-1},\; \gamma_{i+1} - \gamma_i)
```

Unsigned (default): result is in ``[0, \pi]``.  
Signed (2D only): result is in ``(-\pi, \pi]``, positive = CCW.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | Any discrete curve. |
| `i` | `Int` | Vertex index. Wraps for closed curves. |
| `signed` | `Bool` | If `true`, return a signed angle (2D only). Default: `false`. |

### Returns

`T` — turning angle in radians.

### Examples

```julia
c = circle_curve(8)

turning_angle(c, 1)              # → π/4 ≈ 0.785 (regular octagon)
turning_angle(c, 1; signed=true) # → +π/4 (CCW circle, positive)

# Right-angle polygon
sq = square_curve(n=4)
turning_angle(sq, 1)    # → π/2

# Gauss-Bonnet check (sum of all turning angles = 2π for simple closed CCW curve)
c  = fourier_curve(256)
θs = turning_angles(c; signed=true)
sum(θs)   # → ≈ 2π
```

### Notes

For open curves, vertex indices 1 and `nvertices(c)` have one adjacent edge
missing. The library uses one-sided differences at these boundary vertices.
The `turning_angles` function (plural) uses `edge_windows` and therefore only
returns values for interior vertices of open curves.

### See Also

[`turning_angles`](@ref), [`curvature`](@ref)

---

## turning\_angles

```julia
turning_angles(c; signed = false) → Vector{T}
```

All turning angles in a single pass.

- **Closed curves**: returns `nvertices(c)` values (all vertices).
- **Open curves**: returns `nvertices(c) - 2` values (interior vertices only;
  boundary vertices have no centred 3-point window).

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | Any discrete curve. |
| `signed` | `Bool` | Return signed angles (2D only). Default: `false`. |

### Returns

`Vector{T}` of turning angles.

### Examples

```julia
c = regular_polygon(6)    # regular hexagon
θs = turning_angles(c)    # all ≈ π/3 ≈ 1.047

# Gauss-Bonnet identity
c  = fourier_curve(512; seed=7)
θs = turning_angles(c; signed=true)
sum(θs)   # → ≈ 2π  (closed, CCW)

c_cw = orient(c; to=CW)
sum(turning_angles(c_cw; signed=true))   # → ≈ -2π

# Histogram of turning angles (detect corners vs smooth regions)
using Statistics
println("Max turning angle: $(maximum(θs)) rad")
println("Mean: $(mean(θs)) rad")
```

### See Also

[`turning_angle`](@ref), [`curvature`](@ref), [`curvatures`](@ref)

---

## curvature

```julia
curvature(c, i, model = TurningAngle()) → T
```

Discrete curvature at a single vertex `i` using the specified model.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve` | The curve. |
| `i` | `Int` | Vertex index. Wraps for closed curves. |
| `model` | Curvature model | One of the model structs above. Default: `TurningAngle()`. |

### Returns

`T` — scalar curvature value (always non-negative, except `SignedCurvature2D`).

### Examples

```julia
c = circle_curve(256)

curvature(c, 1)                        # ≈ 1.0  (κ^A, unit circle)
curvature(c, 1, SteinerCurvature())    # ≈ 1.0  (κ^B)
curvature(c, 1, OsculatingCircle())    # ≈ 1.0  (κ^D = 1/R, exact for circles)
curvature(c, 1, TangentCurvature())    # ≈ 1.0  (κ^C)
curvature(c, 1, SignedCurvature2D())   # ≈ +1.0 (CCW circle → positive)

# CW circle → negative signed curvature
c_cw = orient(circle_curve(256); to=CW)
curvature(c_cw, 1, SignedCurvature2D())   # ≈ -1.0

# Ellipse: curvature varies between a/b² and b/a²
c = ellipse_curve(512; a=2.0, b=1.0)
κ_min = minimum(curvature(c, i) for i in 1:nvertices(c))   # ≈ b/a² = 0.25
κ_max = maximum(curvature(c, i) for i in 1:nvertices(c))   # ≈ a/b² = 2.0
```

### See Also

[`curvatures`](@ref), [`curvature_vector`](@ref), [`turning_angle`](@ref)

---

## curvatures

```julia
curvatures(c, model = TurningAngle()) → Vector{T}
```

All curvature values in one pass. This is significantly faster than calling
`curvature(c, i, model)` in a loop because it fuses the window iteration and
the curvature computation into a single cache-friendly sweep.

- **Closed curves**: length = `nvertices(c)`.
- **Open curves**: length = `nvertices(c) - 2` (interior vertices only).

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | The curve. |
| `model` | Curvature model | Default: `TurningAngle()`. |

### Returns

`Vector{T}` of curvature values.

### Examples

```julia
c = fourier_curve(512; seed=3)

κA = curvatures(c)                          # TurningAngle (default)
κB = curvatures(c, SteinerCurvature())      # Steiner
κC = curvatures(c, TangentCurvature())      # Tangent
κD = curvatures(c, OsculatingCircle())      # Osculating (1/R)
κS = curvatures(c, SignedCurvature2D())     # Signed (2D only)

# Total absolute curvature
sum(κA)          # for closed curves: ≈ 2π / mean_dual_length

# Bending energy E = ∫ κ² ds ≈ ∑ κᵢ² ℓ̄ᵢ
ℓ = edge_lengths(c)
# dual lengths
ℓ_d = [(ℓ[mod1(i-1, nedges(c))] + ℓ[i])/2 for i in 1:nvertices(c)]
E = sum(κA[i]^2 * ℓ_d[i] for i in eachindex(κA))

# Attach as vertex data
vd = attach_vertex_data(c, curvatures(c))
b  = bundle(c)
b[:curvature] = curvatures(c)

# Compare models on a non-circle
c = astroid(512)
println("Model agreement: κ^A vs κ^D difference = ",
        maximum(abs.(curvatures(c) .- curvatures(c, OsculatingCircle()))))
```

### See Also

[`curvature`](@ref), [`curvature_vectors`](@ref), [`attach_vertex_data`](@ref)

---

## curvature\_vector

```julia
curvature_vector(c, i, model = TurningAngle()) → SVector{N,T}
```

Discrete curvature **vector** at vertex `i`: ``\kappa_i \hat{N}_i``.
Points toward the centre of curvature. The magnitude equals the scalar curvature.

For `SteinerCurvature`, the curvature vector is exactly
``(\hat{T}_{\text{out}} - \hat{T}_{\text{in}})``, the negative gradient of
arc length.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{N,T}` | The curve. |
| `i` | `Int` | Vertex index. |
| `model` | Curvature model | Default: `TurningAngle()`. |

### Returns

`SVector{N,T}` — the curvature vector at vertex `i`.

### Examples

```julia
c  = circle_curve(128)
kv = curvature_vector(c, 1)
norm(kv)        # ≈ 1.0  (κ = 1 for unit circle)
kv / norm(kv)   # ≈ -vertex(c, 1)  (points toward origin = centre of curvature)

# Verify: curvature vector is parallel to inward normal
n = inward_normal(c, 1)
dot(normalize(kv), n)   # → ≈ 1.0

# Curvature vector as a flow velocity (curve shortening)
for i in 1:nvertices(c)
    velocity = curvature_vector(c, i, SteinerCurvature())
    # velocity is exactly -∇_{γᵢ} arc_length
end
```

### See Also

[`curvature_vectors`](@ref), [`curvature`](@ref), [`inward_normal`](@ref)

---

## curvature\_vectors

```julia
curvature_vectors(c, model = TurningAngle()) → Vector{SVector{N,T}}
```

All curvature vectors in one pass.

- **Closed curves**: length = `nvertices(c)`.
- **Open curves**: length = `nvertices(c) - 2`.

### Examples

```julia
c   = fourier_curve(256; seed=5)
kvs = curvature_vectors(c)

# Magnitudes = scalar curvatures
κs = norm.(kvs)

# Flow velocity buffer (Steiner model = exact gradient of arc length)
vel = curvature_vectors(c, SteinerCurvature())

# Visualise: plot arrows at each vertex
using GLMakie
curveplot(c) + vectorplot(curvature_vectors(c); anchors=c,
                          arrow_color=:magnitude, normalize=true)

# Verify zero-sum property for closed curves (Steiner model)
sum(curvature_vectors(c, SteinerCurvature()))   # ≈ SVector(0.0, 0.0)
```

### See Also

[`curvature_vector`](@ref), [`curvatures`](@ref), [`CurvatureFlow`](@ref)