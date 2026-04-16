module Flows
 
using LinearAlgebra
using StaticArrays
using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators, ..Lengths
using ..Curvature: _turning_angle_at, _kn_B, _curvature_at,
                   TurningAngle, SteinerCurvature, TangentCurvature, OsculatingCircle,
                   SignedCurvature2D, LengthWeightedCurvature, CustomCurvature

"""
    CurvatureFlow(model=SteinerCurvature(); tangential=false)
 
Curvature-driven flow.  Normal velocity: `vᵢ = κᵢ(model) · N̂ᵢ`.
 
**`tangential=false` (default):** pure normal flow. Vertices move only inward/outward.
  - Efficient, no reparametrisation overhead.
  - Edge lengths become unequal during long flows → use `resample_every` to fix.
 
**`tangential=true`:** adds a tangential velocity component that keeps vertices
  equally spaced along arc length.  This is the *DeTurck trick* (a.k.a.
  arc-length-equalising tangential reparametrisation):
  - Adds `vᵢᵀ = (ℓᵢ − ℓᵢ₋₁)/(ℓᵢ + ℓᵢ₋₁) · Tᵢ` where Tᵢ is the unit tangent.
  - Does NOT change the geometric shape of the curve (same orbits, different speed).
  - Makes the scheme better conditioned for long flows (no mesh degeneration).
  - Slightly more expensive: one extra normalisation per vertex.
 
```julia
CurvatureFlow()                               # κ^B, normal only
CurvatureFlow(TurningAngle())                 # κ^A, normal only
CurvatureFlow(SteinerCurvature(); tangential=true)  # κ^B + arc-length reparam
```
"""
struct CurvatureFlow{M}
    model      :: M
    tangential :: Bool
end
CurvatureFlow(model=SteinerCurvature(); tangential::Bool=false) =
    CurvatureFlow{typeof(model)}(model, tangential)
"""
    PointwiseFlow(f)
 
Flow where `f(curve, i) → SVector{N,T}` gives the velocity at vertex `i`.
`f` receives the full curve, so it may reference any vertex or global quantity.
 
```julia
# Attraction toward centroid
PointwiseFlow((c, i) -> centroid(c) - vertex(c, i))
 
# Repulsion from the origin
PointwiseFlow((c, i) -> vertex(c, i) / norm(vertex(c, i))^2)
```
"""
struct PointwiseFlow{F}
    f :: F   # (AbstractDiscreteCurve, Int) → SVector
end

"""
    BulkFlow(f)
 
Flow where `f(vel, curve) → nothing` fills the entire velocity buffer at once.
Preferred when all velocities are computed in a single coordinated pass.
 
```julia
BulkFlow() do vel, c
    κs = curvature_vectors(c)
    for i in eachindex(κs)
        vel[i] = κs[i]
    end
end
```
"""
struct BulkFlow{F}
    f :: F   # (Vector{SVector}, AbstractDiscreteCurve) → nothing
end

"""
    GradientFlow(f; h=1e-6)
 
Gradient-descent flow minimising scalar functional `f(curve) → T`.
Velocity: `vᵢ = −∇_γᵢ f`  via central finite differences with step `h`.
 
Cost: O(2·N·n) functional evaluations per step. Fine for n ≤ 1000 with cheap f.
For arc-length gradient, use `CurvatureFlow(SteinerCurvature())` instead —
they are theoretically equivalent but CurvatureFlow is exact and O(n).
 
```julia
# Minimise bending energy
GradientFlow(total_squared_curvature)
 
# Minimise Hausdorff distance to a target curve
GradientFlow(c -> hausdorff_distance(c, target))
 
# Minimise enclosed area (collapses the curve)
GradientFlow(enclosed_area; h=1e-5)
```
"""
struct GradientFlow{F, T <: AbstractFloat}
    f :: F
    h :: T
end
GradientFlow(f; h::Real=1e-6) = GradientFlow(f, Float64(h))
GradientFlow(f, h::T) where T <: AbstractFloat = GradientFlow{typeof(f),T}(f, h)


end