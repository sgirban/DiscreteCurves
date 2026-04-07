module Curvature
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators

    export curvature, curvature_vector, turning_angle, turning_angles, CustomCurvature, TurningAngle, OsculatingCircle, SteinerCurvature, LengthWeightedCurvature, SignedCurvature2D

"""κ = θ / ℓ̄  — turning angle divided by mean half-length."""
struct TurningAngle end

"""κ = 1/R  — inverse circumradius of the osculating circle."""
struct OsculatingCircle end

"""κ = 2·sin(θ/2) / ℓ̄  — DEC/Steiner curvature, O(h²) convergence."""
struct SteinerCurvature end

"""κ = θ  — raw turning angle, the integrated curvature per vertex."""
struct LengthWeightedCurvature end

"""κ = signed 2D turning angle / ℓ̄  — positive for left turns. 2D only."""
struct SignedCurvature2D end

"""
    CustomCurvature(f)
 
Wrap a user-provided 3-point kernel `f(prev, curr, next) → T` into the
curvature dispatch system. The library handles topology, iteration, and
pre-allocation. The user only writes the mathematical kernel.
 
```julia
# Define your own curvature: e.g. the projection onto the z-axis normal
my_model = CustomCurvature() do p, q, r
    v1 = q - p
    v2 = r - q
    # … your formula here …
    norm(v2 - v1) / (norm(v1) * norm(v2))
end
 
κs = curvatures(c, my_model)
```
"""
struct CustomCurvature{F}
    f :: F
end
# ─────────────────────────────────────────────────────────────────────────────
# Turning angle
# ─────────────────────────────────────────────────────────────────────────────
 
"""
    turning_angle(c, i; signed=false)  →  T
 
Turning angle at interior vertex `i`.
Default (signed=false): unsigned angle, range [0, π].
If signed=true (2D only): signed angle, range [-π, π]. Positive = CCW.
Uses atan2 internally for maximum numerical precision.
 
```julia
θ = turning_angle(c, 5)           # Unsigned
θ_s = turning_angle(c, 5; signed=true)  # Signed (2D only)
```
"""
function turning_angle(c::AbstractDiscreteCurve, i::Int; signed::Bool=false)
    θ = _turning_angle_at(vertex(c, i-1), vertex(c, i), vertex(c, i+1))
    return signed ? θ : abs(θ)
end

"""
    turning_angles(c; signed=false)  →  Vector{T}
 
All turning angles. Closed curves: length n. Open curves: length n-2 (interior only).
Default (signed=false): unsigned angles, range [0, π].
If signed=true (2D only): signed angles, range [-π, π]. Sum ≈ ±2π for closed curves.
Uses atan2 for maximum numerical precision.
Pre-allocated, @inbounds, zero realloc. O(n).

```julia
θ = turning_angles(c)                    # Unsigned [0, π]
θ_s = turning_angles(c; signed=true)     # Signed [-π, π] (2D only)
```
"""
function turning_angles(c::AbstractDiscreteCurve{N,T}; signed::Bool=false) where {N,T}
    n_win = length(edge_windows(c; k=3))
    out   = Vector{T}(undef, n_win)
    if signed && N == 2
        # 2D: use signed angles directly
        @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
            out[k] = _turning_angle_at(prev, curr, next)
        end
    else
        # Default: unsigned (apply abs)
        @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
            out[k] = abs(_turning_angle_at(prev, curr, next))
        end
    end
    out
end

@inline function _turning_angle_at(p::SVector{2,T}, q::SVector{2,T}, r::SVector{2,T}) where {T<:AbstractFloat}
    # 2D: return SIGNED angle using atan2 for maximum precision
    t1 = q - p
    t2 = r - q
    n1 = norm(t1)
    n2 = norm(t2)
    (n1 < eps(T) || n2 < eps(T)) && return zero(T)
    # Normalized tangent vectors
    u1 = t1 / n1
    u2 = t2 / n2
    # Signed angle via atan2: atan2(cross_product_z, dot_product)
    # Range: [-π, π]. Positive = CCW turn, Negative = CW turn
    atan(u1[1] * u2[2] - u1[2] * u2[1], dot(u1, u2))
end

@inline function _turning_angle_at(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}) where {N,T<:AbstractFloat}
    # 3D+: return SIGNED angle using atan2(magnitude(cross), dot_product) for stability
    # Taking abs() gives unsigned angle [0, π]
    t1 = q - p
    t2 = r - q
    n1 = norm(t1)
    n2 = norm(t2)
    (n1 < eps(T) || n2 < eps(T)) && return zero(T)
    # Use atan2 for superior numerical stability vs acos
    # Cross product magnitude gives sin(θ) ≥ 0, so result is in [0, π]
    c = dot(t1, t2) / (n1 * n2)
    s = norm(t1 × t2) / (n1 * n2)  # Cross product magnitude for sine
    atan(s, c)  # Numerically stable; result [0, π] (always positive since s ≥ 0)
end

# ─────────────────────────────────────────────────────────────────────────────
# Curvature dispatch system
# ─────────────────────────────────────────────────────────────────────────────

"""
    curvature(c::AbstractDiscreteCurve, model::TurningAngle) → Vector

Turning angle divided by mean half-length: κᵢ = θᵢ / ℓ̄ᵢ
Default: unsigned angles [0, π] for all dimensions.
For signed curvature (2D only), use curvature(c, SignedCurvature2D()) instead.
"""
function curvature(c::AbstractDiscreteCurve{N,T}, ::TurningAngle) where {N,T}
    out = turning_angles(c)
    @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
        ℓ₁ = norm(curr - prev)
        ℓ₂ = norm(next - curr)
        ℓ̄ = (ℓ₁ + ℓ₂) / 2
        out[k] = ℓ̄ > eps(T) ? out[k] / ℓ̄ : zero(T)
    end
    out
end

"""
    curvature(c::AbstractDiscreteCurve{2}, model::SignedCurvature2D) → Vector

Signed curvature κᵢ = θᵢ / ℓ̄ᵢ (2D only).
Positive = counterclockwise turn. Sum ≈ 2π (CCW) or -2π (CW).
Uses high-precision atan2-based computation.
"""
function curvature(c::AbstractDiscreteCurve{2,T}, ::SignedCurvature2D) where {T<:AbstractFloat}
    out = turning_angles(c; signed=true)  # Get signed angles
    @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
        ℓ₁ = norm(curr - prev)
        ℓ₂ = norm(next - curr)
        ℓ̄ = (ℓ₁ + ℓ₂) / 2
        out[k] = ℓ̄ > eps(T) ? out[k] / ℓ̄ : zero(T)
    end
    out
end

"""
    curvature(c::AbstractDiscreteCurve, model::SteinerCurvature) → Vector

Steiner/DEC curvature: κᵢ = 2·sin(|θᵢ|/2) / ℓ̄ᵢ.
O(h²) convergence, numerically stable even for small angles.
Uses abs() for 2D to get unsigned curvature.
"""
function curvature(c::AbstractDiscreteCurve{N,T}, ::SteinerCurvature) where {N,T}
    out = turning_angles(c)
    @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
        ℓ₁ = norm(curr - prev)
        ℓ₂ = norm(next - curr)
        ℓ̄ = (ℓ₁ + ℓ₂) / 2
        out[k] = ℓ̄ > eps(T) ? 2 * sin(abs(out[k]) / 2) / ℓ̄ : zero(T)
    end
    out
end

"""
    curvature(c::AbstractDiscreteCurve, model::LengthWeightedCurvature) → Vector

Raw turning angle (integrated curvature per vertex): κᵢ = θᵢ.
No length weighting. Returns unsigned angles [0, π].
For signed angles (2D only), use curvature(c, SignedCurvature2D()) instead.
"""
function curvature(c::AbstractDiscreteCurve{N,T}, ::LengthWeightedCurvature) where {N,T}
    turning_angles(c; signed=false)
end

"""
    curvature(c::AbstractDiscreteCurve, model::CustomCurvature) → Vector

Apply user-provided 3-point kernel to all interior vertices.
"""
function curvature(c::AbstractDiscreteCurve{N,T}, model::CustomCurvature) where {N,T}
    n_win = length(edge_windows(c; k=3))
    out   = Vector{T}(undef, n_win)
    @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
        out[k] = model.f(prev, curr, next)
    end
    out
end

"""
    total_curvature(c; model=TurningAngle()) → T

Sum of all curvatures (integrated total curvature).
For closed curves with LengthWeightedCurvature, should ≈ 2π (CCW) or -2π (CW).
"""
function total_curvature(c::AbstractDiscreteCurve; model=TurningAngle())
    sum(curvature(c, model))
end

"""
    bending_energy(c; model=TurningAngle()) → T

H = Σ κᵢ² · ℓᵢ — integrated square of curvature.
Measures smoothness; lower values indicate smoother curves.
"""
function bending_energy(c::AbstractDiscreteCurve{N,T}; model=TurningAngle()) where {N,T}
    κs = curvature(c, model)
    ℓs = edge_lengths(c)
    energy = zero(T)
    @inbounds for i in 1:length(κs)
        edge_ℓ = (ℓs[i] + ℓs[i+1]) / 2
        energy += κs[i]^2 * edge_ℓ
    end
    energy
end

end