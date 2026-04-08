module Curvature
    using LinearAlgebra
    using StaticArrays
    using ..AbstractTypes, ..CurveTypes, ..CurveTopology, ..Iterators

    export curvature, curvatures, curvature_vector, turning_angle, turning_angles, CustomCurvature, TurningAngle, OsculatingCircle, TangentCurvature, SteinerCurvature, LengthWeightedCurvature, SignedCurvature2D

"""κ = θ / ℓ̄  — turning angle divided by mean half-length."""
struct TurningAngle end

"""κ = 1/R  — inverse circumradius of the osculating circle."""
struct OsculatingCircle end

"""
κ^C: `κᵢ = 2tan(|θᵢ|/2) / ℓ̄ᵢ`
 
From Steiner parallel-offset formula with edge-extension gaps (eq 5, option C).
Preserves the **Kirchhoff elastic-curve analogy**.
Warning: blows up near straight-angle vertices (|θ| → π), which is physical
(κ → ∞ at a cusp).
"""
struct TangentCurvature end
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
@inline function _turning_angle_at(p::SVector{2,T}, q::SVector{2,T},
                                    r::SVector{2,T}) where T
    t1   = q - p
    t2   = r - q
    n1sq = t1[1]*t1[1] + t1[2]*t1[2] 
    n2sq = t2[1]*t2[1] + t2[2]*t2[2]
    (n1sq ≤ eps(T)^2 || n2sq ≤ eps(T)^2) && return zero(T)
    atan(t1[1]*t2[2] - t1[2]*t2[1],
         t1[1]*t2[1] + t1[2]*t2[2])
end
@inline function _turning_angle_at(p::SVector{3,T}, q::SVector{3,T},
                                    r::SVector{3,T}) where T
    t1   = q - p
    t2   = r - q
    n1sq = dot(t1, t1)
    n2sq = dot(t2, t2)
    (n1sq ≤ eps(T)^2 || n2sq ≤ eps(T)^2) && return zero(T)
    atan(norm(t1 × t2), dot(t1, t2))
end

@inline function _turning_angle_at(p::SVector{N,T}, q::SVector{N,T},
                                    r::SVector{N,T}) where {N,T}
    t1 = q - p
    t2 = r - q
    a2 = dot(t1, t1)
    b2 = dot(t2, t2)
    (a2 ≤ eps(T)^2 || b2 ≤ eps(T)^2) && return zero(T)
    ab   = dot(t1, t2)
    gram = sqrt(max(zero(T), a2*b2 - ab^2))
    atan(gram, ab)
end

"""
    turning_angle(c, i; signed=false)  →  T
 
Turning angle at vertex `i`.
- `signed=false` (default): unsigned angle in [0, π].
- `signed=true` (2D only): signed angle in (−π, π]. Positive = CCW.
 
```julia
θ   = turning_angle(c, 5)
θ_s = turning_angle(c, 5; signed=true)
```
"""
function turning_angle(c::AbstractDiscreteCurve, i::Int; signed::Bool=false)
    θ = _turning_angle_at(vertex(c, i-1), vertex(c, i), vertex(c, i+1))
    signed ? θ : abs(θ)
end

"""
    turning_angles(c; signed=false)  →  Vector{T}
 
All turning angles. Length = nvertices for closed curves, nvertices-2 for open curves.
"""
function turning_angles(c::AbstractDiscreteCurve{N,T}; signed::Bool=false) where {N,T}
    n_win = length(edge_windows(c; k=3))
    out   = Vector{T}(undef, n_win)
    if signed
        @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
            out[k] = _turning_angle_at(prev, curr, next)
        end
    else
        @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
            out[k] = abs(_turning_angle_at(prev, curr, next))
        end
    end
    out
end
# ─────────────────────────────────────────────────────────────────────────────
# Single-vertex curvature
# ─────────────────────────────────────────────────────────────────────────────
"""
    curvature(c, i, model=TurningAngle())  →  T
 
Discrete curvature at vertex `i` using the specified model.
 
```julia
curvature(c, 5)                            # κ^A: TurningAngle (default)
curvature(c, 5, SteinerCurvature())        # κ^B
curvature(c, 5, TangentCurvature())        # κ^C
curvature(c, 5, OsculatingCircle())        # κ^D = 1/R
curvature(c, 5, SignedCurvature2D())       # signed (2D only)
curvature(c, 5, LengthWeightedCurvature()) # raw |θ|
curvature(c, 5, CustomCurvature(f))        # user kernel
```
"""
@inline function curvature(c::AbstractDiscreteCurve, i::Int, model=TurningAngle())
    _curvature_at(vertex(c, i-1), vertex(c, i), vertex(c, i+1), model)
end

"""
    curvatures(c, model=TurningAngle())  →  Vector{T}
 
All curvature values in one pass.
Length = nvertices (closed) or nvertices-2 (open). 
```julia
curvatures(c)                              # κ^A
curvatures(c, SteinerCurvature())         # κ^B
curvatures(c, TangentCurvature())         # κ^C
curvatures(c, OsculatingCircle())         # κ^D = 1/R
curvatures(c, SignedCurvature2D())        # signed (2D only)
curvatures(c, CustomCurvature(my_f))      # custom kernel
```
"""
function curvatures(c::AbstractDiscreteCurve{N,T}, model=TurningAngle()) where {N,T}
    n_win = length(edge_windows(c; k=3))
    out   = Vector{T}(undef, n_win)
    @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
        out[k] = _curvature_at(prev, curr, next, model)
    end
    out
end

"""
    curvature_vector(c, i)  →  SVector{N,T}
 
Discrete curvature VECTOR at vertex `i`:  `κ · N`.
Points toward the centre of curvature. Magnitude = κ.
 
```julia
v = curvature_vector(c, 5)
κ = norm(v)               # κ (unsigned) at vertex 5
N   = v / norm(v)           # principal normal direction
```
"""
function curvature_vector(c::AbstractDiscreteCurve, i::Int, model=TurningAngle())
    _kn_B(vertex(c, i-1), vertex(c, i), vertex(c, i+1), model)
end

"""
    curvature_vectors(c)  →  Vector{SVector{N,T}}
 
All κ · N vectors.
"""
function curvature_vectors(c::AbstractDiscreteCurve{N,T}, model=TurningAngle()) where {N,T}
    n_win = length(edge_windows(c; k=3))
    out   = Vector{SVector{N,T}}(undef, n_win)
    @inbounds for (k, (prev, curr, next)) in enumerate(edge_windows(c; k=3))
        out[k] = _kn_B(prev, curr, next, model)
    end
    out
end
@inline _dual_length(p, q, r) = (norm(q - p) + norm(r - q)) / 2

@inline function _safe_norm(v::SVector{N,T}) where {N,T}
    n = norm(v)
    n ≤ eps(T) && return zero(T)
    n
end

# canonical unit-difference vector (unnormalised)
@inline function _tangent_diff(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}) where {N,T}
    d1 = q - p; d2 = r - q
    n1 = norm(d1); n2 = norm(d2)
    (n1 < eps(T) || n2 < eps(T)) && return zero(SVector{N,T})
    d1/n1 - d2/n2
end
@inline function _kn_B(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}, ::SteinerCurvature) where {N,T}
    _tangent_diff(p,q,r)
end

@inline function _kn_B(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}, ::TurningAngle) where {N,T}
    v = _tangent_diff(p,q,r)
    mag = _safe_norm(v)
    (mag == zero(T)) && return zero(SVector{N,T})
    θ = abs(_turning_angle_at(p,q,r))
    ℓ̄ = _dual_length(p,q,r)
    κ = θ / ℓ̄
    return v * (κ / mag)
end

@inline function _kn_B(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}, ::TangentCurvature) where {N,T}
    v = _tangent_diff(p,q,r); mag = _safe_norm(v)
    (mag == zero(T)) && return zero(SVector{N,T})
    θ = abs(_turning_angle_at(p,q,r)); ℓ̄ = _dual_length(p,q,r)
    κ = 2 * tan(θ/2) / ℓ̄
    v * (κ / mag)
end

@inline function _kn_B(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}, ::OsculatingCircle) where {N,T}
    d = r - p
    w = norm(d)
    (w ≤ eps(T)) && return zero(SVector{N,T})
    v = _tangent_diff(p,q,r); mag = _safe_norm(v)
    (mag == zero(T)) && return zero(SVector{N,T})
    θ = abs(_turning_angle_at(p,q,r))
    κ = 2 * sin(θ) / w
    v * (κ / mag)
end

@inline function _kn_B(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}, ::LengthWeightedCurvature) where {N,T}
    v = _tangent_diff(p,q,r); mag = _safe_norm(v)
    (mag == zero(T)) && return zero(SVector{N,T})
    θ = abs(_turning_angle_at(p,q,r))
    κ = θ
    v * (κ / mag)
end

@inline function _kn_B(p::SVector{2,T}, q::SVector{2,T}, r::SVector{2,T}, ::SignedCurvature2D) where {T}
    v = _tangent_diff(p,q,r); mag = _safe_norm(v)
    (mag == zero(T)) && return zero(SVector{2,T})
    θ = _turning_angle_at(p,q,r)
    ℓ̄ = _dual_length(p,q,r)
    κ = θ / ℓ̄
    v * (κ / mag)
end

@inline function _kn_B(p::SVector{N,T}, q::SVector{N,T}, r::SVector{N,T}, model::CustomCurvature) where {N,T}
    v = _tangent_diff(p,q,r); mag = _safe_norm(v)
    (mag == zero(T)) && return zero(SVector{N,T})
    κ = model.f(p,q,r)
    v * (κ / mag)
end

@inline function _curvature_at(p, q, r, ::TurningAngle)
    T  = eltype(p)
    θ  = abs(_turning_angle_at(p, q, r))
    ℓ  = _dual_length(p, q, r)
    ℓ > eps(T) ? θ / ℓ : zero(T)
end

@inline function _curvature_at(p, q, r, ::SteinerCurvature)
    T  = eltype(p)
    θ  = abs(_turning_angle_at(p, q, r))
    ℓ  = _dual_length(p, q, r)
    ℓ > eps(T) ? 2sin(θ/2) / ℓ : zero(T)
end

@inline function _curvature_at(p, q, r, ::TangentCurvature)
    T  = eltype(p)
    θ  = abs(_turning_angle_at(p, q, r))
    ℓ  = _dual_length(p, q, r)
    ℓ > eps(T) ? 2tan(θ/2) / ℓ : zero(T)
end
 
# ── κ^D — OsculatingCircle ────────────────────────────────────────────────────
@inline function _curvature_at(p, q, r, ::OsculatingCircle)
    T  = eltype(p)
    θ  = abs(_turning_angle_at(p, q, r))
    w  = norm(r - p)
    w > eps(T) ? 2sin(θ) / w : zero(T)
end
 
# ── SignedCurvature2D ─────────────────────────────────────────────────────────
# Signed version of κ^A. Only defined for N=2; errors for other dimensions.
@inline function _curvature_at(p::SVector{2,T}, q::SVector{2,T}, r::SVector{2,T},
                                ::SignedCurvature2D) where T
    θ = _turning_angle_at(p, q, r)
    ℓ = _dual_length(p, q, r)
    ℓ > eps(T) ? θ / ℓ : zero(T)
end
 
function _curvature_at(::SVector{N,T}, _, _, ::SignedCurvature2D) where {N,T}
    error("SignedCurvature2D requires a 2D curve (N=2), got N=$N")
end
 
@inline _curvature_at(p, q, r, ::LengthWeightedCurvature) =
    abs(_turning_angle_at(p, q, r))
 
@inline _curvature_at(p, q, r, m::CustomCurvature) = m.f(p, q, r)
end