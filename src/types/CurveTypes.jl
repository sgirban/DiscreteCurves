module CurveTypes

using StaticArrays
using ..AbstractTypes

export DiscreteCurve, ClosedCurve, OpenCurve, ParametricCurve
# ── DiscreteCurve{N,T} ─────────────────────────────────────────────────────────
#
# The primary concrete curve type.  Parametric on:
#   N – ambient dimension (compile-time constant, so sizeof(SVector{N,T}) is known)
#   T – scalar type (Float64 default; Float32 for GPU; BigFloat for exact arithmetic)
#

struct DiscreteCurve{N, T <: AbstractFloat} <: AbstractDiscreteCurve{N,T}
    points::Vector{SVector{N, T}}
    closed::Bool

    function DiscreteCurve{N,T}(points::Vector{SVector{N, T}}, closed::Bool=false) where {N, T}
        N ≥ 1 || throw(ArgumentError("Ambient dimension N must be at least 1, but got N=$N."))
        T <: AbstractFloat || throw(ArgumentError("Scalar type T must be a subtype of AbstractFloat, but got T=$(T)."))
        new{N, T}(points, closed)
    end
end

AbstractTypes.topology_trait(c::DiscreteCurve) = c.closed ? ClosedTopology() : OpenTopology()

# ── User-facing constructors  ───────────────────────────────
#
# We offer multiple forms so the user can be as terse as they like.
# Type inference fills in N and T automatically in all cases.

"""
    DiscreteCurve(points; closed=false)
    DiscreteCurve(points, closed)
 
Construct a discrete curve from a collection of points.
 
# Examples
```julia
# 2D open curve from a list of tuples
c = DiscreteCurve([(0.0,0.0),(1.0,0.0),(1.0,1.0)])
 
# 3D closed curve from SVectors
pts = [SVector(cos(θ), sin(θ), 0.0) for θ in range(0, 2π, length=64)[1:end-1]]
c = DiscreteCurve(pts; closed=true)
 
# Using the ClosedCurve / OpenCurve aliases
c = ClosedCurve(pts)
 
# From a matrix (columns = points)
M = rand(3, 100)
c = DiscreteCurve(M)
```
"""
function DiscreteCurve(points::Vector{SVector{N, T}}; closed::Bool=false) where {N, T}
    DiscreteCurve{N, T}(points, closed)
end
function DiscreteCurve(pts::Vector{SVector{N,T}}, closed::Bool) where {N,T}
    DiscreteCurve{N,T}(pts, closed)
end

function DiscreteCurve(pts; closed::Bool=false)
    svecs = _to_svectors(pts)
    DiscreteCurve(svecs; closed)
end
function DiscreteCurve(M::AbstractMatrix{T}; closed::Bool=false) where T<:AbstractFloat
    N, n = size(M)
    pts  = [SVector{N,T}(view(M, :, i)) for i in 1:n]
    DiscreteCurve{N,T}(pts, closed)
end

"""
    ClosedCurve(pts)
 
Shorthand for `DiscreteCurve(pts; closed=true)`.  The last vertex implicitly
connects back to the first.
 
```julia
c = ClosedCurve([SVector(cos(θ), sin(θ)) for θ in 0:0.1:2π-0.1])
```
"""
ClosedCurve(pts) = DiscreteCurve(pts; closed=true)
 
"""
    OpenCurve(pts)
 
Shorthand for `DiscreteCurve(pts; closed=false)` (default).
"""
OpenCurve(pts)   = DiscreteCurve(pts; closed=false)

function _to_svectors(pts::AbstractVector{<:AbstractVector{T}}) where T<:AbstractFloat
    N = length(first(pts))
    [SVector{N,T}(p) for p in pts]
end

function _to_svectors(pts::AbstractVector{<:Tuple})
    T = promote_type(map(eltype, pts[1])...)
    N = length(first(pts))
    [SVector{N,T}(p...) for p in pts]
end
 
function _to_svectors(pts::AbstractVector{<:SVector})
    collect(pts)
end

# ── Parametric constructors ────────────────────────────────────────────────────

"""
    ParametricCurve(f, t0, t1, n; closed=false, T=Float64)
 
Sample a parametric curve `f: ℝ → ℝᴺ` at `n` equally spaced parameter values
in `[t0, t1]` (open) or `[t0, t1)` (closed).
 
`f` can return:
  - an `SVector{N,T}` (fastest, type-stable)
  - a `Tuple` of length N (auto-converted)
  - a plain `Vector` of length N (auto-converted)
 
# Examples
```julia
# Unit circle (2D, closed)
c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)
 
# Helix (3D, open)
c = ParametricCurve(t -> SVector(cos(t), sin(t), t), 0, 4π, 256)
 
# Lissajous figure
c = ParametricCurve(t -> SVector(sin(3t), sin(2t)), 0, 2π, 512; closed=true)
 
# Pass component functions separately
c = ParametricCurve((sin, cos), 0, 2π, 128; closed=true)
 
# Fine control: 32-bit precision (useful for GPU)
c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0f0, 2f0*π, 128;
                    closed=true, T=Float32)
```
"""
function ParametricCurve(f, t0::Real, t1::Real, n::Int;
                          closed::Bool=false, T::Type{<:AbstractFloat}=Float64)
    n ≥ 2 || throw(ArgumentError("n must be ≥ 2, got $n"))
 
    ts = closed ?
        range(T(t0), T(t1); length=n+1)[1:end-1] :
        range(T(t0), T(t1); length=n)
 
    first_val = f(first(ts))
    svecs = _param_to_svectors(f, ts, first_val, T)
    DiscreteCurve(svecs; closed)
end

function ParametricCurve(fs::NTuple{N, Function}, t0::Real, t1::Real, n::Int;
                          closed::Bool=false, T::Type{<:AbstractFloat}=Float64) where N
    f = t -> SVector{N,T}(ntuple(d -> T(fs[d](t)), Val(N)))
    ParametricCurve(f, t0, t1, n; closed, T)
end

function ParametricCurve(f, t0::Real, t1::Real;
                          step::Real, closed::Bool=false,
                          T::Type{<:AbstractFloat}=Float64)
    n = round(Int, (t1 - t0) / step)
    ParametricCurve(f, t0, t1, max(2, n); closed, T)
end

function _param_to_svectors(f, ts, ::SVector{N,T}, ::Type{T}) where {N,T}
    [f(t) for t in ts]
end
 
function _param_to_svectors(f, ts, first_val::Tuple, ::Type{T}) where T
    N = length(first_val)
    [SVector{N,T}(f(t)...) for t in ts]
end
 
function _param_to_svectors(f, ts, first_val::AbstractVector, ::Type{T}) where T
    N = length(first_val)
    [SVector{N,T}(f(t)...) for t in ts]
end

DiscreteCurve(f::Function, t0::Real, t1::Real, n::Int;
                       closed::Bool=false,
                       T::Type{<:AbstractFloat}=Float64) =
    ParametricCurve(f, t0, t1, n; closed, T)

end # module CurveTypes
