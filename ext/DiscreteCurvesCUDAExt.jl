module DiscreteCurvesCUDAExt

using DiscreteCurves
using DiscreteCurves.AbstractTypes: AbstractDiscreteCurve, ClosedTopology, OpenTopology
using DiscreteCurves.CurveTypes: DiscreteCurve
using DiscreteCurves.CurveTopology: isclosed, nvertices, nedges
using DiscreteCurves.Curvature: SteinerCurvature, TurningAngle

using CUDA
using CUDA: CuArray, @cuda, threadIdx, blockIdx, blockDim
using LinearAlgebra
using StaticArrays

"""
    GPUNotSupportedError

Thrown when a `GPUCurve` operation has no GPU implementation.
The error message explains available alternatives.
"""
struct GPUNotSupportedError <: Exception
    msg :: String
end
Base.showerror(io::IO, e::GPUNotSupportedError) =
    print(io, "GPUNotSupportedError: ", e.msg)

export GPUNotSupportedError

struct GPUVertexArray{N, T <: AbstractFloat}
    data :: CuMatrix{T}

    function GPUVertexArray{N,T}(data::CuMatrix{T}) where {N,T}
        size(data, 1) == N || throw(DimensionMismatch(
            "GPUVertexArray{$N}: matrix has $(size(data,1)) rows, expected $N"))
        new{N,T}(data)
    end
end

GPUVertexArray{N}(data::CuMatrix{T}) where {N,T} = GPUVertexArray{N,T}(data)

Base.length(a::GPUVertexArray)                     = size(a.data, 2)
Base.eltype(::GPUVertexArray{N,T}) where {N,T}    = SVector{N,T}
Base.copy(a::GPUVertexArray{N,T}) where {N,T}     = GPUVertexArray{N,T}(copy(a.data))

function Base.getindex(a::GPUVertexArray{N,T}, i::Int) where {N,T}
    col = Array(view(a.data, :, i)) 
    SVector{N,T}(ntuple(d -> col[d], Val(N)))
end

"""
    GPUCurve{N, T} <: AbstractDiscreteCurve{N,T}

A discrete curve with GPU-resident vertices.

Construct from a CPU curve with `to_gpu`. Download with `to_cpu`.
All standard `DiscreteCurves` geometry and flow operations dispatch
to GPU kernels automatically — identical call syntax to CPU curves.

```julia
gc = to_gpu(fourier_curve(4096))

arc_length(gc)                # GPU kernel
curvatures(gc)                # GPU kernel → Vector{Float64}
evolve(gc, CurvatureFlow())   # full GPU loop → FlowResult
```
"""
struct GPUCurve{N, T <: AbstractFloat} <: AbstractDiscreteCurve{N,T}
    points :: GPUVertexArray{N,T}
    closed :: Bool
end

DiscreteCurves.AbstractTypes.topology_trait(g::GPUCurve) =
    g.closed ? ClosedTopology() : OpenTopology()

function Base.show(io::IO, g::GPUCurve{N,T}) where {N,T}
    topo = g.closed ? "closed" : "open"
    mem  = Base.format_bytes(sizeof(T) * N * nvertices(g))
    print(io, "GPUCurve{$(N)D,$(T)}  $topo  $(nvertices(g)) vertices  [$mem on GPU]")
end

export GPUCurve

"""
    to_gpu(c::DiscreteCurve{N,T}) → GPUCurve{N,T}

Upload a CPU curve to the GPU. One H→D transfer. Allocates `N × nvertices` `T`
values on the device.
"""
function to_gpu(c::DiscreteCurve{N,T}) where {N,T}
    n = nvertices(c)
    mat = Matrix{T}(undef, N, n)
    @inbounds for j in 1:n, d in 1:N
        mat[d, j] = c.points[j][d]
    end
    GPUCurve{N,T}(GPUVertexArray{N,T}(CuArray(mat)), c.closed)
end

"""
    to_cpu(g::GPUCurve{N,T}) → DiscreteCurve{N,T}

Download a GPU curve to the CPU.
"""
function to_cpu(g::GPUCurve{N,T}) where {N,T}
    mat = Array(g.points.data)        
    n   = size(mat, 2)
    pts = Vector{SVector{N,T}}(undef, n)
    @inbounds for j in 1:n
        pts[j] = SVector{N,T}(ntuple(d -> mat[d, j], Val(N)))
    end
    DiscreteCurve{N,T}(pts, g.closed)
end

export to_gpu, to_cpu

@inline _mod1_gpu(i::Int32, n::Int32) = ((i - Int32(1)) % n + n) % n + Int32(1)

@inline _adj_gpu(i::Int32, n::Int32, closed::Bool) =
    closed ? _mod1_gpu(i, n) : clamp(i, Int32(1), n)

@inline _gidx() =
    Int32(threadIdx().x) + Int32(blockDim().x) * (Int32(blockIdx().x) - Int32(1))

# Launch a 1-D kernel covering `n` work items
@inline function _launch1d(kernel!, n::Int, args...)
    n == 0 && return
    threads = min(n, 256)
    blocks  = cld(n, threads)
    @cuda threads=threads blocks=blocks kernel!(args...)
end

function _edge_lengths_kernel!(out, pts, n::Int32, N::Int32, closed::Bool)
    i = _gidx()
    ne = closed ? n : n - Int32(1)
    i > ne && return nothing
    T = eltype(pts)
    j = _adj_gpu(i + Int32(1), n, closed)
    sq = zero(T)
    @inbounds for d in Int32(1):N
        dv = pts[d, j] - pts[d, i]
        sq = fma(dv, dv, sq)
    end
    out[i] = sqrt(sq)
    return nothing
end

function _edge_sq_kernel!(out, pts, n::Int32, N::Int32, closed::Bool)
    i = _gidx()
    ne = closed ? n : n - Int32(1)
    i > ne && return nothing
    T = eltype(pts)
    j = _adj_gpu(i + Int32(1), n, closed)
    sq = zero(T)
    @inbounds for d in Int32(1):N
        dv = pts[d, j] - pts[d, i]
        sq = fma(dv, dv, sq)
    end
    out[i] = sq
    return nothing
end

function _area_terms_kernel!(out, pts, n::Int32, closed::Bool)
    i = _gidx()
    ne = closed ? n : n - Int32(1)
    i > ne && return nothing
    T = eltype(pts)
    j = _adj_gpu(i + Int32(1), n, closed)
    out[i] = pts[1, i] * pts[2, j] - pts[1, j] * pts[2, i]
    return nothing
end

function _curvatures_ta_kernel!(out, pts, n::Int32, N::Int32, closed::Bool)
    i = _gidx()
    nw = closed ? n : n - Int32(2)
    i > nw && return nothing
    T  = eltype(pts)
    iv = closed ? i : i + Int32(1)
    ip = _adj_gpu(iv - Int32(1), n, closed)
    iq = _adj_gpu(iv + Int32(1), n, closed)

    n1sq = zero(T); n2sq = zero(T); dotv = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, iv] - pts[d, ip]
        e2 = pts[d, iq] - pts[d, iv]
        n1sq = fma(e1, e1, n1sq)
        n2sq = fma(e2, e2, n2sq)
        dotv = fma(e1, e2, dotv)
    end
    n1 = sqrt(n1sq); n2 = sqrt(n2sq)
    eps_T = eps(T)
    if n1 ≤ eps_T || n2 ≤ eps_T
        out[i] = zero(T); return nothing
    end
    cross_sq = max(zero(T), n1sq * n2sq - dotv * dotv)
    θ        = atan(sqrt(cross_sq), dotv)
    ℓ̄        = (n1 + n2) / 2
    out[i]   = ℓ̄ > eps_T ? abs(θ) / ℓ̄ : zero(T)
    return nothing
end

function _curvatures_steiner_kernel!(out, pts, n::Int32, N::Int32, closed::Bool)
    i = _gidx()
    nw = closed ? n : n - Int32(2)
    i > nw && return nothing
    T  = eltype(pts)
    iv = closed ? i : i + Int32(1)
    ip = _adj_gpu(iv - Int32(1), n, closed)
    iq = _adj_gpu(iv + Int32(1), n, closed)

    n1sq = zero(T); n2sq = zero(T); dotv = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, iv] - pts[d, ip]
        e2 = pts[d, iq] - pts[d, iv]
        n1sq = fma(e1, e1, n1sq)
        n2sq = fma(e2, e2, n2sq)
        dotv = fma(e1, e2, dotv)
    end
    n1 = sqrt(n1sq); n2 = sqrt(n2sq)
    eps_T = eps(T)
    if n1 ≤ eps_T || n2 ≤ eps_T
        out[i] = zero(T); return nothing
    end
    cross_sq = max(zero(T), n1sq * n2sq - dotv * dotv)
    θ        = atan(sqrt(cross_sq), dotv)
    ℓ̄        = (n1 + n2) / 2
    out[i]   = ℓ̄ > eps_T ? 2 * sin(abs(θ) / 2) / ℓ̄ : zero(T)
    return nothing
end

function _curvature_vecs_kernel!(out, pts, n::Int32, N::Int32, closed::Bool)
    i = _gidx()
    nw = closed ? n : n - Int32(2)
    i > nw && return nothing
    T  = eltype(pts)
    iv = closed ? i : i + Int32(1)
    ip = _adj_gpu(iv - Int32(1), n, closed)
    iq = _adj_gpu(iv + Int32(1), n, closed)

    n1sq = zero(T); n2sq = zero(T); dotv = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, iv] - pts[d, ip]
        e2 = pts[d, iq] - pts[d, iv]
        n1sq = fma(e1, e1, n1sq)
        n2sq = fma(e2, e2, n2sq)
        dotv = fma(e1, e2, dotv)
    end
    n1 = sqrt(n1sq); n2 = sqrt(n2sq)
    eps_T = eps(T)
    if n1 ≤ eps_T || n2 ≤ eps_T
        @inbounds for d in Int32(1):N; out[d, i] = zero(T); end
        return nothing
    end

    cross_sq = max(zero(T), n1sq * n2sq - dotv * dotv)
    θ        = atan(sqrt(cross_sq), dotv)
    ℓ̄        = (n1 + n2) / 2
    κ        = ℓ̄ > eps_T ? abs(θ) / ℓ̄ : zero(T)

    dir_sq = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, iv] - pts[d, ip]
        e2 = pts[d, iq] - pts[d, iv]
        v  = e2/n2 - e1/n1
        dir_sq = fma(v, v, dir_sq)
    end
    dir_norm = sqrt(dir_sq)
    if dir_norm ≤ eps_T
        @inbounds for d in Int32(1):N; out[d, i] = zero(T); end
        return nothing
    end
    scale = κ / dir_norm
    @inbounds for d in Int32(1):N
        e1 = pts[d, iv] - pts[d, ip]
        e2 = pts[d, iq] - pts[d, iv]
        out[d, i] = (e2/n2 - e1/n1) * scale
    end
    return nothing
end

function DiscreteCurves.arc_length(g::GPUCurve{N,T}) where {N,T}
    ne = Int(isclosed(g) ? nvertices(g) : nvertices(g) - 1)
    ne == 0 && return zero(T)
    buf = CUDA.zeros(T, ne)
    _launch1d(_edge_lengths_kernel!, ne, buf, g.points.data,
              Int32(nvertices(g)), Int32(N), isclosed(g))
    sum(buf)                               # CUDA.jl parallel reduction
end

function DiscreteCurves.edge_lengths(g::GPUCurve{N,T}) where {N,T}
    ne = Int(isclosed(g) ? nvertices(g) : nvertices(g) - 1)
    ne == 0 && return T[]
    buf = CUDA.zeros(T, ne)
    _launch1d(_edge_lengths_kernel!, ne, buf, g.points.data,
              Int32(nvertices(g)), Int32(N), isclosed(g))
    Array(buf)                             
end

function DiscreteCurves.cumulative_length(g::GPUCurve{N,T}) where {N,T}
    ls = DiscreteCurves.edge_lengths(g)
    n  = nvertices(g)
    s  = Vector{T}(undef, n)
    s[1] = zero(T)
    @inbounds for i in 1:min(length(ls), n - 1)
        s[i + 1] = s[i] + ls[i]
    end
    s
end

function DiscreteCurves.cfl_timestep(g::GPUCurve{N,T};
                                      safety::Real = T(0.25)) where {N,T}
    ne = Int(isclosed(g) ? nvertices(g) : nvertices(g) - 1)
    ne == 0 && return T(safety)
    buf = CUDA.zeros(T, ne)
    _launch1d(_edge_sq_kernel!, ne, buf, g.points.data,
              Int32(nvertices(g)), Int32(N), isclosed(g))
    T(safety) * minimum(buf)              
end

function DiscreteCurves.curvatures(g::GPUCurve{N,T}) where {N,T}
    DiscreteCurves.curvatures(g, TurningAngle())
end

function DiscreteCurves.curvatures(g::GPUCurve{N,T}, ::TurningAngle) where {N,T}
    nw = isclosed(g) ? nvertices(g) : nvertices(g) - 2
    nw ≤ 0 && return T[]
    buf = CUDA.zeros(T, nw)
    _launch1d(_curvatures_ta_kernel!, nw, buf, g.points.data,
              Int32(nvertices(g)), Int32(N), isclosed(g))
    Array(buf)
end

function DiscreteCurves.curvatures(g::GPUCurve{N,T}, ::SteinerCurvature) where {N,T}
    nw = isclosed(g) ? nvertices(g) : nvertices(g) - 2
    nw ≤ 0 && return T[]
    buf = CUDA.zeros(T, nw)
    _launch1d(_curvatures_steiner_kernel!, nw, buf, g.points.data,
              Int32(nvertices(g)), Int32(N), isclosed(g))
    Array(buf)
end

function DiscreteCurves.curvature_vectors(g::GPUCurve{N,T}) where {N,T}
    DiscreteCurves.curvature_vectors(g, TurningAngle())
end

function DiscreteCurves.curvature_vectors(g::GPUCurve{N,T}, ::TurningAngle) where {N,T}
    nw = isclosed(g) ? nvertices(g) : nvertices(g) - 2
    nw ≤ 0 && return SVector{N,T}[]
    buf = CUDA.zeros(T, N, nw)
    _launch1d(_curvature_vecs_kernel!, nw, buf, g.points.data,
              Int32(nvertices(g)), Int32(N), isclosed(g))
    mat = Array(buf)
    [SVector{N,T}(ntuple(d -> mat[d, j], Val(N))) for j in 1:nw]
end

function DiscreteCurves.signed_area(g::GPUCurve{2,T}) where {T}
    ne = Int(isclosed(g) ? nvertices(g) : nvertices(g) - 1)
    ne == 0 && return zero(T)
    buf = CUDA.zeros(T, ne)
    _launch1d(_area_terms_kernel!, ne, buf, g.points.data,
              Int32(nvertices(g)), isclosed(g))
    sum(buf) / 2
end

function _steiner_vel_kernel!(vel, pts, n::Int32, N::Int32, closed::Bool, tangential::Bool)
    i = _gidx()
    i > n && return nothing
    T = eltype(pts)

    if !closed && (i == Int32(1) || i == n)
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    ip  = _adj_gpu(i - Int32(1), n, closed)
    iq2 = _adj_gpu(i + Int32(1), n, closed)

    n1sq = zero(T); n2sq = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, i] - pts[d, ip]
        e2 = pts[d, iq2] - pts[d, i]
        n1sq = fma(e1, e1, n1sq)
        n2sq = fma(e2, e2, n2sq)
    end
    n1 = sqrt(n1sq); n2 = sqrt(n2sq)
    eps_T = eps(T)

    if n1 ≤ eps_T || n2 ≤ eps_T
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    @inbounds for d in Int32(1):N
        e1 = pts[d, i] - pts[d, ip]
        e2 = pts[d, iq2] - pts[d, i]
        t2d = e2 / n2
        v   = t2d - e1 / n1
        vel[d, i] = tangential ? v + (n2 - n1) / 2 * t2d : v
    end
    return nothing
end


function _turning_vel_kernel!(vel, pts, n::Int32, N::Int32, closed::Bool, tangential::Bool)
    i = _gidx()
    i > n && return nothing
    T = eltype(pts)

    if !closed && (i == Int32(1) || i == n)
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    ip  = _adj_gpu(i - Int32(1), n, closed)
    iq2 = _adj_gpu(i + Int32(1), n, closed)

    n1sq = zero(T); n2sq = zero(T); dotv = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, i] - pts[d, ip]
        e2 = pts[d, iq2] - pts[d, i]
        n1sq = fma(e1, e1, n1sq)
        n2sq = fma(e2, e2, n2sq)
        dotv = fma(e1, e2, dotv)
    end
    n1 = sqrt(n1sq); n2 = sqrt(n2sq)
    eps_T = eps(T)

    if n1 ≤ eps_T || n2 ≤ eps_T
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    cross_sq = max(zero(T), n1sq * n2sq - dotv * dotv)
    θ        = atan(sqrt(cross_sq), dotv)
    ℓ̄        = (n1 + n2) / 2
    κ        = ℓ̄ > eps_T ? abs(θ) / ℓ̄ : zero(T)

    dir_sq = zero(T)
    @inbounds for d in Int32(1):N
        e1 = pts[d, i] - pts[d, ip]
        e2 = pts[d, iq2] - pts[d, i]
        v  = e2/n2 - e1/n1
        dir_sq = fma(v, v, dir_sq)
    end
    dir_norm = sqrt(dir_sq)
    if dir_norm ≤ eps_T
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    scale = κ / dir_norm
    @inbounds for d in Int32(1):N
        e1  = pts[d, i] - pts[d, ip]
        e2  = pts[d, iq2] - pts[d, i]
        t2d = e2 / n2
        v   = (t2d - e1/n1) * scale
        vel[d, i] = tangential ? v + (n2 - n1) / 2 * t2d : v
    end
    return nothing
end


function _euler_step_kernel!(pts, vel, n::Int32, N::Int32, τ::T, max_v_sq) where {T}
    i = _gidx()
    i > n && return nothing
    v_sq = zero(T)
    @inbounds for d in Int32(1):N
        vdi = vel[d, i]
        v_sq = fma(vdi, vdi, v_sq)
        pts[d, i] = fma(τ, vdi, pts[d, i])
    end
    CUDA.atomic_max!(pointer(max_v_sq, 1), v_sq)
    return nothing
end


function _gpu_step!(vel::CuMatrix{T}, g::GPUCurve{N,T},
                    law::DiscreteCurves.Flows.CurvatureFlow{SteinerCurvature},
                    τ::T) where {N,T}
    n = Int32(nvertices(g))
    _launch1d(_steiner_vel_kernel!, Int(n), vel, g.points.data,
              n, Int32(N), isclosed(g), law.tangential)
    max_v_sq = CUDA.zeros(T, 1)
    _launch1d(_euler_step_kernel!, Int(n), g.points.data, vel,
              n, Int32(N), τ, max_v_sq)
    sqrt(Array(max_v_sq)[1])
end

function _gpu_step!(vel::CuMatrix{T}, g::GPUCurve{N,T},
                    law::DiscreteCurves.Flows.CurvatureFlow{TurningAngle},
                    τ::T) where {N,T}
    n = Int32(nvertices(g))
    _launch1d(_turning_vel_kernel!, Int(n), vel, g.points.data,
              n, Int32(N), isclosed(g), law.tangential)
    max_v_sq = CUDA.zeros(T, 1)
    _launch1d(_euler_step_kernel!, Int(n), g.points.data, vel,
              n, Int32(N), τ, max_v_sq)
    sqrt(Array(max_v_sq)[1])
end

function _gpu_step!(vel, g::GPUCurve, law, τ)
    throw(GPUNotSupportedError(
        "No GPU kernel for $(typeof(law)).\n" *
        "Supported: CurvatureFlow{SteinerCurvature}, CurvatureFlow{TurningAngle}.\n" *
        "To use another law, evolve on CPU: evolve(to_cpu(gc), law; ...).\n" *
        "To add GPU support: specialise _gpu_step! in DiscreteCurvesCUDAExt."))
end


"""
    evolve(gc, law; kwargs...) → FlowResult

GPU-accelerated curve evolution with identical syntax to the CPU `evolve`.

**Supported laws** (native GPU kernels):
- `CurvatureFlow(SteinerCurvature(); tangential=true/false)`
- `CurvatureFlow(TurningAngle();    tangential=true/false)`

Other laws throw `GPUNotSupportedError`.

**Notes**
- `resample_every` is not yet supported for `GPUCurve` (throws if > 0).
- Tracking callbacks receive the *live* `GPUCurve`; any GPU-overridden function
  (e.g. `arc_length`, `signed_area`) runs as a GPU kernel inside the loop.
- Snapshots (`save_every`) download to CPU as `DiscreteCurve` (one D→H per save).
- The final `.curve` in `FlowResult` is a CPU `DiscreteCurve`.

```julia
gc = to_gpu(fourier_curve(4096; seed=1))

r = evolve(gc, CurvatureFlow(SteinerCurvature(); tangential=true);
           nsteps     = 5000,
           save_every = 100,
           track      = [:L => arc_length,     # GPU kernel each tracked step
                         :A => signed_area,
                         :κ => c -> maximum(curvatures(c))])

r.curve       # final DiscreteCurve
r[:L]         # Vector of arc lengths
plot(r[100])  # curve at step 100
```
"""
function DiscreteCurves.evolve(g0::GPUCurve{N,T},
                                law;
                                τ              :: Union{Nothing,Real} = nothing,
                                cfl_safety     :: Real                = T(0.25),
                                nsteps         :: Int                 = 200,
                                tol            :: Real                = T(eps(T)^(1/3)),
                                resample_every :: Int                 = 0,
                                save_every     :: Int                 = 0,
                                track          :: AbstractVector      = Pair{Symbol}[],
                                track_every    :: Int                 = 1,
                                callback                              = nothing) where {N,T}

    resample_every > 0 && throw(GPUNotSupportedError(
        "resample_every is not yet supported for GPUCurve. " *
        "Options: (a) evolve on CPU with resample_every, " *
        "(b) evolve on GPU without it, " *
        "(c) use DiscreteCurves.resample(to_cpu(gc), n) |> to_gpu between evolve calls."))

    g   = GPUCurve{N,T}(copy(g0.points), g0.closed)
    vel = CUDA.zeros(T, N, nvertices(g))

    τ_cur  = τ === nothing ?
             T(DiscreteCurves.cfl_timestep(g; safety = T(cfl_safety))) : T(τ)
    τ_cur > zero(T) || throw(ArgumentError("τ must be positive, got $τ_cur"))

    trackers = [(Symbol(k), f) for (k, f) in track]
    history  = Dict{Symbol, AbstractVector}()
    for (sym, f) in trackers
        s = f(g)
        history[sym] = typeof(s)[s]
    end

    max_v_hist   = T[]
    snapshots    = DiscreteCurve{N,T}[]
    snap_times   = T[]
    total_time   = zero(T)
    actual_steps = 0
    tol_T        = T(tol)

    if save_every > 0                    
        push!(snapshots,  to_cpu(g))
        push!(snap_times, zero(T))
    end

    for σ in 1:nsteps
        actual_steps  = σ
        max_v         = _gpu_step!(vel, g, law, τ_cur) 
        total_time   += τ_cur
        push!(max_v_hist, max_v)

        if save_every > 0 && σ % save_every == 0
            push!(snapshots,  to_cpu(g))
            push!(snap_times, total_time)
        end

        if !isempty(trackers) && σ % track_every == 0
            for (sym, f) in trackers
                push!(history[sym], f(g)) 
            end
        end

        callback !== nothing && callback(σ, g, max_v) && break
        max_v < tol_T && break
    end

    history[:max_v] = max_v_hist
    DiscreteCurves.Flows.FlowResult{N,T}(
        to_cpu(g), actual_steps, total_time,
        max_v_hist, history, snapshots, snap_times)
end

function DiscreteCurves.evolve(c::DiscreteCurve{N,T},
                                law::DiscreteCurves.Flows.CurvatureFlow; kwargs...) where {N,T}
    DiscreteCurves.evolve(to_gpu(c), law; kwargs...)
end

end