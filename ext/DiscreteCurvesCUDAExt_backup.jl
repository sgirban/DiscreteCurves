"""
    DiscreteCurvesCUDAExt

GPU-accelerated curve flows and geometry using CUDA.jl.

## Design

All GPU operations use a *column-major matrix* layout:

    pts_gpu :: CuMatrix{T}   # size (N, n_vertices)

This is cache-friendly for CUDA (threads read consecutive memory along the first
axis = the N spatial coordinates of one vertex) and avoids the
`Vector{SVector{…}}` indirection that is not GPU-compatible.

## Public API (added by this extension)

```julia
using CUDA
using DiscreteCurves

c   = fourier_curve(1024)
gr  = to_gpu(c)           # GPUCurve{2,Float64}
c2  = to_cpu(gr)          # back to DiscreteCurve
gr2 = gpu_evolve(gr, SteinerCurvature(); nsteps=500)
```

Zygote differentiability is provided through the companion
`DiscreteCurvesChainRulesExt` — load both by `using CUDA, Zygote`.
"""
module DiscreteCurvesCUDAExt

using DiscreteCurves
using DiscreteCurves.AbstractTypes: AbstractDiscreteCurve
using DiscreteCurves.CurveTypes: DiscreteCurve
using DiscreteCurves.CurveTopology: isclosed, nvertices, nedges
using DiscreteCurves.Curvature: SteinerCurvature, TurningAngle, TangentCurvature,
                                 OsculatingCircle, LengthWeightedCurvature,
                                 SignedCurvature2D, CustomCurvature
using CUDA
using CUDA: CuArray, @cuda, threadIdx, blockIdx, blockDim, CUDAKernelError
using LinearAlgebra
using StaticArrays

"""
    GPUCurve{N, T}

A discrete curve whose vertices live on the GPU.

Fields:
- `pts :: CuMatrix{T}` — layout is `(N, n_vertices)` (column = one vertex).
- `closed :: Bool`

Construct via `to_gpu(c::DiscreteCurve)`. Convert back with `to_cpu(g)`.

!!! note
    `GPUCurve` does **not** subtype `AbstractDiscreteCurve` because that type
    requires `SVector` indexing which is meaningless on device memory.  Instead
    it exposes a minimal compatible interface (`nvertices`, `nedges`, `isclosed`).
"""
struct GPUCurve{N, T <: AbstractFloat}
    pts    :: CuMatrix{T}
    closed :: Bool

    function GPUCurve{N,T}(pts::CuMatrix{T}, closed::Bool) where {N,T}
        size(pts, 1) == N || throw(DimensionMismatch(
            "pts has $(size(pts,1)) rows but N=$N"))
        new{N,T}(pts, closed)
    end
end

DiscreteCurves.CurveTopology.nvertices(g::GPUCurve) = size(g.pts, 2)
DiscreteCurves.CurveTopology.nedges(g::GPUCurve)    =
    isclosed(g) ? nvertices(g) : nvertices(g) - 1
DiscreteCurves.CurveTopology.isclosed(g::GPUCurve)  = g.closed

Base.eltype(::GPUCurve{N,T}) where {N,T} = T
Base.ndims(::GPUCurve{N})    where  N    = N

function Base.show(io::IO, g::GPUCurve{N,T}) where {N,T}
    topo = g.closed ? "closed" : "open"
    print(io, "GPUCurve{$N,$T}($topo, $(nvertices(g)) vertices on GPU)")
end

export GPUCurve

"""
    to_gpu(c::DiscreteCurve{N,T}) → GPUCurve{N,T}

Upload a CPU curve to the GPU.  Allocates one contiguous `CuMatrix{T}` of
size `(N, nvertices(c))`.

```julia
c  = fourier_curve(1024)
gc = to_gpu(c)
```
"""
function to_gpu(c::DiscreteCurve{N,T}) where {N,T}
    n   = nvertices(c)
    cpu = Matrix{T}(undef, N, n)
    @inbounds for j in 1:n, d in 1:N
        cpu[d, j] = c.points[j][d]
    end
    GPUCurve{N,T}(CuArray(cpu), c.closed)
end

"""
    to_cpu(g::GPUCurve{N,T}) → DiscreteCurve{N,T}

Download a GPU curve to the CPU.

```julia
c2 = to_cpu(gc)
```
"""
function to_cpu(g::GPUCurve{N,T}) where {N,T}
    cpu = Array(g.pts)                    # single synchronised D→H transfer
    n   = size(cpu, 2)
    pts = Vector{SVector{N,T}}(undef, n)
    @inbounds for j in 1:n
        pts[j] = SVector{N,T}(ntuple(d -> cpu[d,j], Val(N)))
    end
    DiscreteCurve{N,T}(pts, g.closed)
end

export to_gpu, to_cpu

# ─────────────────────────────────────────────────────────────────────────────
#  GPU kernels
# ─────────────────────────────────────────────────────────────────────────────
#
# All kernels share the same indexing convention:
#   i  ∈ 1:n   — vertex index (1-based, standard Julia)
#   d  ∈ 1:N   — spatial dimension
#   pts[d, i]  — d-th coordinate of vertex i
#
# Topology helpers operate on the thread-local i — no shared memory needed
# because each vertex only needs its two neighbours.

# ── Modular index (1-based, closed topology) ────────────────────────────────
@inline _gpu_mod1(i::Int32, n::Int32) = ((i - Int32(1)) % n + n) % n + Int32(1)

@inline function _gpu_adj(i::Int32, n::Int32, closed::Bool)
    closed ? _gpu_mod1(i, n) : clamp(i, Int32(1), n)
end

# ── 2-norm of column j of pts ───────────────────────────────────────────────
@inline function _gpu_col_norm(pts, N::Int32, j::Int32)
    s = zero(eltype(pts))
    @inbounds for d in Int32(1):N
        v = pts[d, j]
        s = fma(v, v, s)
    end
    sqrt(s)
end

# ── Unit difference: d2/|d2| - d1/|d1| (the κ·N̂ direction vector) ─────────
@inline function _gpu_tangent_diff!(out, pts, N::Int32,
                                    ip::Int32, iq::Int32, iq2::Int32)
    # d1 = pts[:,iq] - pts[:,ip],  d2 = pts[:,iq2] - pts[:,iq]
    n1 = zero(eltype(pts)); n2 = zero(eltype(pts))
    @inbounds for d in Int32(1):N
        d1 = pts[d, iq] - pts[d, ip]
        d2 = pts[d, iq2] - pts[d, iq]
        n1 = fma(d1, d1, n1)
        n2 = fma(d2, d2, n2)
        out[d] = d1   # store raw d1 temporarily
    end
    n1 = sqrt(n1); n2 = sqrt(n2)
    eps_T = eps(eltype(pts))^2
    if n1 ≤ eps_T || n2 ≤ eps_T
        @inbounds for d in Int32(1):N; out[d] = zero(eltype(pts)); end
        return false
    end
    @inbounds for d in Int32(1):N
        d1 = pts[d, iq]  - pts[d, ip]
        d2 = pts[d, iq2] - pts[d, iq]
        out[d] = d2/n2 - d1/n1
    end
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
#  Kernel: SteinerCurvature velocities
#  v_i = t̂_out - t̂_in  (the κ^B * N̂ vector, no division by ℓ̄)
# ─────────────────────────────────────────────────────────────────────────────
function _steiner_vel_kernel!(vel, pts, n::Int32, N::Int32, closed::Bool,
                               tangential::Bool)
    i = Int32(threadIdx().x) + Int32(blockDim().x) * (Int32(blockIdx().x) - Int32(1))
    i > n && return nothing

    if !closed && (i == Int32(1) || i == n)
        @inbounds for d in Int32(1):N; vel[d, i] = zero(eltype(pts)); end
        return nothing
    end

    ip  = _gpu_adj(i - Int32(1), n, closed)
    iq2 = _gpu_adj(i + Int32(1), n, closed)

    T   = eltype(pts)
    buf = MArray{Tuple{8}, T}(undef)   # max N=8 (stack-allocated)

    ok = _gpu_tangent_diff!(buf, pts, N, ip, i, iq2)
    if !ok
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    if tangential
        # DeTurck tangential term: (|d2| - |d1|)/2 * t̂_out
        n1 = zero(T); n2 = zero(T)
        @inbounds for d in Int32(1):N
            d1v = pts[d, i]   - pts[d, ip]
            d2v = pts[d, iq2] - pts[d, i]
            n1  = fma(d1v, d1v, n1)
            n2  = fma(d2v, d2v, n2)
        end
        n1 = sqrt(n1); n2 = sqrt(n2)
        targ = (n2 - n1) / 2
        @inbounds for d in Int32(1):N
            d2v = pts[d, iq2] - pts[d, i]
            vel[d, i] = buf[d] + (n2 > eps(T) ? targ * d2v / n2 : zero(T))
        end
    else
        @inbounds for d in Int32(1):N
            vel[d, i] = buf[d]
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
#  Kernel: TurningAngle velocities
#  κ^A * N̂:   κ = |θ| / ℓ̄,   direction = normalised tangent-diff
# ─────────────────────────────────────────────────────────────────────────────
function _turning_vel_kernel!(vel, pts, n::Int32, N::Int32, closed::Bool,
                               tangential::Bool)
    i = Int32(threadIdx().x) + Int32(blockDim().x) * (Int32(blockIdx().x) - Int32(1))
    i > n && return nothing

    T = eltype(pts)

    if !closed && (i == Int32(1) || i == n)
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    ip  = _gpu_adj(i - Int32(1), n, closed)
    iq2 = _gpu_adj(i + Int32(1), n, closed)

    buf = MArray{Tuple{8}, T}(undef)
    ok  = _gpu_tangent_diff!(buf, pts, N, ip, i, iq2)
    if !ok
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    # |buf| = κ^B * magnitude; we need to rescale to κ^A = |θ|/ℓ̄
    buf_norm = zero(T)
    @inbounds for d in Int32(1):N; v = buf[d]; buf_norm = fma(v,v,buf_norm); end
    buf_norm = sqrt(buf_norm)
    buf_norm ≤ eps(T) && begin
        @inbounds for d in Int32(1):N; vel[d, i] = zero(T); end
        return nothing
    end

    # Compute ℓ̄ = (|e_{i-1}| + |e_i|) / 2
    ell1 = zero(T); ell2 = zero(T)
    @inbounds for d in Int32(1):N
        v1 = pts[d,i] - pts[d,ip];   ell1 = fma(v1,v1,ell1)
        v2 = pts[d,iq2] - pts[d,i];  ell2 = fma(v2,v2,ell2)
    end
    ell_bar = (sqrt(ell1) + sqrt(ell2)) / 2

    # Compute turning angle: sin(θ) ≈ |d1×d2|/(|d1||d2|), cos(θ) = d1·d2/(|d1||d2|)
    # For N=2: atan(cross_z, dot)
    # For N≥3: atan(|cross|, dot)
    dot_val = zero(T); cross_sq = zero(T)
    @inbounds for d in Int32(1):N
        t1d = pts[d,i]   - pts[d,ip]
        t2d = pts[d,iq2] - pts[d,i]
        dot_val = fma(t1d, t2d, dot_val)
    end
    # |t1×t2|² = |t1|²|t2|² - (t1·t2)²
    n1_sq = ell1; n2_sq = ell2
    cross_sq = max(zero(T), n1_sq * n2_sq - dot_val * dot_val)
    θ = atan(sqrt(cross_sq), dot_val)

    κ = ell_bar > eps(T) ? abs(θ) / ell_bar : zero(T)
    scale = κ / buf_norm

    if tangential
        n1 = sqrt(ell1); n2 = sqrt(ell2)
        targ = (n2 - n1) / 2
        @inbounds for d in Int32(1):N
            d2v = pts[d, iq2] - pts[d, i]
            tang = n2 > eps(T) ? targ * d2v / n2 : zero(T)
            vel[d, i] = buf[d] * scale + tang
        end
    else
        @inbounds for d in Int32(1):N
            vel[d, i] = buf[d] * scale
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
#  Kernel: Euler step  pts ← pts + τ * vel, return max |vel|² (atomicMax)
# ─────────────────────────────────────────────────────────────────────────────
function _euler_step_kernel!(pts, vel, n::Int32, N::Int32, τ::T,
                              max_v_sq::CuDeviceArray{T}) where T
    i = Int32(threadIdx().x) + Int32(blockDim().x) * (Int32(blockIdx().x) - Int32(1))
    i > n && return nothing

    v_sq = zero(T)
    @inbounds for d in Int32(1):N
        vdi   = vel[d, i]
        v_sq  = fma(vdi, vdi, v_sq)
        pts[d, i] = fma(τ, vdi, pts[d, i])
    end
    CUDA.atomic_max!(pointer(max_v_sq, 1), v_sq)
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
#  GPU-side CFL timestep  (min edge-length² × safety)
# ─────────────────────────────────────────────────────────────────────────────
function _min_edge_sq_kernel!(out::CuDeviceArray{T}, pts,
                               n::Int32, N::Int32, closed::Bool) where T
    i = Int32(threadIdx().x) + Int32(blockDim().x) * (Int32(blockIdx().x) - Int32(1))
    ne = closed ? n : n - Int32(1)
    i > ne && return nothing

    j  = _gpu_adj(i + Int32(1), n, closed)
    sq = zero(T)
    @inbounds for d in Int32(1):N
        dv = pts[d,j] - pts[d,i]
        sq = fma(dv, dv, sq)
    end
    CUDA.atomic_min!(pointer(out, 1), sq)
    return nothing
end

function gpu_cfl_timestep(g::GPUCurve{N,T}; safety::Real = T(0.25)) where {N,T}
    n      = Int32(nvertices(g))
    Ni     = Int32(N)
    closed = isclosed(g)
    ne     = closed ? n : n - Int32(1)
    buf    = CuArray(fill(typemax(T), 1))
    threads = min(Int(ne), 256)
    blocks  = cld(Int(ne), threads)
    @cuda threads=threads blocks=blocks _min_edge_sq_kernel!(buf, g.pts, n, Ni, closed)
    T(safety) * Array(buf)[1]
end

# ─────────────────────────────────────────────────────────────────────────────
#  Dispatch helper — launch the right kernel for a given curvature model
# ─────────────────────────────────────────────────────────────────────────────

function _launch_vel_kernel!(vel, g::GPUCurve{N,T},
                              ::SteinerCurvature, tangential::Bool) where {N,T}
    n  = Int32(nvertices(g)); Ni = Int32(N)
    threads = min(Int(n), 256)
    blocks  = cld(Int(n), threads)
    @cuda threads=threads blocks=blocks _steiner_vel_kernel!(
        vel, g.pts, n, Ni, isclosed(g), tangential)
end

function _launch_vel_kernel!(vel, g::GPUCurve{N,T},
                              ::TurningAngle, tangential::Bool) where {N,T}
    n  = Int32(nvertices(g)); Ni = Int32(N)
    threads = min(Int(n), 256)
    blocks  = cld(Int(n), threads)
    @cuda threads=threads blocks=blocks _turning_vel_kernel!(
        vel, g.pts, n, Ni, isclosed(g), tangential)
end

# Fallback: download to CPU, run CPU law, upload result.
# This handles CustomCurvature, PointwiseFlow, GradientFlow, etc.
function _launch_vel_kernel_cpu_fallback!(vel_gpu, g::GPUCurve{N,T}, law,
                                           tangential::Bool) where {N,T}
    c_cpu   = to_cpu(g)
    n       = nvertices(c_cpu)
    vel_cpu = Vector{SVector{N,T}}(undef, n)

    if law isa DiscreteCurves.Flows.CurvatureFlow
        DiscreteCurves.Flows.compute_velocities!(vel_cpu, c_cpu, law)
    else
        DiscreteCurves.Flows.compute_velocities!(vel_cpu, c_cpu, law)
    end

    # Pack vel_cpu (Vector{SVector}) → Matrix → CuArray
    mat = Matrix{T}(undef, N, n)
    @inbounds for j in 1:n, d in 1:N; mat[d,j] = vel_cpu[j][d]; end
    copyto!(vel_gpu, CuArray(mat))
end

# ─────────────────────────────────────────────────────────────────────────────
#  gpu_evolve_step!
# ─────────────────────────────────────────────────────────────────────────────

"""
    gpu_evolve_step!(vel, g, law, τ) → max_speed::T

One explicit-Euler GPU step. Mutates `g.pts` in place.
`vel` is a pre-allocated `CuMatrix{T}` of the same size as `g.pts`.

Returns the maximum vertex speed (√max|vᵢ|²) of this step.
"""
function gpu_evolve_step!(vel::CuMatrix{T}, g::GPUCurve{N,T},
                           law, τ::T) where {N,T}
    # 1. Compute velocities
    if law isa DiscreteCurves.Flows.CurvatureFlow{SteinerCurvature}
        _launch_vel_kernel!(vel, g, SteinerCurvature(), law.tangential)
    elseif law isa DiscreteCurves.Flows.CurvatureFlow{TurningAngle}
        _launch_vel_kernel!(vel, g, TurningAngle(), law.tangential)
    else
        _launch_vel_kernel_cpu_fallback!(vel, g, law, false)
    end

    # 2. Euler step + max |v|
    n       = Int32(nvertices(g))
    Ni      = Int32(N)
    max_v_sq = CuArray(zeros(T, 1))
    threads  = min(Int(n), 256)
    blocks   = cld(Int(n), threads)
    @cuda threads=threads blocks=blocks _euler_step_kernel!(
        g.pts, vel, n, Ni, τ, max_v_sq)

    sqrt(Array(max_v_sq)[1])
end

export gpu_evolve_step!

# ─────────────────────────────────────────────────────────────────────────────
#  gpu_evolve — high-level driver mirroring CPU evolve API
# ─────────────────────────────────────────────────────────────────────────────

"""
    gpu_evolve(g, law; kwargs…) → GPUFlowResult

GPU-accelerated analogue of `evolve`.  Accepts either a `GPUCurve` or a
`DiscreteCurve` (auto-uploaded).

# Supported laws (native GPU kernels)
- `CurvatureFlow(SteinerCurvature(); tangential=false/true)`
- `CurvatureFlow(TurningAngle();    tangential=false/true)`
- Any other law: automatically falls back to CPU per-step (slower but correct).

# Keyword arguments

| kwarg           | default    | meaning                                          |
|-----------------|------------|--------------------------------------------------|
| `τ`             | auto CFL   | fixed timestep (`nothing` = auto)                |
| `cfl_safety`    | 0.25       | CFL safety factor                                |
| `nsteps`        | 200        | maximum number of steps                          |
| `tol`           | eps^(1/3)  | stop when max speed < tol                        |
| `save_every`    | 0          | download+save CPU snapshot every k steps         |
| `track`         | []         | `[:L => (g -> arc_length_gpu(g)), …]`            |
| `callback`      | nothing    | `(σ, g, max_v) → Bool`; return true to stop     |

# Returns `GPUFlowResult`
- `.curve`  — final `GPUCurve` (still on GPU)
- `.cpu_curve` — `to_cpu` of final curve
- `.nsteps`, `.time`, `.max_v` — as in CPU `FlowResult`
- `.snapshots` — `Vector{DiscreteCurve}` of saved CPU snapshots

# Example

```julia
using CUDA, DiscreteCurves

c  = fourier_curve(2048; seed=1)
gc = to_gpu(c)

r  = gpu_evolve(gc, CurvatureFlow(SteinerCurvature(); tangential=true);
                nsteps = 2000, tol = 1e-7,
                save_every = 100)

r.cpu_curve   # final result back on CPU
r.snapshots   # Vector of 20 DiscreteCurve snapshots
```
"""
struct GPUFlowResult{N, T}
    curve     :: GPUCurve{N,T}
    cpu_curve :: DiscreteCurve{N,T}
    nsteps    :: Int
    time      :: T
    max_v     :: Vector{T}
    history   :: Dict{Symbol, Vector}
    snapshots :: Vector{DiscreteCurve{N,T}}
    snap_times:: Vector{T}
end

function Base.show(io::IO, r::GPUFlowResult{N,T}) where {N,T}
    print(io, "GPUFlowResult{$N,$T}  σ=$(r.nsteps)  t=$(round(r.time;sigdigits=4))",
          "  snapshots=$(length(r.snapshots))")
end

export GPUFlowResult

function gpu_evolve(g0::GPUCurve{N,T},
                    law;
                    τ           :: Union{Nothing,Real} = nothing,
                    cfl_safety  :: Real = T(0.25),
                    nsteps      :: Int  = 200,
                    tol         :: Real = T(eps(T)^(1/3)),
                    save_every  :: Int  = 0,
                    track       :: AbstractVector = Pair{Symbol}[],
                    callback            = nothing) where {N,T}

    # Clone pts so we mutate a working copy
    g      = GPUCurve{N,T}(copy(g0.pts), g0.closed)
    vel    = CUDA.zeros(T, N, nvertices(g))
    τ_cur  = τ === nothing ? T(gpu_cfl_timestep(g; safety=T(cfl_safety))) : T(τ)
    τ_cur > zero(T) || throw(ArgumentError("τ must be positive, got $τ_cur"))

    trackers = [(Symbol(k), f) for (k,f) in track]
    history  = Dict{Symbol, Vector}()
    for (sym, f) in trackers
        s = f(g)
        history[sym] = [s]
    end

    max_v_hist  = T[]
    snapshots   = DiscreteCurve{N,T}[]
    snap_times  = T[]
    total_time  = zero(T)
    actual_steps = 0

    if save_every > 0
        push!(snapshots,  to_cpu(g))
        push!(snap_times, zero(T))
    end

    for σ in 1:nsteps
        actual_steps = σ
        max_v   = gpu_evolve_step!(vel, g, law, τ_cur)
        total_time += τ_cur
        push!(max_v_hist, max_v)

        if save_every > 0 && σ % save_every == 0
            push!(snapshots,  to_cpu(g))
            push!(snap_times, total_time)
        end

        for (sym, f) in trackers
            push!(history[sym], f(g))
        end

        callback !== nothing && callback(σ, g, max_v) && break
        max_v < T(tol) && break
    end

    cpu_final = to_cpu(g)
    GPUFlowResult{N,T}(g, cpu_final, actual_steps, total_time,
                       max_v_hist, history, snapshots, snap_times)
end

# Accept CPU DiscreteCurve directly
function gpu_evolve(c::DiscreteCurve{N,T}, law; kwargs...) where {N,T}
    gpu_evolve(to_gpu(c), law; kwargs...)
end

export gpu_evolve

# ─────────────────────────────────────────────────────────────────────────────
#  GPU geometry helpers (for use in track= callbacks)
# ─────────────────────────────────────────────────────────────────────────────

"""
    arc_length_gpu(g::GPUCurve{N,T}) → T

Total arc length computed on CPU after a one-time sync.  Use sparingly in
`track=` callbacks — it serialises the GPU pipeline.
"""
arc_length_gpu(g::GPUCurve) = DiscreteCurves.arc_length(to_cpu(g))

"""
    max_curvature_gpu(g) → T

Maximum discrete curvature (TurningAngle model) after a CPU sync.
"""
max_curvature_gpu(g::GPUCurve) = maximum(DiscreteCurves.curvatures(to_cpu(g)))

export arc_length_gpu, max_curvature_gpu

# ─────────────────────────────────────────────────────────────────────────────
#  Differentiable GPU arc-length (via CUDA.jl array operations)
#  Works with Zygote if loaded; Zygote can differentiate through CuArray ops.
# ─────────────────────────────────────────────────────────────────────────────

"""
    gpu_arc_length(pts::CuMatrix{T}) → T

Differentiable arc-length from a raw `(N, n)` CuMatrix.
Uses `sum(norm.(diff(cols)))` so Zygote can differentiate through it.

```julia
using Zygote, CUDA, DiscreteCurves
pts = to_gpu(fourier_curve(256)).pts
grad = Zygote.gradient(gpu_arc_length, pts)[1]   # (N, n) CuMatrix
```
"""
function gpu_arc_length(pts::CuMatrix{T}) where T
    # pts: (N, n)  — each column is a vertex
    # Compute squared differences column-wise, sum rows, sqrt, sum over edges
    n   = size(pts, 2)
    dp  = @view(pts[:, 2:n]) .- @view(pts[:, 1:(n-1)])   # (N, n-1)
    sum(sqrt.(sum(dp .^ 2; dims=1)))
end

export gpu_arc_length

end # module DiscreteCurvesCUDAExt