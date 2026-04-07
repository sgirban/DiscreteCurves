# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  test/runtests.jl                                                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

using Test
using DiscreteCurves
using StaticArrays
using LinearAlgebra

# ── helpers ───────────────────────────────────────────────────────────────────
ngon(n, R=1.0) = ClosedCurve([SVector(R*cos(2π*i/n), R*sin(2π*i/n)) for i in 0:n-1])

@testset "DiscreteCurves.jl" begin

# ────────────────────────────────────────────────────────────────────────────
@testset "Construction — explicit" begin
    pts = [SVector(0.0,0.0), SVector(1.0,0.0), SVector(1.0,1.0)]
    c = DiscreteCurve(pts)
    @test nvertices(c) == 3
    @test nedges(c)    == 2
    @test isopen(c)
    @test !isclosed(c)

    c2 = ClosedCurve(pts)
    @test isclosed(c2)
    @test nedges(c2) == 3

    # From tuples
    c3 = DiscreteCurve([(0.0,0.0),(1.0,0.0),(2.0,1.0)])
    @test nvertices(c3) == 3

    # From Matrix (columns = points)
    M = Float64[0 1 2; 0 1 0]
    c4 = DiscreteCurve(M)
    @test nvertices(c4) == 3

    # 3D
    c5 = OpenCurve([SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(0.0,0.0,1.0)])
    @test nvertices(c5) == 3

    @test_throws ArgumentError DiscreteCurve([SVector(1.0,2.0)])  # too few vertices
end

# ────────────────────────────────────────────────────────────────────────────
@testset "ParametricCurve — functional construction" begin
    # Basic 2D circle
    c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)
    @test nvertices(c) == 128
    @test isclosed(c)
    @test arc_length(c) ≈ 2π*sin(π/128)*128  atol=1e-4

    # From tuple-returning function
    c2 = ParametricCurve(t -> (cos(t), sin(t)), 0, 2π, 64; closed=true)
    @test nvertices(c2) == 64

    # 3D helix
    c3 = ParametricCurve(t -> SVector(cos(t), sin(t), t/(2π)), 0, 2π, 128)
    @test nvertices(c3) == 128

    # NTuple of functions
    c4 = ParametricCurve((sin, cos), 0, 2π, 64; closed=true)
    @test nvertices(c4) == 64
    @test vertex(c4, 1) ≈ SVector(sin(0.0), cos(0.0))  atol=1e-10

    # Step-based constructor
    c5 = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π; step=0.1)
    @test nvertices(c5) > 2

    # Float32 precision
    c6 = ParametricCurve(t -> SVector(cos(t), sin(t)), 0f0, 2f0*Float32(π), 64;
                          closed=true, T=Float32)
    @test eltype(c6.points[1]) == Float32
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Topology" begin
    c = ngon(6)

    # Modular indexing
    @test vertex(c, 7)  ≈ vertex(c, 1)
    @test vertex(c, 0)  ≈ vertex(c, 6)
    @test vertex(c, 13) ≈ vertex(c, 1)

    # Open curve: no wraparound
    c_open = OpenCurve(c.points)
    @test c_open[1] ≈ c.points[1]
    @test c_open[6] ≈ c.points[6]

    @test nedges(c)       == 6
    @test nedges(c_open)  == 5
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Lengths" begin
    sq = ClosedCurve([(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)])
    @test arc_length(sq) ≈ 4.0

    ℓ = edge_lengths(sq)
    @test length(ℓ) == 4
    @test all(ℓ .≈ 1.0)

    s = cumulative_length(sq)
    @test s[1] ≈ 0.0
    @test s[end] ≈ 4.0
    @test length(s) == 4
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Curvature" begin
    n = 256
    R = 2.0
    c = ParametricCurve(t -> SVector(R*cos(t), R*sin(t)), 0, 2π, n; closed=true)

    for model in (TurningAngle(), OsculatingCircle(), SteinerCurvature())
        κ = curvatures(c, model)
        @test length(κ) == n
        @test all(abs.(κ .- 1/R) .< 0.01)
    end

    # Gauss–Bonnet: ∑θ = 2π for a closed curve winding once
    θ = turning_angles(c)
    @test sum(θ) ≈ 2π  atol=0.01

    # Straight line → curvature = 0
    line = OpenCurve([SVector(i*1.0, 0.0) for i in 0:9])
    @test all(abs.(curvatures(line)) .< 1e-10)
end

# ────────────────────────────────────────────────────────────────────────────
@testset "CurveData — VertexData / EdgeData" begin
    c = ngon(8)
    temps = rand(8)

    vd = attach_vertex_data(c, temps)
    @test vd[1] ≈ temps[1]
    @test vd[9] ≈ temps[1]   # mod wraparound

    # From function
    vd2 = attach_vertex_data(c, p -> norm(p))
    @test all(abs.(collect(vd2) .- 1.0) .< 1e-10)

    # EdgeData
    ed = attach_edge_data(c, (p,q) -> norm(q-p))
    @test length(ed) == nedges(c)

    # Broadcasting works
    @test (vd2 .* 2) == collect(vd2) .* 2

    # Iteration
    total = sum(v for v in vd2)
    @test total ≈ 8.0  atol=1e-10

    # zip with vertices
    r_sum = sum(norm(p) + d for (p,d) in zip(vertices(c), vd2))
    @test r_sum ≈ 16.0  atol=1e-10

    @test_throws DimensionMismatch attach_vertex_data(c, rand(5))
end

# ────────────────────────────────────────────────────────────────────────────
@testset "LazyVertexField" begin
    c = ngon(8)

    lf = lazy_vertex_field(c, p -> norm(p))
    @test lf[1] ≈ 1.0  atol=1e-10
    @test length(lf) == 8

    # Iteration
    @test all(abs.(collect(lf) .- 1.0) .< 1e-10)

    # Materialise
    vd = collect(lf)
    @test vd isa VertexData
    @test length(vd) == 8
end

# ────────────────────────────────────────────────────────────────────────────
@testset "CurveDataBundle" begin
    c = ngon(8)
    b = bundle(c)

    b[:curvature]   = curvatures(c)
    b[:temperature] = rand(8)
    b[:edge_len]    = edge_lengths(c)   # edge channel (8 edges for closed octagon)
    b[:norm]        = p -> norm(p)      # lazy channel

    # Access styles
    @test b[:curvature, 1] ≈ b[:curvature][1]
    @test vertex_data(b, 1, :curvature) ≈ b[:curvature, 1]
    @test b.curvature[1]   ≈ b[:curvature, 1]

    # Channel types
    @test channel_type(b, :curvature)   == :vertex
    @test channel_type(b, :edge_len)    == :edge
    @test channel_type(b, :norm)        == :lazy

    # snapshot
    s = snapshot(b, 1)
    @test haskey(s, :curvature)
    @test haskey(s, :norm)
    @test !haskey(s, :edge_len)   # edges excluded from snapshot

    # Iteration over bundle channels
    channel_count = 0
    for (k, ch) in b
        channel_count += 1
    end
    @test channel_count == 4

    # Extract channel for tight loop (no Dict overhead per vertex)
    temp_ch = b[:temperature]
    total = sum(temp_ch[i] for i in 1:nvertices(c))
    @test total > 0  # sanity
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Iterators" begin
    c = ngon(5)

    # vertex iterator: length and values
    vs = collect(vertices(c))
    @test length(vs) == 5
    @test vs[1] ≈ c.points[1]

    # edge iterator: length, last edge wraps (closed)
    es = collect(edges(c))
    @test length(es) == 5
    @test es[end].src ≈ c.points[5]
    @test es[end].dst ≈ c.points[1]   # closing edge

    # arc_length via edge iterator (zero alloc path)
    L = sum(norm(e.dst - e.src) for e in edges(c))
    @test L ≈ arc_length(c)

    # window iterator (k=3): closed → 5 windows
    ws = collect(edge_windows(c; k=3))
    @test length(ws) == 5

    # window iterator: open → interior only
    c_open = OpenCurve(c.points)
    ws_open = collect(edge_windows(c_open; k=3))
    @test length(ws_open) == 3

    # Each window is the right size
    @test length(ws[1]) == 3

    # Co-iterate vertex + data (zip)
    temps = attach_vertex_data(c, rand(5))
    pairs = [(p, d) for (p,d) in zip(vertices(c), temps)]
    @test length(pairs) == 5
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Macros — @curve" begin
    # Basic 2D circle via macro
    c = @curve (cos(t), sin(t))  t ∈ 0..2π  n=128  closed
    @test nvertices(c) == 128
    @test isclosed(c)
    @test arc_length(c) ≈ 2π  atol=0.01

    # 3D helix
    c2 = @curve (cos(t), sin(t), t/(2π))  t ∈ 0..2π  n=64
    @test nvertices(c2) == 64

    # Custom variable name
    r = 2.0
    c3 = @curve (r*cos(θ), r*sin(θ))  θ ∈ 0..2π  n=32  closed
    @test nvertices(c3) == 32
    @test arc_length(c3) ≈ 2π*r  atol=0.1
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Macros — @where, @any_vertex, @all_vertices, @map_verts" begin
    c = @curve (cos(t), sin(t))  t ∈ 0..2π  n=128  closed

    # Filter: upper half-plane
    upper = @where c p  p[2] > 0
    @test all(p[2] > 0 for p in upper)
    @test length(upper) < nvertices(c)

    # Existential
    @test @any_vertex c p  p[2] > 0.9
    @test !(@any_vertex c p  norm(p) > 2.0)

    # Universal
    @test @all_vertices c p  abs(norm(p) - 1.0) < 1e-3

    # Map: project each vertex onto unit sphere (already unit circle here)
    c2 = @map_verts c p  p / norm(p)
    @test all(abs(norm(vertex(c2,i)) - 1.0) < 1e-10 for i in 1:nvertices(c))
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Macros — @functional" begin
    @functional avg_radius c = sum(norm(p) for p in vertices(c)) / nvertices(c)

    c = @curve (cos(t), sin(t))  t ∈ 0..2π  n=128  closed
    @test abs(avg_radius(c) - 1.0) < 1e-3

    @functional total_turning c = sum(turning_angles(c))
    @test abs(total_turning(c) - 2π) < 0.01
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Utils" begin
    c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 64; closed=true)

    c2 = translate(c, SVector(1.0, 0.0))
    @test vertex(c2,1) ≈ vertex(c,1) + SVector(1.0,0.0)

    c3 = scale_curve(c, 3.0)
    @test norm(vertex(c3,1)) ≈ 3.0  atol=1e-10

    c4 = map_vertices(p -> SVector(-p[2], p[1]), c)  # 90° rotation
    @test norm(vertex(c4,1)) ≈ 1.0  atol=1e-10

    c5 = resample(c, 200)
    @test nvertices(c5) == 200
    @test arc_length(c5) ≈ arc_length(c)  atol=1e-4
end

# ────────────────────────────────────────────────────────────────────────────
@testset "Type stability" begin
    c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 64; closed=true)
    @inferred vertex(c, 1)
    @inferred arc_length(c)
    @inferred edge_lengths(c)
    @inferred tangent(c, 1)
    @inferred unit_tangent(c, 1)
    @inferred curvature(c, 2)
    @inferred normal(c, 2)
end

end # @testset
println("\n✓  All tests passed.")