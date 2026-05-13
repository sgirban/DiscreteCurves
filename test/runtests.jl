# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  test/runtests.jl - DiscreteCurves Unit Test Suite                          ║
# ║  Comprehensive, well-organized test suite for DiscreteCurves.jl             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

using Test
using DiscreteCurves
using StaticArrays
using LinearAlgebra
using Statistics

# ──────────────────────────────────────────────────────────────────────────────
# TEST HELPERS
# ──────────────────────────────────────────────────────────────────────────────

"""Regular n-gon inscribed in circle of radius R."""
ngon(n, R=1.0) = ClosedCurve([SVector(R*cos(2π*i/n), R*sin(2π*i/n)) for i in 0:n-1])

"""Create a simple square."""
make_square(closed=true) = ClosedCurve([
    SVector(0.0, 0.0),
    SVector(1.0, 0.0),
    SVector(1.0, 1.0),
    SVector(0.0, 1.0),
])

"""Create a simple line segment."""
make_line(n=11) = OpenCurve([SVector(Float64(i-1)/(n-1), 0.0) for i in 1:n])

"""Create a unit circle."""
make_circle(n=64) = ClosedCurve([SVector(cos(2π*i/n), sin(2π*i/n)) for i in 0:n-1])

# ──────────────────────────────────────────────────────────────────────────────
# MAIN TEST SUITE
# ──────────────────────────────────────────────────────────────────────────────

@testset "DiscreteCurves.jl" begin

# ────────────────────────────────────────────────────────────────────────────
@testset "1. CONSTRUCTION" begin
    @testset "1.1 Explicit construction from points" begin
        pts = [SVector(0.0,0.0), SVector(1.0,0.0), SVector(1.0,1.0)]
        c = DiscreteCurve(pts)
        @test nvertices(c) == 3
        @test nedges(c)    == 2
        @test isopen(c)
        @test !isclosed(c)
    end
    
    @testset "1.2 Open vs closed curves" begin
        pts = [SVector(0.0,0.0), SVector(1.0,0.0), SVector(1.0,1.0)]
        
        c_open = OpenCurve(pts)
        @test nedges(c_open) == 2
        @test isopen(c_open)
        
        c_closed = ClosedCurve(pts)
        @test nedges(c_closed) == 3
        @test isclosed(c_closed)
    end

    @testset "1.3 Construction from different input formats" begin
        # From tuples
        c1 = DiscreteCurve([(0.0,0.0),(1.0,0.0),(2.0,1.0)])
        @test nvertices(c1) == 3

        # From matrix (columns = points)
        M = Float64[0 1 2; 0 1 0]
        c2 = DiscreteCurve(M)
        @test nvertices(c2) == 3
        @test c2[1][1] ≈ 0.0
        @test c2[2][1] ≈ 1.0
    end

    @testset "1.4 3D curves" begin
        pts_3d = [SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(0.0,0.0,1.0)]
        c = OpenCurve(pts_3d)
        @test nvertices(c) == 3
        @test length(c[1]) == 3
    end

    @testset "1.5 Error handling" begin
        @test_throws ArgumentError DiscreteCurve([SVector(1.0,2.0)])  # too few vertices
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "1.6 Parametric curve construction" begin
    @testset "1.6.1 Basic 2D circle from parametric function" begin
        c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)
        @test nvertices(c) == 128
        @test isclosed(c)
        @test arc_length(c) ≈ 2π*sin(π/128)*128  atol=1e-4
    end

    @testset "1.6.2 Tuple-returning function" begin
        c = ParametricCurve(t -> (cos(t), sin(t)), 0, 2π, 64; closed=true)
        @test nvertices(c) == 64
    end

    @testset "1.6.3 3D helix" begin
        c = ParametricCurve(t -> SVector(cos(t), sin(t), t/(2π)), 0, 2π, 128)
        @test nvertices(c) == 128
        @test length(c[1]) == 3
    end

    @testset "1.6.4 Function tuple (sin, cos)" begin
        c = ParametricCurve((sin, cos), 0, 2π, 64; closed=true)
        @test nvertices(c) == 64
        @test vertex(c, 1) ≈ SVector(sin(0.0), cos(0.0))  atol=1e-10
    end

    @testset "1.6.5 Step-based constructor" begin
        c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π; step=0.1)
        @test nvertices(c) > 2
    end

    @testset "1.6.6 Float32 precision" begin
        c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0f0, 2f0*Float32(π), 64;
                              closed=true, T=Float32)
        @test eltype(c.points[1]) == Float32
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "2. TOPOLOGY" begin
    @testset "2.1 Modular indexing (closed curves)" begin
        c = ngon(6)
        
        # Wraparound indexing
        @test vertex(c, 7)  ≈ vertex(c, 1)
        @test vertex(c, 0)  ≈ vertex(c, 6)
        @test vertex(c, 13) ≈ vertex(c, 1)
        @test vertex(c, -1) ≈ vertex(c, 5)
    end

    @testset "2.2 Edge counting (open vs closed)" begin
        pts = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0)]
        
        c_open = OpenCurve(pts)
        @test nedges(c_open) == 2
        
        c_closed = ClosedCurve(pts)
        @test nedges(c_closed) == 3
    end

    @testset "2.3 Open curve: no wraparound" begin
        c_open = make_line()
        @test c_open[1] ≈ c_open.points[1]
        @test c_open[nvertices(c_open)] ≈ c_open.points[end]
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "3. LENGTHS & ARC LENGTH" begin
    @testset "3.1 Unit square arc length" begin
        sq = make_square()
        @test arc_length(sq) ≈ 4.0
    end

    @testset "3.2 Edge lengths" begin
        sq = make_square()
        ℓ = edge_lengths(sq)
        @test length(ℓ) == 4
        @test all(ℓ .≈ 1.0)
    end

    @testset "3.3 Cumulative length" begin
        sq = make_square()
        s = cumulative_length(sq)
        @test s[1] ≈ 0.0
        @test s[end] ≈ 4.0
        @test length(s) == 4
        # Cumulative should be monotonically increasing
        @test all(diff(s) .>= 0)
    end

    @testset "3.4 Open vs closed arc length" begin
        c_open = make_line()
        c_closed = ClosedCurve(c_open.points)
        
        # Closed adds one more edge
        @test arc_length(c_closed) > arc_length(c_open)
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "4. CURVATURE" begin
    @testset "4.1 Curvature of a circle" begin
        n = 256
        R = 2.0
        c = ParametricCurve(t -> SVector(R*cos(t), R*sin(t)), 0, 2π, n; closed=true)

        for model in (TurningAngle(), OsculatingCircle(), SteinerCurvature())
            κ = curvatures(c, model)
            @test length(κ) == n
            # For circle: curvature = 1/R everywhere
            @test all(abs.(κ .- 1/R) .< 0.01)
        end
    end

    @testset "4.2 Turning angles sum to 2π (Gauss-Bonnet)" begin
        c = make_circle(128)
        θ = turning_angles(c)
        @test sum(θ) ≈ 2π  atol=0.05
    end

    @testset "4.3 Straight line has zero curvature" begin
        line = make_line(20)
        κ = curvatures(line, SteinerCurvature())
        @test all(abs.(κ) .< 1e-10)
    end

    @testset "4.4 Curvature vector properties" begin
        c = make_circle(64)
        κv = curvature_vectors(c, SteinerCurvature())
        @test length(κv) == nvertices(c)
        
        # For a circle, curvature vectors should point inward
        for i in 1:nvertices(c)
            v = c[i]
            κvec = κv[i]
            # Dot product with position should be negative (pointing in)
            @test dot(κvec, v) < 0
        end
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "5. GEOMETRIC PROPERTIES" begin
    @testset "5.1 Signed area (orientation)" begin
        # CCW square
        sq_ccw = make_square()
        area_ccw = signed_area(sq_ccw)
        @test area_ccw ≈ 1.0
        
        # CW square (reversed vertices)
        sq_cw = ClosedCurve(reverse(sq_ccw.points))
        area_cw = signed_area(sq_cw)
        @test area_cw ≈ -1.0
    end

    @testset "5.2 Bounding box" begin
        c = make_square()
        bbox = bounding_box(c)
        @test bbox.min ≈ SVector(0.0, 0.0)
        @test bbox.max ≈ SVector(1.0, 1.0)
    end

    @testset "5.3 Centroid" begin
        sq = make_square()
        ctr = centroid(sq)
        @test ctr ≈ SVector(0.5, 0.5)  atol=0.01
    end

    @testset "5.4 Isoperimetric quotient" begin
        # For a circle, isoperimetric quotient = 1
        circle = make_circle(256)
        q = isoperimetric_quotient(circle)
        @test 0.95 < q ≤ 1.0
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "6. TANGENT & NORMAL VECTORS" begin
    @testset "6.1 Tangent vectors (unit circle)" begin
        c = make_circle(64)
        tvecs = tangent_vectors(c)
        
        @test length(tvecs) == nvertices(c)
        # All should be unit length (approximately)
        @test all(abs.(norm.(tvecs) .- 1.0) .< 0.1)
    end

    @testset "6.2 Normal vectors (orthogonal to tangent)" begin
        c = make_circle(64)
        tvecs = tangent_vectors(c)
        nvecs = normal_vectors(c)
        
        @test length(nvecs) == nvertices(c)
        # Tangent and normal should be orthogonal
        for i in 1:nvertices(c)
            @test abs(dot(tvecs[i], nvecs[i])) < 1e-10
        end
    end

    @testset "6.3 Inward vs outward normals" begin
        sq = make_square()
        n_in = inward_normals(sq)
        n_out = outward_normals(sq)
        
        @test length(n_in) == 4
        @test length(n_out) == 4
        
        # They should point in opposite directions
        for i in 1:4
            @test n_in[i] ≈ -n_out[i]  atol=1e-10
        end
    end

    @testset "6.4 Normals point toward/away from centroid" begin
        c = make_circle(64)
        n_in = inward_normals(c)
        ctr = centroid(c)
        
        # Inward normals should have positive dot product with centroid direction
        for i in 1:nvertices(c)
            to_center = ctr - c[i]
            @test dot(n_in[i], to_center) > 0
        end
    end
end

# ────────────────────────────────────────────────────────────────────────────
@testset "7. DATA ATTACHMENT" begin
    @testset "7.1 Vertex data from array" begin
        c = ngon(8)
        temps = rand(8)

        vd = attach_vertex_data(c, temps)
        @test vd[1] ≈ temps[1]
        @test vd[9] ≈ temps[1]   # wraparound for closed curves
    end

    @testset "7.2 Vertex data from function" begin
        c = ngon(8)
        vd = attach_vertex_data(c, p -> norm(p))
        @test all(abs.(collect(vd) .- 1.0) .< 1e-10)
    end

    @testset "7.3 Edge data" begin
        c = ngon(8)
        ed = attach_edge_data(c, (p,q) -> norm(q-p))
        @test length(ed) == nedges(c)
    end

    @testset "7.4 Vertex data broadcasting" begin
        c = ngon(8)
        vd = attach_vertex_data(c, p -> norm(p))
        @test (vd .* 2) == collect(vd) .* 2
    end

    @testset "7.5 Vertex data iteration" begin
        c = ngon(8)
        vd = attach_vertex_data(c, p -> norm(p))
        total = sum(v for v in vd)
        @test total ≈ 8.0  atol=1e-10
    end

    @testset "7.6 Vertex data zip with vertices" begin
        c = ngon(8)
        vd = attach_vertex_data(c, p -> norm(p))
        r_sum = sum(norm(p) + d for (p,d) in zip(vertices(c), vd))
        @test r_sum ≈ 16.0  atol=1e-10
    end

    @testset "7.7 Vertex data dimension checking" begin
        c = ngon(8)
        @test_throws DimensionMismatch attach_vertex_data(c, rand(5))
    end

    @testset "7.8 Lazy vertex field" begin
        c = ngon(8)
        lf = lazy_vertex_field(c, p -> norm(p))
        
        @test lf[1] ≈ 1.0  atol=1e-10
        @test length(lf) == 8
        @test all(abs.(collect(lf) .- 1.0) .< 1e-10)
        
        # Materialize
        vd = collect(lf)
        @test vd isa VertexData
    end
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