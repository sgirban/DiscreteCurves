module Generators
using LinearAlgebra
using StaticArrays
using ..CurveTypes
using ..AbstractTypes
using ..CurveTopology
using Random

export regular_polygon, inscribed_polygon, circle_curve, ellipse_curve, rectangle_curve,
       square_curve, triangle_curve, equilateral_triangle, right_triangle,
       parallelogram_curve, rhombus_curve, star_curve,
       cardioid, limacon, rose_curve, lemniscate, epicycloid, hypocycloid,
       astroid, deltoid, nephroid, trefoil_knot, figure_eight, heart_curve,
       kidney_curve, butterfly_curve, spiral_curve, logarithmic_spiral
export random_convex_polygon, random_concave_polygon, fourier_curve, add_noise,
       rough_curve, noisy_polygon, perlin_curve
# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers
# ─────────────────────────────────────────────────────────────────────────────
@inline _angles(n::Int) = range(0.0, 2π; length=n+1)[1:n]
 
function _transform(pts::Vector{SVector{2,T}},
                    center::SVector{2,T},
                    angle::Real) where T
    cosα, sinα = T(cos(angle)), T(sin(angle))
    [SVector{2,T}(cosα*p[1] - sinα*p[2] + center[1],
                  sinα*p[1] + cosα*p[2] + center[2]) for p in pts]
end

function _sample_closed(f, n::Int, t0::Real=0.0, t1::Real=2π, ::Type{T}=Float64) where T
    ts  = range(T(t0), T(t1); length=n+1)[1:n]
    pts = [SVector{2,T}(f(t)) for t in ts]
    ClosedCurve(pts)
end

"""
    regular_polygon(n; center, radius, angle, T)
 
Regular `n`-gon inscribed in a circle of the given `radius`.
 
```julia
regular_polygon(6)                          # unit hexagon
regular_polygon(4; radius=2.0)             # square, side = √2 × 2
regular_polygon(5; center=SVector(1,1), angle=π/10) # regular pentagon, rotated by π/10, centered at (1,1)
```
"""
function regular_polygon(n::Int;
                          center :: SVector{2,Float64} = SVector(0.0, 0.0),
                          radius :: Real               = 1.0,
                          angle  :: Real               = 0.0,
                          T      :: Type{<:AbstractFloat} = Float64)
    n ≥ 3 || throw(ArgumentError("n must be ≥ 3, got $n"))
    r = T(radius)
    pts = [SVector{2,T}(r*cos(2π*i/n), r*sin(2π*i/n)) for i in 0:n-1]
    ClosedCurve(_transform(pts, SVector{2,T}(center), angle))
end

"""
    inscribed_polygon(n; center, radius, irregularity, T)
 
`n` points sampled on a circle with controllable irregularity.
 
- `irregularity=0.0` — equally spaced (= `regular_polygon`).
- `irregularity=1.0` — fully random spacing (Voronoi-style random polygon).
 
```julia
inscribed_polygon(8)                        # regular octagon
inscribed_polygon(8; irregularity=0.5)     # slightly wobbly
inscribed_polygon(20; irregularity=1.0)    # fully random
```
"""
function inscribed_polygon(n::Int;
                            center      :: SVector{2,Float64} = SVector(0.0, 0.0),
                            radius      :: Real               = 1.0,
                            irregularity:: Real               = 0.0,
                            seed        :: Int                = 0,
                            T           :: Type{<:AbstractFloat} = Float64)
    n ≥ 3 || throw(ArgumentError("n must be ≥ 3, got $n"))
    rng   = _seeded_rng(seed)
    irr   = clamp(T(irregularity), zero(T), one(T))
    step  = 2π / n
    max_δ = irr * step / 2
    angles = T[2π*i/n + max_δ*(2rand(rng, T)-1) for i in 0:n-1]
    sort!(angles)
    r  = T(radius)
    pts = [SVector{2,T}(r*cos(θ), r*sin(θ)) for θ in angles]
    ClosedCurve(_transform(pts, SVector{2,T}(center), 0.0))
end

# ─────────────────────────────────────────────────────────────────────────────
#  Standard shapes
# ─────────────────────────────────────────────────────────────────────────────
 
"""
    circle_curve(n=128; center, radius, T)
 
Uniformly sampled circle.
"""
circle_curve(n::Int=128;
             center=SVector(0.0,0.0), radius::Real=1.0, T::Type{<:AbstractFloat}=Float64) =
    regular_polygon(n; center, radius, T)
"""
    ellipse_curve(n=128; a, b, center, angle, T)
 
Ellipse with semi-axes `a` (horizontal) and `b` (vertical).
"""
function ellipse_curve(n::Int=128;
                        a      :: Real  = 1.0,
                        b      :: Real  = 0.5,
                        center         = SVector(0.0,0.0),
                        angle  :: Real  = 0.0,
                        T      :: Type{<:AbstractFloat} = Float64)
    _sample_closed(n, T) do t
        SVector{2,T}(T(a)*cos(t), T(b)*sin(t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), angle))
end

"""
    rectangle_curve(; width, height, center, angle, n, T)
 
Axis-aligned rectangle with `n` total vertices distributed proportionally.
"""
function rectangle_curve(;
                          width  :: Real  = 2.0,
                          height :: Real  = 1.0,
                          center         = SVector(0.0,0.0),
                          angle  :: Real  = 0.0,
                          n      :: Int   = 64,
                          T      :: Type{<:AbstractFloat} = Float64)
    w, h   = T(width)/2, T(height)/2
    perim  = 2*(width + height)
    n_h = max(2, round(Int, n * width / perim))
    n_v = max(2, round(Int, n * height / perim))
 
    pts = SVector{2,T}[]
    for x in range(-w, w; length=n_h+1)[1:end-1]; push!(pts, SVector(-w+x+w, -h)); end  # nope
    pts = SVector{2,T}[]
    for x in range(-w,  w; length=n_h+1)[1:end-1]; push!(pts, SVector{2,T}( x, -h)); end
    for y in range(-h,  h; length=n_v+1)[1:end-1]; push!(pts, SVector{2,T}( w,  y)); end
    for x in range( w, -w; length=n_h+1)[1:end-1]; push!(pts, SVector{2,T}( x,  h)); end
    for y in range( h, -h; length=n_v+1)[1:end-1]; push!(pts, SVector{2,T}(-w,  y)); end
    ClosedCurve(_transform(pts, SVector{2,T}(center), angle))
end
 
"""    square_curve(; side=2.0, kwargs...)  — axis-aligned square."""
square_curve(; side::Real=2.0, kw...) = rectangle_curve(; width=side, height=side, kw...)
 
"""
    triangle_curve(p1, p2, p3; n, T)
 
Triangle from three vertex positions. `n` points distributed proportionally.
"""
function triangle_curve(p1, p2, p3;
                         n :: Int = 48,
                         T :: Type{<:AbstractFloat} = Float64)
    v1 = SVector{2,T}(p1);  v2 = SVector{2,T}(p2);  v3 = SVector{2,T}(p3)
    ℓ12 = norm(v2-v1);  ℓ23 = norm(v3-v2);  ℓ31 = norm(v1-v3)
    perim = ℓ12 + ℓ23 + ℓ31
    n12 = max(2, round(Int, n*ℓ12/perim))
    n23 = max(2, round(Int, n*ℓ23/perim))
    n31 = n - n12 - n23
 
    pts = SVector{2,T}[]
    for t in range(0,1; length=n12+1)[1:end-1]; push!(pts, (1-t)*v1 + t*v2); end
    for t in range(0,1; length=n23+1)[1:end-1]; push!(pts, (1-t)*v2 + t*v3); end
    for t in range(0,1; length=n31+1)[1:end-1]; push!(pts, (1-t)*v3 + t*v1); end
    ClosedCurve(pts)
end
 
"""    equilateral_triangle(; side=2.0, center, angle, n, T)"""
function equilateral_triangle(; side::Real=2.0, center=SVector(0.0,0.0),
                                angle::Real=π/2, n::Int=48, T::Type{<:AbstractFloat}=Float64)
    r = T(side) / sqrt(T(3))
    p1 = SVector{2,T}(r*cos(π/2),      r*sin(π/2))
    p2 = SVector{2,T}(r*cos(π/2+2π/3), r*sin(π/2+2π/3))
    p3 = SVector{2,T}(r*cos(π/2+4π/3), r*sin(π/2+4π/3))
    c  = triangle_curve(p1, p2, p3; n, T)
    ClosedCurve(_transform(c.points, SVector{2,T}(center), angle))
end
 
"""    right_triangle(; a=1.0, b=1.0, n=48, T)  — right angle at origin."""
right_triangle(; a::Real=1.0, b::Real=1.0, n::Int=48, T::Type{<:AbstractFloat}=Float64) =
    triangle_curve(SVector(0,0), SVector(a,0), SVector(0,b); n, T)
 
"""
    parallelogram_curve(; a=2.0, b=1.0, skew=π/4, n=64, T)
 
Parallelogram with base `a`, side `b`, interior angle `skew` (radians).
"""
function parallelogram_curve(; a::Real=2.0, b::Real=1.0, skew::Real=π/4,
                               center=SVector(0.0,0.0), n::Int=64,
                               T::Type{<:AbstractFloat}=Float64)
    A, B, sk = T(a), T(b), T(skew)
    v1 = SVector{2,T}(zero(T), zero(T))
    v2 = SVector{2,T}(A, zero(T))
    v3 = SVector{2,T}(A + B*cos(sk), B*sin(sk))
    v4 = SVector{2,T}(B*cos(sk), B*sin(sk))
    ctr = (v1+v2+v3+v4)/4
    pts = SVector{2,T}[v1-ctr, v2-ctr, v3-ctr, v4-ctr]
    out = SVector{2,T}[]
    for (va, vb) in [(pts[1],pts[2]),(pts[2],pts[3]),(pts[3],pts[4]),(pts[4],pts[1])]
        np = max(2, n ÷ 4)
        for t in range(0,1; length=np+1)[1:end-1]; push!(out, (1-t)*va + t*vb); end
    end
    ClosedCurve(_transform(out, SVector{2,T}(center), 0.0))
end
 
"""    rhombus_curve(; side=1.0, angle=π/3, n=64, T)"""
rhombus_curve(; side::Real=1.0, angle::Real=π/3, kw...) =
    parallelogram_curve(; a=side, b=side, skew=angle, kw...)
 
"""
    star_curve(n_points; inner_r, outer_r, center, angle, T)
 
`n_points`-pointed star polygon.
 
```julia
star_curve(5)                          # 5-pointed star
star_curve(6; inner_r=0.4)            # Star of David shape
```
"""
function star_curve(n_points::Int=5;
                    inner_r :: Real  = 0.4,
                    outer_r :: Real  = 1.0,
                    center          = SVector(0.0,0.0),
                    angle   :: Real  = -π/2,
                    T       :: Type{<:AbstractFloat} = Float64)
    n_points ≥ 3 || throw(ArgumentError("need ≥ 3 points"))
    ri, ro = T(inner_r), T(outer_r)
    pts = SVector{2,T}[]
    for k in 0:2n_points-1
        r = iseven(k) ? ro : ri
        θ = T(angle) + k*π/n_points
        push!(pts, SVector{2,T}(r*cos(θ), r*sin(θ)))
    end
    ClosedCurve(_transform(pts, SVector{2,T}(center), 0.0))
end

# ─────────────────────────────────────────────────────────────────────────────
#  Parametric curves zoo
# ─────────────────────────────────────────────────────────────────────────────

"""
    cardioid(n=256; a=1.0, center, angle, T)
 
Heart-shaped curve: r(θ) = a(1 − cos θ) in polar.
"""
cardioid(n::Int=256; a::Real=1.0, center=SVector(0.0,0.0),
          angle::Real=0.0, T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        r = T(a)*(1 - cos(t))
        SVector{2,T}(r*cos(t), r*sin(t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), angle))

"""
    limacon(n=256; a=0.5, b=1.0, center, T)
 
Limaçon of Pascal: r(θ) = a + b·cos θ.
- b > a: limaçon with inner loop
- b = a: cardioid
- b < a: convex limaçon
"""
limacon(n::Int=256; a::Real=0.5, b::Real=1.0,
         center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        r = T(a) + T(b)*cos(t)
        SVector{2,T}(r*cos(t), r*sin(t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))

"""
    rose_curve(k, n=512; a=1.0, center, T)
 
Rose curve: r(θ) = a·cos(k·θ).
- Odd k: k petals.
- Even k: 2k petals.
Range: [0, 2π] for odd k, [0, 4π] for even k (closes after double period).
"""
function rose_curve(k::Int, n::Int=512; a::Real=1.0,
                    center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    period = iseven(k) ? 4π : 2π
    _sample_closed(n, T, 0.0, period) do t
        r = T(a)*cos(k*t)
        SVector{2,T}(r*cos(t), r*sin(t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
end
 
"""
    lemniscate(n=256; a=1.0, center, T)
 
Lemniscate of Bernoulli: figure-eight shaped curve.
Polar: r² = a²·cos(2θ). Only real for cos(2θ) ≥ 0.
"""
function lemniscate(n::Int=256; a::Real=1.0,
                    center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    A = T(a)
    _sample_closed(n, T) do t
        denom = 1 + sin(t)^2
        SVector{2,T}(A*cos(t)/denom, A*sin(t)*cos(t)/denom)
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
end

"""
    epicycloid(k, n=512; R=1.0, r, center, T)
 
Epicycloid: small circle of radius `r` rolling outside a circle of radius `R`.
`k` cusps when r = R/k.
"""
function epicycloid(k::Int=3, n::Int=512;
                    R::Real=1.0, r::Real=0.0,
                    center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    Rv, rv = T(R), r == 0.0 ? T(R)/k : T(r)
    _sample_closed(n, T) do t
        SVector{2,T}(
            (Rv+rv)*cos(t) - rv*cos((Rv+rv)/rv*t),
            (Rv+rv)*sin(t) - rv*sin((Rv+rv)/rv*t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
end
 
"""
    hypocycloid(k, n=512; R=1.0, r, center, T)
 
Hypocycloid: small circle of radius `r` rolling inside circle of radius `R`.
k=3: deltoid, k=4: astroid.
"""
function hypocycloid(k::Int=3, n::Int=512;
                     R::Real=1.0, r::Real=0.0,
                     center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    Rv, rv = T(R), r == 0.0 ? T(R)/k : T(r)
    _sample_closed(n, T) do t
        SVector{2,T}(
            (Rv-rv)*cos(t) + rv*cos((Rv-rv)/rv*t),
            (Rv-rv)*sin(t) - rv*sin((Rv-rv)/rv*t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
end
 
"""    astroid(n=256; a=1.0, center, T)  — 4-cusped hypocycloid: x^(2/3)+y^(2/3)=a^(2/3)"""
astroid(n::Int=256; a::Real=1.0, center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        SVector{2,T}(T(a)*cos(t)^3, T(a)*sin(t)^3)
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
 
"""    deltoid(n=256; R=1.0, T)  — 3-cusped hypocycloid"""
deltoid(n::Int=256; R::Real=1.0, T::Type{<:AbstractFloat}=Float64) = hypocycloid(3, n; R, T)
 
"""    nephroid(n=256; a=1.0, T)  — 2-cusped epicycloid (kidney-shaped)"""
nephroid(n::Int=256; a::Real=1.0, center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    epicycloid(2, n; R=T(a), center, T)
 
"""
    trefoil_knot(n=512; a=2.0, b=1.0, T)
 
Trefoil knot projected onto 2D (approximate — actual knot is 3D).
Parametric 2D trefoil:
  x(t) = sin(t) + 2sin(2t)
  y(t) = cos(t) − 2cos(2t)
"""
trefoil_knot(n::Int=512; a::Real=2.0, b::Real=1.0,
              center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        SVector{2,T}(sin(t)+T(a)*sin(2t), cos(t)-T(a)*cos(2t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
 
"""
    figure_eight(n=256; a=1.0, center, T)
 
Figure-eight (lemniscate of Gerono): x = a·sin(t), y = a·sin(t)cos(t).
"""
figure_eight(n::Int=256; a::Real=1.0, center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        SVector{2,T}(T(a)*sin(t), T(a)*sin(t)*cos(t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
 
"""
    heart_curve(n=256; scale=1.0, center, T)
 
Algebraic heart curve (the classic x²+(y-∛x²)²=1 approximation).
Parametric: x = sin³t, y = cos t − cos²t − cos³t (scaled).
"""
heart_curve(n::Int=256; scale::Real=1.0, center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        s = T(scale)
        SVector{2,T}(s*sin(t)^3,
                     s*(cos(t) - cos(2t)/2 - cos(3t)/8 - cos(4t)/16) * 0.85)
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
 
"""
    kidney_curve(n=256; a=1.0, center, T)
 
Kidney-shaped curve (limaçon with a=b, similar to cardioid family).
r(θ) = a(1 + 0.6·cos θ).
"""
kidney_curve(n::Int=256; a::Real=1.0, center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T) do t
        r = T(a)*(1 + T(0.6)*cos(t))
        SVector{2,T}(r*cos(t), r*sin(t))
    end |> c -> ClosedCurve(_transform(c.points, SVector{2,T}(center), 0.0))
 
"""
    butterfly_curve(n=512; a=1.0, T)
 
Butterfly curve (Temple H. Fay, 1989).
r = e^{sin θ} − 2cos(4θ) + sin⁵((2θ-π)/24)
"""
butterfly_curve(n::Int=1024; a::Real=1.0, T::Type{<:AbstractFloat}=Float64) =
    _sample_closed(n, T, 0.0, 4π) do t
        r = T(a)*(exp(sin(t)) - 2cos(4t) + sin((2t - π)/24)^5)
        SVector{2,T}(r*cos(t), r*sin(t))
    end
 
"""
    spiral_curve(turns=3, n=512; r_start=0.1, r_end=1.0, T)
 
Archimedean spiral (open curve).
"""
function spiral_curve(turns::Int=3, n::Int=512;
                       r_start::Real=0.1, r_end::Real=1.0,
                       center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    t_end = 2π*turns
    ts    = range(zero(T), T(t_end); length=n)
    pts   = [SVector{2,T}((T(r_start)+(T(r_end)-T(r_start))*t/T(t_end))*cos(t),
                           (T(r_start)+(T(r_end)-T(r_start))*t/T(t_end))*sin(t)) for t in ts]
    OpenCurve(_transform(pts, SVector{2,T}(center), 0.0))
end
 
"""
    logarithmic_spiral(turns=3, n=512; a=0.1, b=0.2, T)
 
Logarithmic (equiangular) spiral: r = a·e^{b·θ}.
"""
function logarithmic_spiral(turns::Int=3, n::Int=512;
                              a::Real=0.1, b::Real=0.2,
                              center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    t_end = 2π*turns
    ts    = range(zero(T), T(t_end); length=n)
    pts   = [SVector{2,T}(T(a)*exp(T(b)*t)*cos(t), T(a)*exp(T(b)*t)*sin(t)) for t in ts]
    OpenCurve(_transform(pts, SVector{2,T}(center), 0.0))
end
 
"""
    parabola_arc(n=128; a=1.0, x_range=(-2,2), T)
 
Open parabolic arc: y = a·x².
"""
function parabola_arc(n::Int=128; a::Real=1.0,
                       x_range::Tuple{Real,Real}=(-2.0,2.0),
                       center=SVector(0.0,0.0), T::Type{<:AbstractFloat}=Float64)
    xs  = range(T(x_range[1]), T(x_range[2]); length=n)
    pts = [SVector{2,T}(x, T(a)*x^2) for x in xs]
    OpenCurve(_transform(pts, SVector{2,T}(center), 0.0))
end

"""
    random_convex_polygon(n=16; radius=1.0, irregularity=0.5, spikiness=0.2, seed, T)
 
Random convex polygon via Valtr's algorithm (sorted random angle gaps).
 
- `irregularity` ∈ [0,1]: variation in angle spacing.
- `spikiness` ∈ [0,1]: variation in radius per vertex.
 
```julia
random_convex_polygon(20)
random_convex_polygon(8; irregularity=0.1, spikiness=0.0)   # almost regular
random_convex_polygon(50; irregularity=1.0, spikiness=0.5)  # very random
```
"""
function random_convex_polygon(n::Int=16;
                                radius      :: Real = 1.0,
                                irregularity:: Real = 0.5,
                                spikiness   :: Real = 0.2,
                                center              = SVector(0.0,0.0),
                                seed        :: Int  = 0,
                                T           :: Type{<:AbstractFloat} = Float64)
    n ≥ 3 || throw(ArgumentError("n ≥ 3 required"))
    rng  = _seeded_rng(seed)
    irr  = clamp(T(irregularity), zero(T), one(T))
    spk  = clamp(T(spikiness), zero(T), one(T))
    R    = T(radius)
 
    gaps = [rand(rng, T) + irr * rand(rng, T) for _ in 1:n]
    gaps .*= (2π / sum(gaps))
    angles = cumsum(gaps)
 
    pts = [SVector{2,T}((R + spk*R*(rand(rng,T)-T(0.5)))*cos(θ),
                         (R + spk*R*(rand(rng,T)-T(0.5)))*sin(θ)) for θ in angles]
    ClosedCurve(_transform(pts, SVector{2,T}(center), 0.0))
end
 
"""
    random_concave_polygon(n=16; radius=1.0, concavity=0.5, seed, T)
 
Random non-convex (concave) polygon: some vertices pushed inward.
 
- `concavity` ∈ [0,1]: fraction of vertices that are concave.
"""
function random_concave_polygon(n::Int=16;
                                 radius    :: Real = 1.0,
                                 concavity :: Real = 0.5,
                                 center            = SVector(0.0,0.0),
                                 seed      :: Int  = 0,
                                 T         :: Type{<:AbstractFloat} = Float64)
    rng = _seeded_rng(seed)
    R   = T(radius)
    cc  = clamp(T(concavity), zero(T), one(T))
 
    angles = sort!(rand(rng, T, n) .* 2π)
    pts = map(eachindex(angles)) do i
        θ = angles[i]
        r = rand(rng) < cc ? R * (T(0.3) + T(0.4)*rand(rng, T)) : R * (T(0.8) + T(0.2)*rand(rng, T))
        SVector{2,T}(r*cos(θ), r*sin(θ))
    end
    ClosedCurve(_transform(pts, SVector{2,T}(center), 0.0))
end

"""
    fourier_curve(n=256; modes=8, radius=1.0, roughness=0.5, seed, center, T)
 
Smooth random closed curve via random Fourier modes.
 
Each Fourier mode k has amplitude decaying as `1/k^α` where `α = 1/(roughness+ε)`.
- `roughness=0.0`: very smooth (low-frequency modes dominate).
- `roughness=1.0`: rougher, more complex boundary.
- `modes`: number of Fourier modes to include (more = richer shape).
 
```julia
fourier_curve()                           # smooth random blob
fourier_curve(256; modes=4, roughness=0.1)  # very smooth
fourier_curve(256; modes=20, roughness=0.9) # complex, many features
```
"""
function fourier_curve(n::Int=256;
                        modes     :: Int  = 8,
                        radius    :: Real = 1.0,
                        roughness :: Real = 0.5,
                        center           = SVector(0.0,0.0),
                        seed      :: Int  = 0,
                        T         :: Type{<:AbstractFloat} = Float64)
    rng = _seeded_rng(seed)
    R   = T(radius)
    α   = 1 / (T(roughness) + T(0.1))   # higher α = faster amplitude decay
 
   
    a_k = [R * rand(rng, T) / k^α for k in 1:modes]
    φ_k = [2π * rand(rng, T)      for k in 1:modes]
 
    ts  = _angles(n)
    pts = map(ts) do t
        r = R + sum(a_k[k] * cos(k*t + φ_k[k]) for k in 1:modes)
        r = max(r, T(0.01) * R)   # prevent collapse to origin
        SVector{2,T}(r*cos(t), r*sin(t))
    end
    ClosedCurve(_transform(pts, SVector{2,T}(center), 0.0))
end
 
"""
    perlin_curve(n=256; octaves=4, persistence=0.5, radius=1.0, amplitude=0.3, seed, center, T)
 
Smooth random closed curve using Perlin-style gradient noise along the angular parameter.
 
- `octaves`: number of noise octaves (more = more detail).
- `persistence`: amplitude falloff per octave (0.5 = halve each octave).
- `amplitude`: maximum radial perturbation (fraction of radius).
 
```julia
perlin_curve()                            # smooth wobbly circle
perlin_curve(512; octaves=6, amplitude=0.5)   # more complex
```
"""
function perlin_curve(n::Int=256;
                       octaves     :: Int  = 4,
                       persistence :: Real = 0.5,
                       radius      :: Real = 1.0,
                       amplitude   :: Real = 0.3,
                       center             = SVector(0.0,0.0),
                       seed        :: Int  = 0,
                       T           :: Type{<:AbstractFloat} = Float64)
    rng = _seeded_rng(seed)
    R   = T(radius)
    amp = T(amplitude)
 
    tables = [_perlin_table(rng, 256) for _ in 1:octaves]
 
    ts  = _angles(n)
    pts = map(ts) do t
        x   = t / (2π)
        noise = zero(T)
        freq  = one(T)
        pers  = one(T)
        for oct in 1:octaves
            noise += pers * _perlin_1d_periodic(tables[oct], x * freq)
            freq  *= 2
            pers  *= T(persistence)
        end
        r = R + R * amp * noise
        r = max(r, T(0.01)*R)
        SVector{2,T}(r*cos(t), r*sin(t))
    end
    ClosedCurve(_transform(pts, SVector{2,T}(center), 0.0))
end

"""
    add_noise(c; amplitude=0.05, type=:gaussian, seed=0)
 
Add noise to every vertex of a curve.
 
- `:gaussian`  — i.i.d. Gaussian noise in each coordinate.
- `:normal`    — noise along the local normal direction only.
- `:tangential`— noise along the local tangent direction only.
- `:radial`    — noise radially from the centroid.
 
```julia
c_noisy = add_noise(c; amplitude=0.1)
c_noisy = add_noise(c; amplitude=0.05, type=:normal)
```
"""
function add_noise(c::AbstractDiscreteCurve{N,T};
                   amplitude :: Real   = T(0.05),
                   type      :: Symbol = :gaussian,
                   seed      :: Int    = 0) where {N,T}
    rng = _seeded_rng(seed)
    amp = T(amplitude)
    n   = nvertices(c)
    pts = Vector{SVector{N,T}}(undef, n)
 
    if type == :gaussian
        @inbounds for i in 1:n
            noise  = SVector{N,T}(ntuple(_ -> amp*(2rand(rng,T)-1), Val(N)))
            pts[i] = c.points[i] + noise
        end
 
    elseif type == :normal || type == :radial
        ctr = sum(c.points) / n
        @inbounds for i in 1:n
            p  = c.points[i]
            nv = p - ctr;  ℓ = norm(nv)
            nhat = ℓ > eps(T) ? nv/ℓ : SVector{N,T}(ntuple(d->d==1 ? one(T) : zero(T), Val(N)))
            pts[i] = p + amp*(2rand(rng,T)-1) * nhat
        end
 
    elseif type == :tangential
        @inbounds for i in 1:n
            ip   = isclosed(c) ? mod1(i-1,n) : max(i-1,1)
            iq   = isclosed(c) ? mod1(i+1,n) : min(i+1,n)
            tv   = c.points[iq] - c.points[ip]
            ℓ    = norm(tv)
            that = ℓ > eps(T) ? tv/ℓ : SVector{N,T}(ntuple(d->d==1 ? one(T) : zero(T), Val(N)))
            pts[i] = c.points[i] + amp*(2rand(rng,T)-1) * that
        end
 
    else
        throw(ArgumentError("noise type must be :gaussian, :normal, :radial, or :tangential"))
    end
    DiscreteCurve(pts; closed=isclosed(c))
end
"""
    rough_curve(c; amplitude=0.05, octaves=4, seed=0)
 
Add Perlin-noise radial perturbation to an existing curve.
Preserves the overall shape better than Gaussian noise.
"""
function rough_curve(c::AbstractDiscreteCurve{2,T};
                     amplitude :: Real = T(0.05),
                     octaves   :: Int  = 4,
                     seed      :: Int  = 0) where T
    rng = _seeded_rng(seed)
    n   = nvertices(c)
    amp = T(amplitude)
    tables = [_perlin_table(rng, 256) for _ in 1:octaves]
    ctr = sum(c.points) / n
 
    pts = Vector{SVector{2,T}}(undef, n)
    @inbounds for i in 1:n
        p  = c.points[i]
        nv = p - ctr;  ℓ = norm(nv)
        nhat = ℓ > eps(T) ? nv/ℓ : SVector{2,T}(one(T), zero(T))
        x    = T(i-1)/n
        noise = zero(T)
        freq  = one(T);  pers = one(T)
        for oct in 1:octaves
            noise += pers * _perlin_1d_periodic(tables[oct], x*freq)
            freq  *= 2;  pers *= T(0.5)
        end
        pts[i] = p + amp * noise * nhat
    end
    DiscreteCurve(pts; closed=isclosed(c))
end
 
"""
    noisy_polygon(n=32; radius=1.0, amplitude=0.2, type=:perlin, seed, T)
 
Generate a random closed polygon from scratch using noise.
Shorthand for `add_noise(regular_polygon(n; radius), ...)`.
"""
function noisy_polygon(n::Int=32;
                        radius    :: Real   = 1.0,
                        amplitude :: Real   = 0.2,
                        type      :: Symbol = :gaussian,
                        seed      :: Int    = 0,
                        T         :: Type{<:AbstractFloat} = Float64)
    c = inscribed_polygon(n; radius, irregularity=0.3, seed, T)
    add_noise(c; amplitude=T(amplitude)*radius, type, seed=seed+1)
end

function _perlin_table(rng, M::Int=256)
    perm = collect(0:M-1)

    for i in M:-1:2
        j = rand(rng, 1:i)
        perm[i], perm[j] = perm[j], perm[i]
    end
 
    grads = [2*(perm[i]%2) - 1 for i in 1:M]
    (perm=perm, grads=grads, M=M)
end
 
function _perlin_1d_periodic(table, x::T) where T<:AbstractFloat
    M   = table.M

    xw  = x * M
    x0  = Int(floor(xw)) % M
    x1  = (x0 + 1) % M
    t   = xw - floor(xw)

    fade = t * t * t * (t * (t*6 - 15) + 10)

    g0   = T(table.grads[x0+1]) * t
    g1   = T(table.grads[x1+1]) * (t - 1)
    return g0 + fade*(g1 - g0)
end
 
function _seeded_rng(seed::Int)
    seed == 0 ? Random.default_rng() : MersenneTwister(seed)
end
 
function _sample_closed(f, n::Int, ::Type{T}, t0::Real=0.0, t1::Real=2π) where T
    ts  = range(T(t0), T(t1); length=n+1)[1:n]
    ClosedCurve([SVector{2,T}(f(t)) for t in ts])
end
end