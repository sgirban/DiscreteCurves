# Curve Generators

Ready-made constructors for standard geometric shapes, mathematical curves, and
random curves. All generators return closed `DiscreteCurve{2,T}` objects unless
noted otherwise.

---

## Basic Polygons

### regular\_polygon

```julia
regular_polygon(n; center = SVector(0,0), radius = 1.0, angle = 0.0, T = Float64)
```

Regular `n`-gon inscribed in a circle of the given `radius`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n` | `Int` | required | Number of vertices (sides). Must be ≥ 3. |
| `center` | `SVector{2}` | `(0,0)` | Centre of the polygon. |
| `radius` | `Real` | `1.0` | Circumradius (vertex-to-centre distance). |
| `angle` | `Real` | `0.0` | Initial rotation in radians. |
| `T` | `Type` | `Float64` | Scalar precision. |

```julia
regular_polygon(3)                                   # equilateral triangle
regular_polygon(4)                                   # square (circumradius 1)
regular_polygon(6; radius=2.0)                       # hexagon, circumradius 2
regular_polygon(5; center=SVector(1.0,1.0), angle=π/10)   # rotated pentagon
```

---

### inscribed\_polygon

```julia
inscribed_polygon(n; center, radius, irregularity, seed, T)
```

`n` points on a circle with controllable angular irregularity. Use to generate
"wobbly" polygons that are still roughly circular.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n` | required | Number of vertices. |
| `irregularity` | `0.0` | Angular spacing noise ∈ [0,1]. `0` = regular polygon, `1` = fully random. |
| `seed` | `0` | RNG seed for reproducibility. `0` = system RNG. |

```julia
inscribed_polygon(8)                            # regular octagon
inscribed_polygon(8; irregularity=0.5)          # slightly irregular
inscribed_polygon(20; irregularity=1.0, seed=3) # fully random spacing
```

---

### circle\_curve

```julia
circle_curve(n = 128; center, radius, T)
```

Uniformly sampled circle (alias for `regular_polygon`).

```julia
c = circle_curve()            # 128 vertices, unit circle
c = circle_curve(512)         # 512 vertices, finer
c = circle_curve(256; radius=3.0, center=SVector(1.0, -2.0))
arc_length(c)   # → ≈ 2π·radius
```

---

### ellipse\_curve

```julia
ellipse_curve(n = 128; a = 1.0, b = 0.5, center, angle, T)
```

Ellipse with semi-axes `a` (horizontal) and `b` (vertical), uniformly sampled
in parameter (not arc-length).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `a` | `1.0` | Semi-axis along x. |
| `b` | `0.5` | Semi-axis along y. |
| `angle` | `0.0` | Rotation of the major axis (radians). |

```julia
c = ellipse_curve()                              # a=1, b=0.5
c = ellipse_curve(256; a=3.0, b=1.0)            # elongated
c = ellipse_curve(256; a=2.0, b=1.0, angle=π/4) # rotated
```

---

### rectangle\_curve

```julia
rectangle_curve(; width = 2.0, height = 1.0, center, angle, n = 64, T)
```

Axis-aligned rectangle with `n` vertices distributed proportionally to the
side lengths (more vertices on longer sides).

```julia
c = rectangle_curve()                           # 2×1 rectangle
c = rectangle_curve(; width=3.0, height=1.5, n=80)
c = rectangle_curve(; width=1.0, height=1.0)   # square (alias: square_curve)
```

---

### square\_curve

```julia
square_curve(; side = 2.0, kwargs...)
```

Shorthand for `rectangle_curve(; width=side, height=side, kwargs...)`.

```julia
c = square_curve()              # 2×2 square, 64 vertices
c = square_curve(; side=1.0, n=40)
```

---

### triangle\_curve

```julia
triangle_curve(p1, p2, p3; n = 48, T)
```

Triangle from three explicit vertex positions. Vertices are distributed
proportionally to the side lengths.

```julia
c = triangle_curve(SVector(0,0), SVector(1,0), SVector(0.5, 0.866))  # equilateral
c = triangle_curve((0,0), (3,0), (1,2); n=60)                        # scalene
```

---

### equilateral\_triangle

```julia
equilateral_triangle(; side = 2.0, center, angle, n = 48, T)
```

Equilateral triangle centred at the origin by default, with one vertex pointing up.

```julia
c = equilateral_triangle()
c = equilateral_triangle(; side=3.0, angle=π/6)
```

---

### right\_triangle

```julia
right_triangle(; a = 1.0, b = 1.0, n = 48, T)
```

Right triangle with the right angle at the origin, legs of length `a` (x) and
`b` (y).

```julia
c = right_triangle()              # 1-1-√2
c = right_triangle(; a=3.0, b=4.0)   # 3-4-5
```

---

### star\_curve

```julia
star_curve(n_points = 5; inner_r = 0.4, outer_r = 1.0, center, angle, T)
```

Regular star polygon with `n_points` points. Alternates between outer and inner
radii.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_points` | `5` | Number of points. Must be ≥ 3. |
| `inner_r` | `0.4` | Radius of inner vertices. |
| `outer_r` | `1.0` | Radius of outer vertices. |
| `angle` | `-π/2` | Orientation angle. |

```julia
star_curve(5)                              # classic 5-pointed star
star_curve(6; inner_r=0.57)               # Star of David shape
star_curve(8; inner_r=0.7, outer_r=1.2)  # 8-pointed star
```

---

### parallelogram\_curve / rhombus\_curve

```julia
parallelogram_curve(; a = 2.0, b = 1.0, skew = π/4, center, n = 64, T)
rhombus_curve(; side = 1.0, angle = π/3, n = 64, T)
```

Parallelogram with base `a`, side `b`, and interior angle `skew`. `rhombus_curve`
is a parallelogram with equal sides.

---

## Mathematical Curves

All mathematical curves are closed, centred at the origin by default.

### cardioid

```julia
cardioid(n = 256; a = 1.0, center, angle, T)
```

Heart-shaped curve: ``r(\theta) = a(1 - \cos\theta)``.

```julia
c = cardioid()
c = cardioid(512; a=2.0)
```

---

### limacon

```julia
limacon(n = 256; a = 0.5, b = 1.0, center, T)
```

Limaçon of Pascal: ``r(\theta) = a + b\cos\theta``.

- `b > a`: limaçon with inner loop
- `b == a`: cardioid
- `b < a`: convex limaçon (no inner loop)

```julia
limacon()                       # inner loop
limacon(; a=1.0, b=1.0)        # cardioid
limacon(; a=2.0, b=1.0)        # convex
```

---

### rose\_curve

```julia
rose_curve(k, n = 512; a = 1.0, center, T)
```

Rose curve: ``r(\theta) = a\cos(k\theta)``.

- Odd `k`: `k` petals.
- Even `k`: `2k` petals.

```julia
rose_curve(3)     # 3-petaled rose
rose_curve(4)     # 8-petaled rose
rose_curve(5)     # 5-petaled rose
rose_curve(2; a=1.5)   # 4-petaled, larger
```

---

### lemniscate

```julia
lemniscate(n = 256; a = 1.0, center, T)
```

Lemniscate of Bernoulli: a figure-eight shaped curve with
``r^2 = a^2\cos(2\theta)``.

```julia
c = lemniscate()
c = lemniscate(512; a=2.0)
```

---

### epicycloid / hypocycloid

```julia
epicycloid(k = 3, n = 512; R = 1.0, r, center, T)
hypocycloid(k = 3, n = 512; R = 1.0, r, center, T)
```

Roulette curves. `k` controls the number of cusps.

- `epicycloid`: small circle rolling **outside** circle of radius `R`. `k` cusps when `r = R/k`.
- `hypocycloid`: small circle rolling **inside** circle of radius `R`.
  - `k=3`: deltoid (3-cusped)
  - `k=4`: astroid (4-cusped)

```julia
epicycloid(3)     # 3-cusped epicycloid (nephroid-like)
epicycloid(5)     # 5-cusped
hypocycloid(3)    # deltoid
hypocycloid(4)    # astroid
```

---

### astroid / deltoid / nephroid

```julia
astroid(n = 256; a = 1.0, center, T)
deltoid(n = 256; R = 1.0, T)
nephroid(n = 256; a = 1.0, center, T)
```

Named hypocycloid/epicycloid shortcuts.

- **Astroid**: ``x^{2/3} + y^{2/3} = a^{2/3}`` (4-cusped)
- **Deltoid**: 3-cusped hypocycloid
- **Nephroid**: 2-cusped epicycloid (kidney-shaped)

```julia
c = astroid(512)
c = deltoid(256)
c = nephroid(256; a=1.5)
```

---

### trefoil\_knot

```julia
trefoil_knot(n = 512; a = 2.0, center, T)
```

Trefoil knot projected onto 2D:
``x(t) = \sin t + a\sin 2t``, ``y(t) = \cos t - a\cos 2t``.

```julia
c = trefoil_knot()
c = trefoil_knot(512; a=1.5)
```

---

### figure\_eight

```julia
figure_eight(n = 256; a = 1.0, center, T)
```

Figure-eight / lemniscate of Gerono:
``x = a\sin t``, ``y = a\sin t \cos t``.

---

### heart\_curve

```julia
heart_curve(n = 256; scale = 1.0, center, T)
```

Classic algebraic heart curve (parametric form).

```julia
c = heart_curve()
c = heart_curve(512; scale=2.0)
```

---

### butterfly\_curve

```julia
butterfly_curve(n = 1024; a = 1.0, T)
```

Temple H. Fay's butterfly curve:
``r = e^{\sin\theta} - 2\cos(4\theta) + \sin^5\!\left(\frac{2\theta-\pi}{24}\right)``.

Requires `n ≥ 1024` for a clean result (parameter domain is `[0, 4π]`).

---

### kidney\_curve / limacon

```julia
kidney_curve(n = 256; a = 1.0, center, T)
```

Kidney-shaped curve: ``r(\theta) = a(1 + 0.6\cos\theta)``.

---

### spiral\_curve (open)

```julia
spiral_curve(turns = 3, n = 512; r_start = 0.1, r_end = 1.0, center, T)
```

Archimedean spiral. Returns an **open** `DiscreteCurve`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `turns` | `3` | Number of full turns. |
| `r_start` | `0.1` | Inner radius. |
| `r_end` | `1.0` | Outer radius. |

```julia
c = spiral_curve()
c = spiral_curve(5; r_start=0.05, r_end=2.0)
isopen(c)   # → true
```

---

### logarithmic\_spiral (open)

```julia
logarithmic_spiral(turns = 3, n = 512; a = 0.1, b = 0.2, center, T)
```

Equiangular spiral: ``r = ae^{b\theta}``. Returns an **open** `DiscreteCurve`.

---

## Random Curves

### random\_convex\_polygon

```julia
random_convex_polygon(n = 16; radius = 1.0, irregularity = 0.5,
                      spikiness = 0.2, center, seed = 0, T)
```

Random convex polygon via Valtr's algorithm.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `irregularity` | `0.5` | Angular spacing variation ∈ [0,1]. `0` = regular, `1` = fully random. |
| `spikiness` | `0.2` | Radial variation ∈ [0,1]. `0` = circle, `1` = very spiky. |
| `seed` | `0` | RNG seed. `0` = system RNG. |

```julia
random_convex_polygon(20)
random_convex_polygon(16; irregularity=0.8, spikiness=0.4, seed=42)
random_convex_polygon(8; irregularity=0.0, spikiness=0.0)   # regular octagon
```

---

### random\_concave\_polygon

```julia
random_concave_polygon(n = 16; radius = 1.0, concavity = 0.5, center, seed = 0, T)
```

Random non-convex polygon with `concavity` fraction of vertices pushed inward.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `concavity` | `0.5` | Fraction of concave vertices ∈ [0,1]. |

```julia
random_concave_polygon(20)
random_concave_polygon(32; concavity=0.7, seed=7)
```

---

### fourier\_curve

```julia
fourier_curve(n = 256; modes = 8, radius = 1.0, roughness = 0.5, center, seed = 0, T)
```

Smooth random closed curve generated by summing random Fourier modes with
amplitude decaying as ``1/k^\alpha``.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `modes` | `8` | Number of Fourier modes. More = richer shape. |
| `roughness` | `0.5` | High-frequency weight ∈ [0,1]. `0` = very smooth, `1` = complex. |
| `seed` | `0` | RNG seed for reproducibility. |

```julia
fourier_curve()                               # smooth random blob
fourier_curve(512; seed=42)                   # reproducible
fourier_curve(256; modes=4, roughness=0.1)   # very smooth, nearly circular
fourier_curve(256; modes=20, roughness=0.9)  # highly complex boundary
fourier_curve(512; seed=1)                   # one specific shape
fourier_curve(512; seed=2)                   # a different shape
```

This is the most commonly used generator for testing and examples.

---

### perlin\_curve

```julia
perlin_curve(n = 256; octaves = 4, persistence = 0.5, radius = 1.0,
             amplitude = 0.3, center, seed = 0, T)
```

Smooth random closed curve using Perlin-style gradient noise.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `octaves` | `4` | Number of noise octaves. More = more detail. |
| `persistence` | `0.5` | Amplitude falloff per octave (0.5 = halved each time). |
| `amplitude` | `0.3` | Maximum radial perturbation as a fraction of `radius`. |

```julia
perlin_curve()
perlin_curve(512; octaves=6, amplitude=0.5, seed=3)
```

---

## Noise and Perturbation

### add\_noise

```julia
add_noise(c; amplitude = 0.05, type = :gaussian, seed = 0)
```

Add noise to every vertex of a curve. Returns a new `DiscreteCurve`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `amplitude` | `0.05` | Noise magnitude. |
| `type` | `:gaussian` | Noise type: `:gaussian`, `:normal`, `:tangential`, or `:radial`. |
| `seed` | `0` | RNG seed. |

**Noise types:**

| Type | Description |
|------|-------------|
| `:gaussian` | Independent Gaussian noise in each coordinate. |
| `:normal` (= `:radial`) | Noise directed radially from the centroid. |
| `:tangential` | Noise along the local tangent direction (perturbs parametrisation). |

```julia
c = circle_curve(128)
c_noisy = add_noise(c; amplitude=0.05)                     # gaussian
c_bump  = add_noise(c; amplitude=0.1, type=:normal)        # radial bumps
c_shift = add_noise(c; amplitude=0.05, type=:tangential)   # tangential only
```

---

### rough\_curve

```julia
rough_curve(c; amplitude = 0.05, octaves = 4, seed = 0)
```

Add Perlin-noise radial perturbation to an existing curve. Preserves the overall
shape better than Gaussian noise (correlated, fractal-like texture).

```julia
c_smooth = circle_curve(256)
c_rough  = rough_curve(c_smooth; amplitude=0.1, octaves=6)
c_rough2 = rough_curve(fourier_curve(256; seed=1); amplitude=0.05)
```

---

### noisy\_polygon

```julia
noisy_polygon(n = 32; radius = 1.0, amplitude = 0.2, type = :gaussian, seed = 0, T)
```

Generate a random closed polygon from scratch. Shorthand for
`add_noise(inscribed_polygon(n; radius), ...)`.

```julia
c = noisy_polygon()
c = noisy_polygon(64; amplitude=0.3, type=:normal, seed=5)
```