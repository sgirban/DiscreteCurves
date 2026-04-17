# Mathematical Background

## Discrete differential geometry of curves

A **discrete curve** is a piecewise-linear path: an ordered sequence of vertices
``\gamma_1, \ldots, \gamma_n \in \mathbb{R}^N`` connected by straight edges.

For a **closed** curve, vertex ``n+1`` is identified with vertex ``1``.
For an **open** curve, ``\gamma_1`` and ``\gamma_n`` are free endpoints.

## Tangent vectors

The **edge tangent** at edge ``i`` is the chord:

```math
T_i^e = \gamma_{i+1} - \gamma_i
```

The **vertex tangent** at vertex ``i`` is the centred difference:

```math
T_i^v = \frac{\gamma_{i+1} - \gamma_{i-1}}{2}
```

(one-sided differences at open endpoints). This gives **O(h²)** accuracy vs O(h)
for a simple forward difference.

## Normal vectors (2D)

The outward unit normal at vertex ``i`` is the unit vertex tangent rotated 90° CCW:

```math
N_i = \frac{1}{|T_i^v|}\begin{pmatrix} -(T_i^v)_y \\ (T_i^v)_x \end{pmatrix}
```

This is equivalent to ``R(90°) \hat{T}_i`` where ``R`` is the rotation matrix.

## Discrete curvature (Crane & Wardetzky 2017)

From the DDG "game": write several equivalent smooth definitions, apply each
to the discrete setting, and analyse trade-offs.

| Model | Formula | Property |
|-------|---------|---------|
| ``\kappa^A`` (TurningAngle) | ``\theta_i / \bar{\ell}_i`` | ``\sum \kappa_i^A \bar{\ell}_i = 2\pi`` (Gauss-Bonnet) |
| ``\kappa^B`` (Steiner) | ``2\sin(\theta_i/2) / \bar{\ell}_i`` | ``-\partial_{\gamma_i} L = \kappa_i^B N_i`` |
| ``\kappa^C`` (Tangent) | ``2\tan(\theta_i/2) / \bar{\ell}_i`` | Kirchhoff analogy |
| ``\kappa^D`` (Osculating) | ``2\sin(\theta_i)/w_i = 1/R_i`` | Circles are fixed points |

Here ``\theta_i`` is the **turning angle** at vertex ``i`` (angle between consecutive edges),
``\bar{\ell}_i = (\ell_{i-1} + \ell_i)/2`` is the **dual edge length**,
and ``w_i = |\gamma_{i+1} - \gamma_{i-1}|``.

**No free lunch:** no single discrete curvature preserves all properties of the smooth one simultaneously.

## Curvature vector (κᴮ · N)

The vector ``\kappa_i^B N_i = T_{i-1,i} - T_{i,i+1}`` is:

- the discrete gradient of arc length: ``-\partial_{\gamma_i} L``
- the velocity that most quickly reduces arc length
- the sum over all vertices = 0 for closed curves (preserves centroid)

## Curve-shortening flow

The continuous curve-shortening flow moves each point with velocity ``\kappa N``:

```math
\frac{\partial \gamma}{\partial t} = \kappa N
```

Discretised by explicit Euler:

```math
\gamma_i^{k+1} = \gamma_i^k + \tau \cdot v_i^k
```

where ``v_i = \kappa_i^B N_i`` (SteinerCurvature, preserves centroid), and the
**CFL stability condition** is ``\tau \leq h_{\min}^2 / 2`` where ``h_{\min}`` is
the minimum edge length.

## Arc length and area

**Arc length** via the shoelace sum:

```math
L = \sum_{i=1}^n |\gamma_{i+1} - \gamma_i|
```

**Signed area** for 2D closed curves (Green's theorem / shoelace formula):

```math
A = \frac{1}{2} \sum_{i=1}^n (\gamma_i \times \gamma_{i+1})
  = \frac{1}{2} \sum_{i=1}^n (x_i y_{i+1} - x_{i+1} y_i)
```

**Isoperimetric ratio**: ``r = 4\pi A / L^2 \in (0,1]``, equals 1 iff the curve is a circle.

## Turning angle and Gauss-Bonnet

For a smooth closed curve winding once counterclockwise:

```math
\oint \kappa \, ds = 2\pi
```

Discretely (using ``\kappa^A``):

```math
\sum_{i=1}^n \theta_i = 2\pi
```

This is an **exact identity** for ``\kappa^A``, which is why ``\kappa^A`` flow
preserves total curvature.

## References

- Crane, K. & Wardetzky, M. (2017). *A Glimpse into Discrete Differential Geometry*. Notices of the AMS.
- Bobenko, A. & Suris, Y. (2008). *Discrete Differential Geometry*. AMS.
- Keenan Crane's DDG course notes: [ddg.cs.cmu.edu](https://ddg.cs.cmu.edu)
