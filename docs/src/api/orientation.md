# Orientation

Functions and types for querying and setting the orientation (CCW vs CW) of 2D
closed curves.

---

## Orientation

```julia
@enum Orientation CCW CW
const CounterClockwise = CCW
const Clockwise        = CW
```

An enum with two values indicating whether a closed 2D curve is traversed
**counter-clockwise** (CCW, positive mathematical convention) or **clockwise** (CW).

| Value | Alias | Signed area |
|-------|-------|-------------|
| `CCW` | `CounterClockwise` | `signed_area(c) > 0` |
| `CW`  | `Clockwise`        | `signed_area(c) < 0` |

### Examples

```julia
# Use Orientation values directly
c = orient(circle_curve(128); to=CCW)
c = orient(circle_curve(128); to=Clockwise)   # alias

# Pattern-match in a function
function describe(c)
    is_ccw(c) ? "CCW (positive area)" : "CW (negative area)"
end
```

---

## is\_ccw / is\_cw

```julia
is_ccw(c) → Bool
is_cw(c)  → Bool
```

Query the orientation of a 2D curve by sign of `signed_area(c)`.

Works for both open and closed curves. For open curves, the area is computed as
if the endpoints were connected by a straight line.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{2}` | A 2D curve. |

### Returns

`Bool`.

### Examples

```julia
# Default generators produce CCW curves
c = circle_curve(128)
is_ccw(c)   # → true
is_cw(c)    # → false

# CW circle
c_cw = orient(c; to=CW)
is_ccw(c_cw)   # → false
is_cw(c_cw)    # → true

# Signed curvature depends on orientation
κs_ccw = curvatures(circle_curve(128), SignedCurvature2D())   # all > 0
κs_cw  = curvatures(orient(circle_curve(128); to=CW), SignedCurvature2D())  # all < 0

# Guard in functions that require a specific orientation
function compute_winding(c)
    is_ccw(c) || error("Expected CCW curve")
    # ...
end
```

### See Also

[`orient`](@ref), [`orient!`](@ref), [`signed_area`](@ref)

---

## orient

```julia
orient(c; to::Orientation = CCW) → DiscreteCurve
```

Return a **new** curve with the specified orientation. If the curve already has
the requested orientation, a copy is returned. Otherwise, the vertex order is
reversed.

**Non-mutating** (creates a new `DiscreteCurve`). For in-place modification,
use [`orient!`](@ref).

Requires a **closed** 2D curve; throws `ArgumentError` for open curves.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `AbstractDiscreteCurve{2}` | A closed 2D curve. |
| `to` | `Orientation` | Target orientation. Default: `CCW`. |

### Returns

`DiscreteCurve{2,T}` — new curve with the requested orientation.

### Examples

```julia
# Ensure CCW before computing signed area
c     = random_convex_polygon(32)
c_ccw = orient(c; to=CCW)
signed_area(c_ccw)   # → positive

c_cw  = orient(c; to=CW)
signed_area(c_cw)    # → negative

# Normalise orientation before a flow
c = fourier_curve(256; seed=1)
c = orient(c; to=CCW)   # ensure CCW
r = evolve(c, CurvatureFlow())

# Aliases
c1 = orient(c; to=CounterClockwise)
c2 = orient(c; to=Clockwise)

# Identity case (already correct orientation)
c_ccw = circle_curve(128)     # already CCW
c_same = orient(c_ccw; to=CCW)   # returns a copy, no reversal
```

### Notes

Reversing the vertex order also reverses the signed curvature and the winding
number. The unsigned curvature (`TurningAngle`, `SteinerCurvature`, `OsculatingCircle`)
is **unchanged** by orientation reversal.

### See Also

[`orient!`](@ref), [`is_ccw`](@ref), [`is_cw`](@ref)

---

## orient!

```julia
orient!(c::DiscreteCurve{2}; to::Orientation = CCW) → DiscreteCurve
```

**In-place** version of [`orient`](@ref). Reverses the internal `points` vector
if the current orientation does not match `to`. Returns the modified curve.

Requires a **closed** 2D `DiscreteCurve` (not just any `AbstractDiscreteCurve`).
Throws `ArgumentError` for open curves.

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `c` | `DiscreteCurve{2}` | A closed concrete `DiscreteCurve`. |
| `to` | `Orientation` | Target orientation. Default: `CCW`. |

### Returns

`c` (the same object, mutated if needed).

### Examples

```julia
c = random_convex_polygon(64; seed=5)
orient!(c; to=CCW)       # normalise to CCW in-place
is_ccw(c)                # → true

orient!(c; to=CW)        # flip in-place
is_cw(c)                 # → true

# Batch normalise a collection
curves = [fourier_curve(128; seed=k) for k in 1:10]
orient!.(curves; to=CCW)   # broadcast in-place on each curve
all(is_ccw, curves)        # → true
```

### Notes

Prefer `orient!` over `orient` when working with large curves in a loop, to
avoid repeated allocation. The modification is O(n) (pointer reversal).

### See Also

[`orient`](@ref), [`is_ccw`](@ref)