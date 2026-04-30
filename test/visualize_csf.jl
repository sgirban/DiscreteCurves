using Pkg
Pkg.activate(".")
Pkg.instantiate()

using GLMakie
using DiscreteCurves
using Random

Random.seed!(1234)

c0 = random_convex_polygon(16; radius=1.0, irregularity=1.0, spikiness=0.5, T=Float64)

c0 = orient(c0, to = CW)
# Evolve under curve shortening with Steiner curvature
r = evolve(c0, CurvatureFlow(SteinerCurvature());
           nsteps=300,
           cfl_safety=1,
           save_every=10,
           track=[:L => arc_length])

# Collect curves: start + saved snapshots
curves = [c0; r.curves]

# Plot all curves together with distinct palette colors
fig = curveplot(curves;
    curve_color = :palette,
    curve_width = 2.5,
    vertex_visible = true,
    aspect = :equal,
    title = "Random convex polygon curve shortening flow",
    legend = false,)

# cf = curveplot(c0; curve_color=:red, curve_width=3.0, vertex_visible=true, label="Initial curve")
save_figure("visualize_csf2.png", fig)
print("Length variation during evolution:\n")
# for l in r[:L]
#     println("  L = $l")
# end
#println("Saved plot to visualize_csf.png")
