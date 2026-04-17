using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DiscreteCurves
using LinearAlgebra
using GLMakie
using StaticArrays
# c = DiscreteCurve([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)], closed = true)
#c = random_convex_polygon(5, irregularity = 0.7)
# c = regular_polygon(7)
c = fourier_curve(10,roughness = 1.0,radius = 5.0)
println("Curve points: ", c.points)
kws= curvature_vectors(c)
kws = -1 .* curvature_vectors(c)
curveplot(c, tooltip = [:curvature, :turning_angle], 
         skin=:studio) +
         vectorplot(kws; anchors=c, arrow_color=:orange, normalize=false,lengthscale=3.0)