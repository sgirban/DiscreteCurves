using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DiscreteCurves
using LinearAlgebra
using GLMakie
using StaticArrays
# c = DiscreteCurve([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)], closed = true)
#c = random_convex_polygon(5, irregularity = 0.7)
c = regular_polygon(4)
println("Curve points: ", c.points)
curveplot(c, labels = :index, title = "Square Curve", controls = true)
