using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DiscreteCurves
using LinearAlgebra

c = DiscreteCurve([(cos(θ), sin(θ), θ) for θ in LinRange(0, 2π, 100)])
save_curve("testcurve.ply", c)
print(c)