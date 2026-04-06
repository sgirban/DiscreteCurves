using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DiscreteCurves
using LinearAlgebra


c = DiscreteCurve([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)], closed=true)
for (p, q) in edges(c)
    println(q-p)
end