using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DiscreteCurves
using LinearAlgebra


c = DiscreteCurve([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)], closed=true)
for (p, q) in edges(c)
    println(q-p)
end

print("Arc length: ", arc_length(c), "\n")
print("Edge lengths: ", edge_lengths(c), "\n")
print("Cumulative length: ", cumulative_length(c), "\n")

print("Turning angles: ", turning_angles(c, signed=true), "\n")
print("Errors: ", turning_angles(c, signed=true) .- [π/2, π/2, π/2, π/2], "\n")


# d = DiscreteCurve(@select p in c: p[1] > 0, closed=false)

# for p in vertices(d)
#     println(p)
# end