module DiscreteCurves

using LinearAlgebra
using StaticArrays

# Types
include("types/AbstractTypes.jl")
using .AbstractTypes

include("types/CurveTypes.jl")
using .CurveTypes

# Topology
include("CurveTopology.jl")
using .CurveTopology
# Iterators
include("Iterators.jl")
using .Iterators

# ── Public API exports ─────────────────────────────────────────────────────────
export DiscreteCurve, ClosedCurve, OpenCurve, ParametricCurve
export isclosed, isopen, nvertices, nedges, vertex, edge_segment
export vertices, edges, edge_windows
end # module DiscreteCurves
