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

include("CurveData.jl")
using .CurveData

include("Macros.jl")
using .Macros

# Geometry
include("geometry/Lengths.jl"); using .Lengths
include("geometry/Curvature.jl"); using .Curvature

# ── Public API exports ─────────────────────────────────────────────────────────
export DiscreteCurve, ClosedCurve, OpenCurve, ParametricCurve
export isclosed, isopen, nvertices, nedges, vertex, edge_segment
export vertices, edges, edge_windows
export VertexData, EdgeData, attach_vertex_data, attach_edge_data, vertex_data, edge_data
export LazyVertexField, CurveDataBundle, bundle, channel_type, snapshot
export @curve,@where,@any_vertex,@all_vertices,@map_verts,@functional

export arc_length, edge_lengths, cumulative_length
export turning_angle, turning_angles, curvature, curvatures, curvature_vector
export TurningAngle, OsculatingCircle, TangentCurvature, SteinerCurvature, LengthWeightedCurvature, SignedCurvature2D, CustomCurvature
end # module DiscreteCurves
