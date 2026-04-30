module DiscreteCurves

using LinearAlgebra
using StaticArrays


include("types/AbstractTypes.jl")
using .AbstractTypes

include("types/CurveTypes.jl")
using .CurveTypes


include("CurveTopology.jl")
using .CurveTopology

include("Iterators.jl")
using .Iterators

include("CurveData.jl")
using .CurveData

include("Macros.jl")
using .Macros

include("geometry/Lengths.jl");
using .Lengths

include("geometry/GeometricProperties.jl")
using .GeometricProperties

include("geometry/Orientation.jl")
using .Orientation

include("geometry/CurveVectors.jl")
using .CurveVectors
include("geometry/Curvature.jl");
using .Curvature
include("geometry/Generators.jl");
using .Generators
include("geometry/Flows.jl");
using .Flows
include("graphics/Graphics.jl");
using .Graphics

# ── Public API exports ─────────────────────────────────────────────────────────
export DiscreteCurve, ClosedCurve, OpenCurve, ParametricCurve
export isclosed, isopen, nvertices, nedges, vertex, edge_segment
export vertices, edges, edge_windows
export VertexData, EdgeData, attach_vertex_data, attach_edge_data, vertex_data, edge_data
export LazyVertexField, CurveDataBundle, bundle, channel_type, snapshot
export @curve, @where, @any_vertex, @all_vertices, @map_verts, @functional

export arc_length, edge_lengths, cumulative_length
export tangent_vector, tangent_vectors, normal_vector, normal_vectors
export left_normal_vector, left_normal_vectors, right_normal_vector, right_normal_vectors
export inward_normal, outward_normal, inward_normals, outward_normals
export signed_area, centroid
export turning_angle, turning_angles, curvature, curvatures, curvature_vector, curvature_vectors
export TurningAngle, OsculatingCircle, TangentCurvature, SteinerCurvature, LengthWeightedCurvature, SignedCurvature2D, CustomCurvature
export Orientation, CCW, CW, CounterClockwise, Clockwise, orient, orient!, is_ccw, is_cw

export regular_polygon, inscribed_polygon, circle_curve, ellipse_curve, rectangle_curve,
	square_curve, triangle_curve, equilateral_triangle, right_triangle,
	parallelogram_curve, rhombus_curve, star_curve,
	cardioid, limacon, rose_curve, lemniscate, epicycloid, hypocycloid,
	astroid, deltoid, nephroid, trefoil_knot, figure_eight, heart_curve,
	kidney_curve, butterfly_curve, spiral_curve, logarithmic_spiral,
	fourier_curve, perlin_curve, random_convex_polygon, random_concave_polygon, add_noise, rough_curve, noisy_polygon

export CurvatureFlow, PointwiseFlow, BulkFlow, GradientFlow
export ScaledFlow, SumFlow, ConstrainedFlow, TangentialCorrection
export compute_velocities!, cfl_timestep, FlowResult, evolve, evolve_step!
export curveplot, curveplot!, plotcurve, vectorplot, vectorplot!, save_figure, @curveplot, @vectorplot
end
