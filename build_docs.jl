#!/usr/bin/env julia
# Build documentation for DiscreteCurves
using Pkg

# Add the parent directory to the load path for development
push!(LOAD_PATH, dirname(@__FILE__))

# Now include the make script
include("docs/make.jl")

