#!/usr/bin/env julia

using Pkg


Pkg.instantiate()

include("./runtests.jl")
