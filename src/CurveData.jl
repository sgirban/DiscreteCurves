module CurveData
    using StaticArrays
    using ..AbstractTypes
    using ..CurveTypes
    using ..CurveTopology
export VertexData, EdgeData, attach_vertex_data, attach_edge_data, vertex_data, edge_data
export LazyVertexField, CurveDataBundle, bundle, channel_type, snapshot
    
# ─────────────────────────────────────────────────────────────────────────────
#               VertexData
# ─────────────────────────────────────────────────────────────────────────────
 """
    VertexData{D, C}
 
A typed array of values of type `D`, one per vertex of curve `C`.
 
Use for pre-computed fields (curvature, temperature, pressure, signal, …).
Topology-aware: indexing wraps for closed curves.
"""
struct VertexData{D, C <: AbstractDiscreteCurve}
    curve :: C
    data  :: Vector{D}
 
    function VertexData{D}(c::C, d::Vector{D}) where {D, C <: AbstractDiscreteCurve}
        length(d) == nvertices(c) || throw(DimensionMismatch(
            "VertexData: $(length(d)) values but nvertices=$(nvertices(c))"))
        new{D,C}(c, d)
    end
end

"""
    attach_vertex_data(c, data)
    attach_vertex_data(c, f)       # f: SVector → value
 
Attach vertex data to curve `c`.  `data` can be a `Vector` or any callable.
 
```julia
vd = attach_vertex_data(c, rand(nvertices(c)))
vd = attach_vertex_data(c, p -> norm(p))
vd = attach_vertex_data(c, p -> p[1]^2 + p[2]^2)   # radial distance squared
```
"""
function attach_vertex_data(c::AbstractDiscreteCurve, data::Vector{D}) where D
    VertexData{D}(c, copy(data))
end
 
function attach_vertex_data(c::AbstractDiscreteCurve, f)
    T_first = typeof(f(vertex(c, 1)))
    data    = Vector{T_first}(undef, nvertices(c))
    @inbounds for i in 1:nvertices(c)
        data[i] = f(vertex(c, i))
    end
    VertexData{T_first}(c, data)
end

@inline function Base.getindex(vd::VertexData, i::Int)
    @inbounds vd.data[mod1(i, length(vd.data))]
end
@inline function Base.setindex!(vd::VertexData, val, i::Int)
    @inbounds vd.data[mod1(i, length(vd.data))] = val
end
@inline vertex_data(vd::VertexData, i::Int) = vd[i]

Base.length(vd::VertexData)          = length(vd.data)
Base.size(vd::VertexData)            = (length(vd.data),)
Base.eltype(::VertexData{D}) where D = D
Base.firstindex(::VertexData)        = 1
Base.lastindex(vd::VertexData)       = length(vd.data)
Base.iterate(vd::VertexData, i=1)   = i > length(vd) ? nothing : (vd.data[i], i+1)
Base.collect(vd::VertexData)         = copy(vd.data)
Base.broadcastable(vd::VertexData)   = vd.data
 
Base.show(io::IO, vd::VertexData{D}) where D =
    print(io, "VertexData{$D}  $(length(vd)) vertices")
# ── EdgeData ──────────────────────────────────────────────────────────────────
 
"""
    EdgeData{D, C}
 
A typed array of values of type `D`, one per edge of curve `C`.
 
```julia
ed = attach_edge_data(c, rand(nedges(c)))
ed = attach_edge_data(c, (p,q) -> norm(q-p))   # from endpoint pair
```
"""
struct EdgeData{D, C <: AbstractDiscreteCurve}
    curve :: C
    data  :: Vector{D}
 
    function EdgeData{D}(c::C, d::Vector{D}) where {D, C <: AbstractDiscreteCurve}
        length(d) == nedges(c) || throw(DimensionMismatch(
            "EdgeData: $(length(d)) values but nedges=$(nedges(c))"))
        new{D,C}(c, d)
    end
end
 
function attach_edge_data(c::AbstractDiscreteCurve, data::Vector{D}) where D
    EdgeData{D}(c, copy(data))
end
 
function attach_edge_data(c::AbstractDiscreteCurve, f)
    T_first = typeof(f(vertex(c, 1), vertex(c, 2)))
    data    = Vector{T_first}(undef, nedges(c))
    @inbounds for i in 1:nedges(c)
        data[i] = f(vertex(c, i), vertex(c, i+1))
    end
    EdgeData{T_first}(c, data)
end
 
@inline function Base.getindex(ed::EdgeData, i::Int)
    @inbounds ed.data[mod1(i, length(ed.data))]
end
 
@inline function Base.setindex!(ed::EdgeData, val, i::Int)
    @inbounds ed.data[mod1(i, length(ed.data))] = val
end
 
@inline edge_data(ed::EdgeData, i::Int) = ed[i]
 
Base.length(ed::EdgeData)          = length(ed.data)
Base.size(ed::EdgeData)            = (length(ed.data),)
Base.eltype(::EdgeData{D}) where D = D
Base.firstindex(::EdgeData)        = 1
Base.lastindex(ed::EdgeData)       = length(ed.data)
Base.iterate(ed::EdgeData, i=1)   = i > length(ed) ? nothing : (ed.data[i], i+1)
Base.collect(ed::EdgeData)         = copy(ed.data)
Base.broadcastable(ed::EdgeData)   = ed.data
 
Base.show(io::IO, ed::EdgeData{D}) where D =
    print(io, "EdgeData{$D}  $(length(ed)) edges")

# ─────────────────────────────────────────────────────────────────────────────
# LazyVertexField: zero-storage function-based data
# ─────────────────────────────────────────────────────────────────────────────
#
struct LazyVertexField{F,C<: AbstractDiscreteCurve}
    f ::F
    curve:: C 
end
"""
    lazy_vertex_field(c, f)
 
Create a zero-storage field that evaluates `f(vertex(c,i))` on demand.
 
```julia
lf = lazy_vertex_field(c, p -> norm(p))
lf[5]         # evaluates norm(vertex(c,5))
collect(lf)   # materialises into a VertexData
```
"""
lazy_vertex_field(c::AbstractDiscreteCurve, f) = LazyVertexField(f,c)

@inline function Base.getindex(lf::LazyVertexField,i::Int)
    lf.f(vertex(lf.curve,mod1(i,nvertices(lf.curve))))
end

Base.length(lf::LazyVertexField) = nvertices(lf.curve)
Base.iterate(lf::LazyVertexField, i=1) =
    i > nvertices(lf.curve) ? nothing : (lf[i], i+1)
Base.collect(lf::LazyVertexField) = attach_vertex_data(lf.curve, lf.f)
Base.show(io::IO, lf::LazyVertexField) =
    print(io, "LazyVertexField  $(nvertices(lf.curve)) vertices  [computed on access]")
# ─────────────────────────────────────────────────────────────────────────────
#  CurveDataBundle
# ─────────────────────────────────────────────────────────────────────────────
 
"""
    CurveDataBundle(c)
    bundle(c)
 
A named container of multiple data channels on a single curve.
Channels are referenced by `Symbol` keys.
 
### Attach channels
```julia
c = ParametricCurve(t -> SVector(cos(t), sin(t)), 0, 2π, 128; closed=true)
b = bundle(c)
 
b[:curvature]   = curvatures(c)          # Vector{Float64} → auto-wrapped as VertexData
b[:temperature] = rand(nvertices(c))     # same
b[:edge_len]    = edge_lengths(c)        # nedges → auto-wrapped as EdgeData
b[:edge_len,toEdges=true]    = edge_lengths(c)        # nedges → auto-wrapped as EdgeData
b[:norm]        = p -> norm(p)           # function → LazyVertexField
```
 
### Access channels
```julia
b[:curvature]           # the whole VertexData channel
b[:curvature, 5]        # value at vertex 5
vertex_data(b, 5, :curvature)   # explicit accessor (identical)
b.curvature[5]          # dot-notation shortcut
 
snapshot(b, 5)          # NamedTuple of all vertex channels at vertex 5
```
 
### Iterate
```julia
for κ in b[:curvature]          # iterates like a Vector
    println(κ)
end
 
for (p, κ) in zip(vertices(c), b[:curvature])   # co-iterate with geometry
    process(p, κ)
end
```
 
### Performance tip
Extract the channel *once* before a tight loop:
```julia
temp = b[:temperature]           # one Dict lookup
for i in 1:nvertices(c)
    do_work(vertex(c,i), temp[i])   # plain array read — no Dict overhead
end
```
"""
mutable struct CurveDataBundle{C <: AbstractDiscreteCurve}
    curve    :: C
    channels :: Dict{Symbol, Any}
 
    CurveDataBundle(c::C) where C <: AbstractDiscreteCurve =
        new{C}(c, Dict{Symbol, Any}())
end
 
bundle(c::AbstractDiscreteCurve) = CurveDataBundle(c)
 

function Base.setindex!(b::CurveDataBundle, data::Vector{D}, key::Symbol;toEdges=false) where D
    n, nv, ne = length(data), nvertices(b.curve), nedges(b.curve)
    if n == nv && !toEdges
        b.channels[key] = VertexData{D}(b.curve, copy(data))
    elseif n == ne
        b.channels[key] = EdgeData{D}(b.curve, copy(data))
    else
        throw(DimensionMismatch(
            "Vector length $n ≠ nvertices=$nv nor nedges=$ne"))
    end
end
 
Base.setindex!(b::CurveDataBundle, vd::VertexData, key::Symbol) =
    (b.channels[key] = vd; b)
 
Base.setindex!(b::CurveDataBundle, ed::EdgeData, key::Symbol) =
    (b.channels[key] = ed; b)

Base.setindex!(b::CurveDataBundle, f::Function, key::Symbol) =
    (b.channels[key] = LazyVertexField(f, b.curve); b)
 
function Base.getindex(b::CurveDataBundle, key::Symbol)
    haskey(b.channels, key) || throw(KeyError(
        "Channel :$key not found. Available: $(sort!(collect(keys(b.channels))))"))
    b.channels[key]
end

@inline Base.getindex(b::CurveDataBundle, key::Symbol, i::Int) = b[key][i]
@inline vertex_data(b::CurveDataBundle, i::Int, ch::Symbol) = b[ch, i]
@inline edge_data(b::CurveDataBundle,   i::Int, ch::Symbol) = b[ch, i]
 

function Base.getproperty(b::CurveDataBundle, key::Symbol)
    key ∈ (:curve, :channels) && return getfield(b, key)
    b.channels[key]
end
 
Base.propertynames(b::CurveDataBundle) =
    (fieldnames(CurveDataBundle)..., keys(b.channels)...)
Base.keys(b::CurveDataBundle)             = keys(b.channels)
Base.haskey(b::CurveDataBundle, k::Symbol) = haskey(b.channels, k)
Base.length(b::CurveDataBundle)           = length(b.channels)
 
"""
    channel_type(b, key)  →  :vertex | :edge | :lazy | :unknown
"""
function channel_type(b::CurveDataBundle, key::Symbol)
    ch = b[key]
    ch isa VertexData      && return :vertex
    ch isa EdgeData        && return :edge
    ch isa LazyVertexField && return :lazy
    return :unknown
end
 
Base.iterate(b::CurveDataBundle, state...) = iterate(b.channels, state...)
 
"""
    snapshot(b, i)  →  NamedTuple
 
Collect all vertex/lazy channels into a NamedTuple at vertex `i`.
EdgeData channels are skipped.
 
```julia
s = snapshot(b, 5)
# → (curvature = 0.23, temperature = 1.14, norm = 1.0)
s.curvature   # NamedTuple field access
```
"""
function snapshot(b::CurveDataBundle, i::Int)
    ks = Symbol[]
    vs = Any[]
    for (k, ch) in b.channels
        (ch isa VertexData || ch isa LazyVertexField) || continue
        push!(ks, k)
        push!(vs, ch[i])
    end
    NamedTuple{Tuple(ks)}(Tuple(vs))
end
 
 
function Base.show(io::IO, b::CurveDataBundle)
    c = b.curve
    println(io, "CurveDataBundle  $(nvertices(c))v / $(nedges(c))e")
    for (k, ch) in sort!(collect(b.channels); by=x->string(x[1]))
        tag  = ch isa VertexData ? "vertex" :
               ch isa EdgeData   ? "edge"   :
               ch isa LazyVertexField ? "lazy" : "?"
        D    = ch isa LazyVertexField ? "on demand" : string(eltype(ch.data))
        n    = ch isa LazyVertexField ? nvertices(c) : length(ch)
        println(io, "  :$(rpad(k,14)) [$tag, $D, $n entries]")
    end
end
end