module CurveIO

using StaticArrays
using ..AbstractTypes, ..CurveTypes, ..CurveTopology

export save_curve, load_curve, save_curves, load_curves

const DCB_MAGIC_SINGLE = b"DCrv\x01\x00"
const DCB_MAGIC_MULTI  = b"DCrvm\x01"

const T_CODE = Dict{Type, UInt8}(Float32 => 0x00, Float64 => 0x01)
const CODE_T = Dict{UInt8, Type}(0x00 => Float32, 0x01 => Float64)


function _detect_format(path::AbstractString)::Symbol
    ext = lowercase(splitext(path)[2])
    ext == ".dcb"  && return :binary
    ext == ".csv"  && return :csv
    ext == ".tsv"  && return :tsv
    ext == ".xyz"  && return :xyz
    ext == ".ply"  && return :ply
    ext == ".txt"  && return :xyz
    error("Cannot detect format from extension '$ext'. " *
          "Pass format= explicitly: :binary, :csv, :xyz, or :ply")
end


"""
    save_curve(path, c; format=:auto, T=nothing)

Write a single `DiscreteCurve` to `path`.

**Supported formats** (auto-detected from extension when `format=:auto`):

| Extension | Format   | Description                              |
|-----------|----------|------------------------------------------|
| `.dcb`    | `:binary`| Compact binary — fastest, lossless       |
| `.csv`    | `:csv`   | CSV with commented header                |
| `.tsv`    | `:tsv`   | TSV with commented header                |
| `.xyz`    | `:xyz`   | Whitespace-separated, no header          |
| `.ply`    | `:ply`   | ASCII PLY (GIS / point-cloud toolchains) |

`T` optionally downcasts floating-point (e.g. `T=Float32` to save space).

```julia
c = fourier_curve(512)

save_curve("curve.dcb", c)               # fast binary
save_curve("curve.csv", c)               # human-readable
save_curve("cloud.ply", c)               # PLY for MeshLab / CloudCompare
save_curve("f32.dcb",   c; T=Float32)   # half the binary size
```
"""
function save_curve(path::AbstractString,
                    c::AbstractDiscreteCurve;
                    format::Symbol = :auto,
                    T::Union{Nothing, Type{<:AbstractFloat}} = nothing)
    fmt = format === :auto ? _detect_format(path) : format
    c2  = T !== nothing ? _cast_curve(c, T) : c
    open(path, "w") do io
        _write_curve(fmt, io, c2)
    end
    path
end

"""
    load_curve(path; format=:auto, T=Float64, closed=nothing)

Load a single curve from `path`.

`T` controls the output floating-point type (default `Float64`).
`closed` overrides the topology stored in the file (useful for plain XYZ files
that have no topology metadata).

```julia
c = load_curve("curve.dcb")
c = load_curve("points.xyz"; closed=true)
c = load_curve("cloud.ply"; T=Float32)
```
"""
function load_curve(path::AbstractString;
                    format::Symbol  = :auto,
                    T::Type{<:AbstractFloat} = Float64,
                    closed::Union{Nothing, Bool} = nothing)
    fmt = format === :auto ? _detect_format(path) : format
    open(path, "r") do io
        c = _read_curve(fmt, io)
        c = _cast_curve(c, T)
        closed !== nothing && return DiscreteCurve(c.points; closed)
        c
    end
end

"""
    save_curves(path, curves; format=:auto, T=nothing)

Write a collection of curves to a single file.

For binary format (`.dcb`) all curves are packed in one file with a shared
header.  For text formats each curve is separated by a blank line (CSV / TSV /
XYZ) or a named element (PLY).

```julia
r = evolve(c, CurvatureFlow(); save_every=10)
save_curves("flow_snapshots.dcb", r.curves)
save_curves("flow_snapshots.csv", r.curves)
```
"""
function save_curves(path::AbstractString,
                     curves::AbstractVector{<:AbstractDiscreteCurve};
                     format::Symbol = :auto,
                     T::Union{Nothing, Type{<:AbstractFloat}} = nothing)
    isempty(curves) && error("curves is empty — nothing to write")
    fmt = format === :auto ? _detect_format(path) : format
    cs  = T !== nothing ? [_cast_curve(c, T) for c in curves] : collect(curves)
    open(path, "w") do io
        _write_curves(fmt, io, cs)
    end
    path
end

"""
    load_curves(path; format=:auto, T=Float64, closed=nothing)

Load a collection of curves from a file previously written with `save_curves`.

```julia
curves = load_curves("flow_snapshots.dcb")
curves = load_curves("data.csv"; closed=true)
```
"""
function load_curves(path::AbstractString;
                     format::Symbol  = :auto,
                     T::Type{<:AbstractFloat} = Float64,
                     closed::Union{Nothing, Bool} = nothing)
    fmt = format === :auto ? _detect_format(path) : format
    open(path, "r") do io
        cs = _read_curves(fmt, io)
        cs = [_cast_curve(c, T) for c in cs]
        if closed !== nothing
            cs = [DiscreteCurve(c.points; closed) for c in cs]
        end
        cs
    end
end


function _cast_curve(c::AbstractDiscreteCurve{N,T}, ::Type{U}) where {N, T, U}
    T === U && return c
    pts = [SVector{N,U}(p) for p in c.points]
    DiscreteCurve(pts; closed = isclosed(c))
end
_cast_curve(c, ::Nothing) = c

# ─────────────────────────────────────────────────────────────────────────────
#  BINARY FORMAT  (.dcb)
# ─────────────────────────────────────────────────────────────────────────────
#
#  Single-curve layout:
#    [6]  magic  "DCrv\x01\x00"
#    [1]  N      ambient dimension (UInt8)
#    [1]  T_code 0=Float32  1=Float64 (UInt8)
#    [1]  flags  bit0=closed (UInt8)
#    [1]  pad    reserved 0x00
#    [8]  n      nvertices (Int64 LE)
#    [N*n*sizeof(T)] vertex data, row-major (vertex-major: x1 y1 z1 x2 y2 z2 …)
#
#  Multi-curve layout:
#    [6]  magic  "DCrvm\x01"
#    [8]  nc     n_curves (Int64 LE)
#    [repeated nc times: single-curve block without leading magic]
#       [1] N, [1] T_code, [1] flags, [1] pad, [8] n, [data]

function _write_curve(::Val{:binary}, io::IO, c::AbstractDiscreteCurve{N,T}) where {N,T}
    haskey(T_CODE, T) || error("Binary format only supports Float32/Float64, got $T")
    write(io, DCB_MAGIC_SINGLE)
    write(io, UInt8(N))
    write(io, T_CODE[T])
    write(io, UInt8(isclosed(c) ? 0x01 : 0x00))
    write(io, UInt8(0x00))
    write(io, Int64(nvertices(c)))
    for p in c.points, coord in p
        write(io, coord)
    end
end
_write_curve(fmt::Symbol, io, c) = _write_curve(Val(fmt), io, c)

function _write_curves(::Val{:binary}, io::IO,
                       cs::AbstractVector{<:AbstractDiscreteCurve})
    write(io, DCB_MAGIC_MULTI)
    write(io, Int64(length(cs)))
    for c in cs
        T  = eltype(c.points[1])
        N  = length(c.points[1])
        haskey(T_CODE, T) || error("Binary format only supports Float32/Float64, got $T")
        write(io, UInt8(N))
        write(io, T_CODE[T])
        write(io, UInt8(isclosed(c) ? 0x01 : 0x00))
        write(io, UInt8(0x00))
        write(io, Int64(nvertices(c)))
        for p in c.points, coord in p
            write(io, coord)
        end
    end
end
_write_curves(fmt::Symbol, io, cs) = _write_curves(Val(fmt), io, cs)

function _read_header_body(io::IO)
    N_u     = read(io, UInt8)
    t_code  = read(io, UInt8)
    flags   = read(io, UInt8)
    _pad    = read(io, UInt8)
    n       = read(io, Int64)
    closed  = (flags & 0x01) != 0
    T       = haskey(CODE_T, t_code) ? CODE_T[t_code] :
              error("Unknown T_code $t_code in binary curve file")
    return (Int(N_u), T, Int(n), closed)  # N, T, n, closed
end

function _read_curve_body(io::IO, N::Int, T::Type, n::Int, closed::Bool)
    pts = Vector{SVector{N,T}}(undef, n)
    for j in 1:n
        coords = ntuple(_ -> read(io, T), Val(N))
        pts[j] = SVector{N,T}(coords)
    end
    DiscreteCurve(pts; closed)
end

function _read_curve(::Val{:binary}, io::IO)
    magic = read(io, 6)
    if magic == DCB_MAGIC_SINGLE
        N, T, n, closed = _read_header_body(io)
        return _read_curve_body(io, N, T, n, closed)
    elseif magic[1:5] == DCB_MAGIC_MULTI[1:5]
        # Multi-curve file — return first curve
        nc = read(io, Int64)
        nc == 0 && error("Empty multi-curve file")
        N, T, n, closed = _read_header_body(io)
        return _read_curve_body(io, N, T, n, closed)
    end
    error("Not a DCB file (magic mismatch)")
end
_read_curve(fmt::Symbol, io) = _read_curve(Val(fmt), io)

function _read_curves(::Val{:binary}, io::IO)
    magic = read(io, 6)
    if magic == DCB_MAGIC_SINGLE
        # Treat as single-curve file
        N, T, n, closed = _read_header_body(io)
        c = _read_curve_body(io, N, T, n, closed)
        return [c]
    elseif magic[1:5] == DCB_MAGIC_MULTI[1:5]
        nc = read(io, Int64)
        curves = AbstractDiscreteCurve[]
        for _ in 1:nc
            N, T, n, closed = _read_header_body(io)
            push!(curves, _read_curve_body(io, N, T, n, closed))
        end
        return curves
    end
    error("Not a DCB file (magic mismatch)")
end
_read_curves(fmt::Symbol, io) = _read_curves(Val(fmt), io)

# ─────────────────────────────────────────────────────────────────────────────
#  CSV FORMAT  (.csv)
# ─────────────────────────────────────────────────────────────────────────────
#  Header (comment lines):
#    # DiscreteCurves v1
#    # closed: true
#    # N: 2
#  Data rows: x,y[,z]   (comma-separated, no spaces)
#
#  Multi-curve: curves are separated by a blank line.

_csv_sep(::Val{:csv}) = ','
_csv_sep(::Val{:tsv}) = '\t'

function _write_curve_text(sep, io::IO, c::AbstractDiscreteCurve{N,T}) where {N,T}
    println(io, "# DiscreteCurves v1")
    println(io, "# closed: ", isclosed(c))
    println(io, "# N: $N")
    println(io, join(("x","y","z","w")[1:N], sep))
    for p in c.points
        println(io, join(p, sep))
    end
end

function _write_curve(::Val{:csv}, io, c) ; _write_curve_text(',', io, c) ; end
function _write_curve(::Val{:tsv}, io, c) ; _write_curve_text('\t', io, c) ; end

function _write_curves(::Val{:csv}, io, cs)
    for (k, c) in enumerate(cs)
        k > 1 && println(io)     # blank line separator
        _write_curve_text(',', io, c)
    end
end
function _write_curves(::Val{:tsv}, io, cs)
    for (k, c) in enumerate(cs)
        k > 1 && println(io)
        _write_curve_text('\t', io, c)
    end
end

function _parse_text_blocks(io::IO, sep::Char)
    blocks  = Vector{Tuple{Dict{String,String}, Vector{Vector{Float64}}}}()
    meta    = Dict{String,String}()
    rows    = Vector{Vector{Float64}}()
    header_done = false

    for raw_line in eachline(io)
        line = strip(raw_line)
        if startswith(line, "#")
            inner = strip(line[2:end])
            m = match(r"^(\w+):\s*(.*)", inner)
            m !== nothing && (meta[m.captures[1]] = strip(m.captures[2]))
            continue
        end
        if isempty(line)
            if !isempty(rows)
                push!(blocks, (copy(meta), copy(rows)))
                empty!(meta); empty!(rows)
                header_done = false
            end
            continue
        end
        parts = split(line, sep)
        if !header_done && any(c -> isletter(c) || c == '_', line)
            header_done = true
            continue
        end
        header_done = true
        vals = tryparse.(Float64, parts)
        all(!isnothing, vals) || continue
        push!(rows, collect(Float64, vals))
    end
    !isempty(rows) && push!(blocks, (copy(meta), copy(rows)))
    blocks
end

function _block_to_curve(meta, rows)
    isempty(rows) && error("Empty curve block")
    N      = length(rows[1])
    closed = get(meta, "closed", "false") == "true"
    pts    = [SVector{N,Float64}(r) for r in rows]
    DiscreteCurve(pts; closed)
end

function _read_curve(::Val{:csv}, io) ; _read_curve_text(io, ',') ; end
function _read_curve(::Val{:tsv}, io) ; _read_curve_text(io, '\t') ; end

function _read_curve_text(io, sep)
    blocks = _parse_text_blocks(io, sep)
    isempty(blocks) && error("No curve data found")
    _block_to_curve(blocks[1]...)
end

function _read_curves(::Val{:csv}, io) ; _read_curves_text(io, ',') ; end
function _read_curves(::Val{:tsv}, io) ; _read_curves_text(io, '\t') ; end

function _read_curves_text(io, sep)
    blocks = _parse_text_blocks(io, sep)
    isempty(blocks) && error("No curve data found")
    [_block_to_curve(m, r) for (m, r) in blocks]
end

# ─────────────────────────────────────────────────────────────────────────────
#  XYZ FORMAT  (.xyz / .txt)
# ─────────────────────────────────────────────────────────────────────────────
#  Plain whitespace-separated, one vertex per line, no header.
#  Topology (open/closed) is not stored — defaults to open on read.
#  Multi-curve: separated by blank lines.

function _write_curve(::Val{:xyz}, io::IO, c::AbstractDiscreteCurve)
    for p in c.points
        println(io, join(p, ' '))
    end
end

function _write_curves(::Val{:xyz}, io::IO, cs)
    for (k, c) in enumerate(cs)
        k > 1 && println(io)
        _write_curve(Val(:xyz), io, c)
    end
end

function _read_curve(::Val{:xyz}, io)
    blocks = _parse_text_blocks(io, ' ')
    isempty(blocks) && error("No XYZ data found")
    m, rows = blocks[1]
    pts = [SVector{length(r), Float64}(r) for r in rows]
    DiscreteCurve(pts; closed = false)
end

function _read_curves(::Val{:xyz}, io)
    blocks = _parse_text_blocks(io, ' ')
    isempty(blocks) && error("No XYZ data found")
    [let N = length(r[1]); DiscreteCurve([SVector{N,Float64}(v) for v in r]; closed=false)
     end for (_, r) in blocks]
end

# ─────────────────────────────────────────────────────────────────────────────
#  PLY FORMAT  (.ply)  — ASCII only
# ─────────────────────────────────────────────────────────────────────────────
#  Standard ASCII PLY with one polyline element.
#  Each curve is stored as a face element with vertex_indices.

function _write_curve(::Val{:ply}, io::IO, c::AbstractDiscreteCurve{N,T}) where {N,T}
    n    = nvertices(c)
    props = N ≥ 3 ? ("x","y","z") : ("x","y")
    props = props[1:min(N,3)]

    println(io, "ply")
    println(io, "format ascii 1.0")
    println(io, "comment DiscreteCurves v1")
    println(io, "comment closed ", isclosed(c))
    println(io, "element vertex $n")
    for p in props; println(io, "property double $p"); end
    println(io, "element edge 1")
    println(io, "property list int int vertex_indices")
    println(io, "end_header")
    for p in c.points
        println(io, join(p[1:min(N,3)], ' '))
    end
    indices = isclosed(c) ? vcat(0:n-1, [0]) : collect(0:n-1)
    println(io, length(indices), ' ', join(indices, ' '))
end

function _write_curves(::Val{:ply}, io::IO, cs)
    n_curves = length(cs)
    all_pts  = SVector[]
    offset   = 0
    edge_lists = Vector{Vector{Int}}()

    for c in cs
        n = nvertices(c)
        for p in c.points; push!(all_pts, p); end
        idx = isclosed(c) ? vcat(offset:offset+n-1, [offset]) : collect(offset:offset+n-1)
        push!(edge_lists, idx)
        offset += n
    end

    N    = length(cs[1].points[1])
    props = ("x","y","z")[1:min(N,3)]

    println(io, "ply")
    println(io, "format ascii 1.0")
    println(io, "comment DiscreteCurves v1 multi")
    println(io, "comment n_curves $n_curves")
    println(io, "element vertex $(length(all_pts))")
    for p in props; println(io, "property double $p"); end
    println(io, "element edge $n_curves")
    println(io, "property list int int vertex_indices")
    println(io, "end_header")
    for p in all_pts
        println(io, join(p[1:min(N,3)], ' '))
    end
    for el in edge_lists
        println(io, length(el), ' ', join(el, ' '))
    end
end

function _read_curve(::Val{:ply}, io)
    cs = _read_curves(Val(:ply), io)
    isempty(cs) && error("No curves in PLY file")
    cs[1]
end

function _read_curves(::Val{:ply}, io)
    n_verts = 0; n_edges = 0; prop_names = String[]; in_header = true
    closed_flag = false; n_curves_meta = 0

    for line in eachline(io)
        s = strip(line)
        s == "end_header" && (in_header = false; break)
        startswith(s, "comment closed ")  && (closed_flag = strip(s[16:end]) == "true")
        startswith(s, "comment n_curves") && (n_curves_meta = parse(Int, split(s)[end]))
        startswith(s, "element vertex ")  && (n_verts = parse(Int, split(s)[end]))
        startswith(s, "element edge ")    && (n_edges = parse(Int, split(s)[end]))
        startswith(s, "property double ") && push!(prop_names, split(s)[end])
    end

    in_header && error("PLY file missing end_header")
    n_verts == 0 && error("PLY file has no vertex element")
    N = length(prop_names)
    N == 0 && (N = 3)

    all_pts = Vector{SVector{N, Float64}}(undef, n_verts)
    for i in 1:n_verts
        parts = split(strip(readline(io)))
        all_pts[i] = SVector{N, Float64}(parse.(Float64, parts[1:N])...)
    end

    n_edges == 0 && return [DiscreteCurve(all_pts; closed = closed_flag)]

    curves = AbstractDiscreteCurve[]
    for _ in 1:n_edges
        parts = split(strip(readline(io)))
        cnt   = parse(Int, parts[1])
        idxs  = parse.(Int, parts[2:end])
        cl = length(idxs) > cnt && idxs[end] == idxs[1]
        unique_idxs = cl ? idxs[1:end-1] : idxs
        pts = [all_pts[i+1] for i in unique_idxs]  # +1: PLY is 0-based
        push!(curves, DiscreteCurve(pts; closed = cl || closed_flag))
    end
    curves
end

"""
    curve_info(path)  →  String

Print metadata for a curve file without loading the full data.

```julia
curve_info("large_flow.dcb")
```
"""
function curve_info(path::AbstractString)
    fmt = _detect_format(path)
    open(path, "r") do io
        if fmt === :binary
            magic = read(io, 6)
            if magic == DCB_MAGIC_SINGLE
                N, T, n, closed = _read_header_body(io)
                return "DCB single curve: N=$N T=$T n=$n closed=$closed  " *
                       "($(filesize(path)) bytes)"
            elseif magic[1:5] == DCB_MAGIC_MULTI[1:5]
                nc = read(io, Int64)
                return "DCB multi-curve: $nc curves  ($(filesize(path)) bytes)"
            end
            return "DCB file (unknown variant)"
        else
            return "$(uppercase(String(fmt))) file at $path  ($(filesize(path)) bytes)"
        end
    end
end

end