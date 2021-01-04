import Base: show, display, @propagate_inbounds
import Base: getindex, setindex!, similar, eachindex
import Base: eltype, length, size, iterate
#import Base: endof, indices, start, next, done, iteratorsize, iteratoreltype


"""
    SVertex
Struct
"""
struct SVertex{T} <: AbstractDict{T,1}
    n::Int              # length
    nB::UInt32          # triple-to-int parameter
    nB2::UInt32         # triple-to-int parameter
    offset::UInt32      # triple-to-int parameter
    data::Array{T,1}    # independent data segments
    map::Dict{UInt32,UInt32}  # map from full range to data
    SVertex(map::Dict{UInt32,UInt32}, data::Array{T,1}, base::UInt32, offset::UInt32) where {T} =
        new{T}(length(data), base, base*base, offset, data, map)
end

"""
    SVertex(redM::ReduceMap, data, base)
"""
SVertex(redm::ReduceMap, data, base, offset) = SVertex(redm.map, data, UInt32(base), UInt32(offset))

# ------------------ Indexing ------------------
IndexStyle(::SVertex) = IndexLinear()

@propagate_inbounds getindex(sv::SVertex, i::Int, j::Int, k::Int) = 
    get(sv.data, sv.map[tripleToInt(i,j,k,sv.offset,sv.nB,sv.nB2)], missing)
@propagate_inbounds getindex(sv::SVertex, i::UInt32) = 
    get(sv.data, sv.map[i], missing)
@propagate_inbounds getindex(sv::SVertex, i::Number) = 
    get(sv.data, sv.map[UInt32(i)], missing)
@propagate_inbounds getindex(sv::SVertex, I) = [sv[i] for i in I]

@propagate_inbounds setindex!(sv::SVertex{T}, val::T, i::UInt32, j::UInt32, k::UInt32) where T = 
    sv.data[sv.map[tripleToInt(i,j,k,sv.offset,sv.nB,sv.nB2)]] = val
#@propagate_inbounds setindex!(sv::SVertex{T}, val::T, i::UInt32) = sv.data[sv.map[i]]
#@propagate_inbounds setindex!(sv::SVertex{T}, val::T, i::Number) = sv.data[sv.map[UInt32(i)]]
#@propagate_inbounds setindex!(sv::SVertex{T}, val::T, I) = [sv.data[sv.map[UInt32(i)]] for i in I]
# ... setindex!


# ----------------- Iteration ------------------
@propagate_inbounds eachindex(sv::SVertex) = keys(sv.map)
@propagate_inbounds iterate(sv::SVertex) = iterate(sv.data)
@propagate_inbounds iterate(sv::SVertex, i) = iterate(sv.data, i)
eltype(sv::SVertex{T}) where T = T
length(sv::SVertex) = sv.n
size(sv::SVertex) = (sv.n,)

# -------------------- IO ----------------------

# Adapted from dict.jl
function Base.show(io::IO, ::MIME"text/plain", sv::SVertex{T}) where {T}
    recur_io = IOContext(io, :SHOWN_SET => sv,
                             :typeinfo => eltype(sv))
    n = 0
    limit = get(io, :limit, false)::Bool
    println(io, "SparseVertex{$T} with $n elements. Encoded with base $(sv.nB), offset $(sv.offset)")
    for i in eachindex(sv.data)
        println(io, "  ", sv.data[i])
        limit && n >= 10 && (print(io, "…"); break)
    end
end


@fastmath @inline tripleToInt(i::UInt32, j::UInt32, k::UInt32, nB::UInt32, nB2::UInt32)::UInt32 = 
    UInt32(i*nB2 + j*nB + k)
@fastmath @inline tripleToInt(i::Int, j::Int, k::Int, offset::UInt32, nB::UInt32, nB2::UInt32)::UInt32 =
    UInt32((i+offset)*nB2 + (j+offset)*nB + (k+offset))

@fastmath @inline function intToTriple(::Type{T}, z::UInt32) where {T<:Integer}
    r,k = divrem(z,nB)
    i,j = divrem(r,nB)
    return (convert(T,i)-nBh,convert(T,j)-nBh,convert(T,k)-nBh)
end
intToTriple(z::UInt32) = intToTriple(Int64, z)
