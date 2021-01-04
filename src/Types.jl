import Base: show, display
import Base: getindex, setindex!, similar, eachindex
import Base: eltype, length, size
#import Base: endof, indices, start, next, done, iteratorsize, iteratoreltype


"""
    SVertex
Struct
"""
struct SVertex{T} <: AbstractArray{T,1}
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


# ................... direct ...................
getindex(sv::SVertex, i::UInt32, j::UInt32, k::UInt32) = sv.data[sv.map[tripleToInt(i,j,k,sv.offset,sv.nB,sv.nB2)]]
getindex(sv::SVertex, i::UInt32) = sv.data[sv.map[i]]
getindex(sv::SVertex, i::Number) = sv.data[sv.map[UInt32(i)]]
getindex(sv::SVertex, I) = [sv.data[sv.map[UInt32(i)]] for i in I]

# ... setindex!

# .................. abstract ..................

size(sv::SVertex) = (sv.n,)

# -------------------- IO ----------------------

function Base.show(io::IO, ::MIME"text/plain", sv::SVertex{T}) where {T}
    println(io, "SparseVertex{$T} with $n elements. Encoded with base $(sv.nB), offset $(sv.offset)")
           for i in eachindex(v)
               println(io, "  ", v.v[i])
           end
       end
