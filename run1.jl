using SparseVertex

using EquivalenceClassesConstructor
using Test
using JLD
include("../vertexIntTriple.jl")

data = load("test/data/freqList.jld")
expMap = data["ExpandMap"]
redMap = data["ReduceMap"]
base = data["base"]
nFermi = data["nFermi"]
nBose = data["nBose"]
shift = data["shift"]
offset = data["offset"]
testData_full = [i for i in 1:length(redMap)]
testData_red  = [i for i in 1:length(expMap)]
#sv_01 = SVertex(Dict(3=>1, 4=>2), Int64.([10,11]))
sv_01 = SVertex(Dict(0x00000001=>0x00000001, 0x00000002=>0x00000002), Int64.([4,5]),0x00000005,0x00000000)
sv_02 = SVertex(Dict(0x00000005=>0x00000001, 0x00000009=>0x00000002), Int64.([4,5]),0x00000005,0x00000000)
sv_read_1 = SVertex(redMap, testData_red, base, offset)
sv_read_2 = SVertex(redMap, testData_red, base, offset)
#for (i,iel) in enumerate(eachindex(sv_read_1))
#    sv_read_2[iel] = testData_red[i]
#end
