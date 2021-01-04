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
sp_read = SVertex(redMap, testData_red, base, offset)
