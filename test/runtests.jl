using SparseVertex

using EquivalenceClassesConstructor
using Test
using JLD


data = load("data/freqList.jld")
expMap = data["ExpandMap"]
redMap = data["ReduceMap"]
base = data["base"]
nFermi = data["nFermi"]
nBose = data["nBose"]
shift = data["shift"]
offset = data["offset"]
testData_full = [i for i in 1:length(redMap)]
testData_red  = [i for i in 1:length(expMap)]
sv_read_1 = SVertex(redMap, testData_red, base, offset)
sv_read_2 = SVertex(redMap, testData_red, base, offset)
sv_01 = SVertex(Dict(0x00000001=>0x00000001, 0x00000002=>0x00000002), Int64.([4,5]),0x00000005,0x00000000)
sv_02 = SVertex(Dict(0x00000005=>0x00000001, 0x00000009=>0x00000002), Int64.([4,5]),0x00000005,0x00000000)

@testset "SparseVertex" begin


    @testset "Vertex" begin
        include("Types.jl")   
    end
    @testset "IO" begin
        include("IO.jl")   
    end
    @testset "SVertexTools" begin
        include("SVertexTools.jl")   

end
