using SparseVertex

using EquivalenceClassesConstructor
using Test
using JLD

@testset "SparseVertex" begin

    data = load("data/freqList.jld")
    expMap = data["ExpandMap"]
    redMap = data["ReduceMap"]
    nFermi = data["nFermi"]
    nBose = data["nBose"]
    shift = data["shift"]
    testData_full  = [i for i in 1:length(redMap)]
    testData_red  = [i for i in 1:length(expMap)]
    sv_01 = SVertex(Dict(3=>1, 4=>2), Int64.([10,11]))
    sp_read = SVertex{Int64,Float64}(redMap, testData_red)


    @testset "Vertex" begin
        include("Types.jl")   
    end
    @testset "IO" begin
        include("IO.jl")   
    end

end
