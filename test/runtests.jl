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
sv_01 = SVertex(Dict(0x00000001=>0x00000001, 0x00000002=>0x00000002), Int64.([4,5]),0x00000005,0x00000000,0,0)
sv_02 = SVertex(Dict(0x00000001=>0x00000001, 0x00000002=>0x00000002, 0x00000005=>0x00000001, 0x00000009=>0x00000002), Int64.([4,5]),0x00000005,0x00000000,0,0)

fullMap = identityMap(-5, 5, -5, 4, 0, offset, base)
freq_red, Fup_red_data, Fdo_red_data  = read_vert_chi("data/vert_chi_red")   
freq_full, Fup_full_data, Fdo_full_data = read_vert_chi("data/vert_chi_full")   
gImp = read_gm_wim(nBose+2*nFermi, "data/gm_wim", storedInverse=false)
sv_up_test_red = SVertex(redMap, Fup_red_data, base, offset, nBose, nFermi)
sv_do_test_red = SVertex(redMap, Fdo_red_data, base, offset, nBose, nFermi)
sv_up_test_full = SVertex(fullMap, Fup_full_data, base, offset, nBose, nFermi)
sv_do_test_full = SVertex(fullMap, Fdo_full_data, base, offset, nBose, nFermi)

@testset "SparseVertex" begin
    @testset "Vertex" begin
        include("SVertex.jl")   
    end
    @testset "IO" begin
        include("IO.jl")   
    end
    @testset "SVertexTools" begin
        include("SVertexTools.jl")   
    end

end
