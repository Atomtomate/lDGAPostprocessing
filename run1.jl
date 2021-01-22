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
sv_01 = SVertex(Dict(0x00000001=>0x00000001, 0x00000002=>0x00000002), Int64.([4,5]),0x00000005,0x00000000,0,0)
sv_02 = SVertex(Dict(0x00000005=>0x00000001, 0x00000009=>0x00000002), Int64.([4,5]),0x00000005,0x00000000,0,0)
sv_read_1 = SVertex(redMap, testData_red, base, offset, nBose, nFermi)
sv_read_2 = SVertex(redMap, testData_red, base, offset, nBose, nFermi)


fullMap = identityMap(-5, 5, -5, 4, 0, offset, base)
freq_red, Fup_red_data, Fdo_red_data = read_vert_chi("test/data/vert_chi_red")   
freq_full, Fup_full_data, Fdo_full_data = read_vert_chi("test/data/vert_chi_full")   

sv_up_test_red = SVertex(redMap, Fup_red_data, base, offset, nBose, nFermi)
sv_do_test_red = SVertex(redMap, Fdo_red_data, base, offset, nBose, nFermi)
sv_up_test_full = SVertex(fullMap, Fup_full_data, base, offset, nBose, nFermi)
sv_do_test_full = SVertex(fullMap, Fdo_full_data, base, offset, nBose, nFermi)

gImp = read_gm_wim(nBose+2*nFermi, "test/data/gm_wim", storedInverse=false)
χ0 = SparseVertex.computeχ0(-5:5, -5:4, gImp, 25.0)
F_den, F_mag = SparseVertex.computeF(freq_full, sv_up_test_full, sv_do_test_full, χ0)


ind = indices(sv_up_test_full) 
