using Pkg
Pkg.activate(@__DIR__)
using SparseVertex

using EquivalenceClassesConstructor
using Test
using DataStructures
using JLD2
include("./vertexIntTriple.jl")

gridPath = ARGS[1]  # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
dataPath = ARGS[2] # "/scratch/usr/hhpstobb/lDGA/tests/non_qubic/ed_vertex"
#dataPath = "test/data"
@load gridPath*"/freqList.jld2" freqRed_map freqList freqList_min parents ops nFermi nBose shift base offset
println("Expanding Vertex for nFermi=",nFermi,", nBose=",nBose,", shift=",shift)
base2 = base*base
freq_red, Fup_red_data, Fdo_red_data = read_vert_chi(dataPath*"/vert_chi") 
#fullMap = identityMap(-nBose, nBose, -nFermi, nFermi-1, 0, offset, base)
#freq_full, Fup_full_data, Fdo_full_data = read_vert_chi("test/data/vert_chi_full")   
#freq_full, Fup_full_data, Fdo_full_data = read_vert_chi(dataPath*"/../50_full_nshift/vert_chi")   
#sv_up_test_full = SVertex(fullMap, Fup_full_data, base, offset, nBose, nFermi)
#sv_do_test_full = SVertex(fullMap, Fdo_full_data, base, offset, nBose, nFermi)
#χ_full_data = 

gImp = read_gm_wim(nBose+2*nFermi, dataPath*"/gm_wim", storedInverse=false)
g0 = 1.0 ./ read_gm_wim(nBose+2*nFermi, dataPath*"/g0mand", storedInverse=false)
χ0_full = SparseVertex.computeχ0(-nBose:nBose, -(nFermi+shift*nBose):(nFermi+shift*nBose)-1, gImp, 25.0)
#F_den, F_mag = SparseVertex.computeF(freq_full, sv_up_test_full, sv_do_test_full, χ0_full)
#ind = indices(sv_up_test_full) 

transform_f(F_up::Array{T,1}, F_do::Array{T,1}, prev::Int, next::Int, ops::Array{Int,1}) where T = 
                               if ops[next] == 1 
                                   F_up[next] = conj(F_up[prev])
                                   F_do[next] = conj(F_do[prev])
                               elseif ops[next] == 3 || ops[next] == 4
                                   F_up[next] = -F_up[prev]
                                   F_do[next] = F_do[prev] - F_up[prev]
                               else
                                   F_up[next] = F_up[prev]
                                   F_do[next] = F_do[prev]
                               end

#Fup_exp_t, Fdo_exp_t = SparseVertex.expand_test(Fup_red_data, Fdo_red_data, freqRed_map, freqList, parents, ops, nBose, nFermi);
χup_exp,χdo_exp = expand(Fup_red_data, Fdo_red_data, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);
χDMFTch = χup_exp .+ χdo_exp
χDMFTsp = χup_exp .- χdo_exp
Γch = -1.0 .* SparseVertex.computeχ(freqList, χDMFTch, χ0_full,nBose,nFermi)
Γsp = -1.0 .* SparseVertex.computeχ(freqList, χDMFTsp, χ0_full,nBose,nFermi)
#SparseVertex.write_fort_dir("chi", freqList, χch, χsp, dataPath*"/chi_dir", nBose, nFermi)
#SparseVertex.write_fort_dir("gamma", freqList, Fup_exp, Fdo_exp, dataPath*"/gamma_dir", nBose, nFermi)
SparseVertex.write_fort_dir("gamma", freqList, Γch, Γsp, dataPath*"/gamma_dir", nBose, nFermi)
SparseVertex.write_fort_dir("chi", freqList, χup_exp, χdo_exp, dataPath*"/chi_dir", nBose, nFermi)
χDMFTch = permutedims(reshape(χDMFTch, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTsp = permutedims(reshape(χDMFTsp, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])

@save dataPath*"/ED_out.jld2" Γch Γsp χDMFTch χDMFTsp gImp g0
