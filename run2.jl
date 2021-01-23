using Pkg
Pkg.activate(".")
using SparseVertex

using EquivalenceClassesConstructor
using Test
using DataStructures
using JLD2
include("../vertexIntTriple.jl")

@load "test/data/freqList_new_2.jld2" freqRed_map freqList freqList_min parents ops nFermi nBose shift base offset
base2 = base*base

fullMap = identityMap(-nBose, nBose, -nFermi, nFermi-1, 0, offset, base)
freq_red, Fup_red_data, Fdo_red_data = read_vert_chi("test/data/vert_chi_red_2") 
freq_full, Fup_full_data, Fdo_full_data = read_vert_chi("test/data/vert_chi_full_2")   
sv_up_test_full = SVertex(fullMap, Fup_full_data, base, offset, nBose, nFermi)
sv_do_test_full = SVertex(fullMap, Fdo_full_data, base, offset, nBose, nFermi)
#χ_full_data = 

gImp = read_gm_wim(nBose+2*nFermi, "test/data/gm_wim", storedInverse=false)
χ0_full = SparseVertex.computeχ0(-nBose:nBose, -nFermi:nFermi-1, gImp, 25.0)
F_den, F_mag = SparseVertex.computeF(freq_full, sv_up_test_full, sv_do_test_full, χ0_full)
ind = indices(sv_up_test_full) 

@inline transform_f(val, op) = if op == 1 
                                return conj(val)
                               elseif op == 3
                                return -val
                               elseif op == 4
                                return -val
                               else
                                return val
                               end

Fup_exp_t, Fdo_exp_t = SparseVertex.expand_test(Fup_red_data, Fdo_red_data, freqRed_map, freqList, parents, ops, nBose, nFermi);
Fup_exp, Fdo_exp = expand(Fup_red_data, Fdo_red_data, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);
chi_den = SparseVertex.computeχ(Fup_exp .+ Fdo_exp,χ0_full,nBose,nFermi)
chi_mag = SparseVertex.computeχ(Fup_exp .- Fdo_exp,χ0_full,nBose,nFermi)

