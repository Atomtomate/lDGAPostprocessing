using Pkg
Pkg.activate(@__DIR__)
using SparseVertex

using EquivalenceClassesConstructor
using Test
using DataStructures
using DelimitedFiles, JLD2
include("./vertexIntTriple.jl")

gridPath = ARGS[1]  # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
dataPath = ARGS[2] # "/scratch/usr/hhpstobb/lDGA/tests/non_qubic/ed_vertex"
β = parse(Float64, ARGS[3])

println("opening: ", gridPath*"/freqList.jld2")
f = jldopen(gridPath*"/freqList.jld2", "r")
for k in keys(f)
    s=symbol(k)
    @eval (($s) = ($(f[k])))
end

#@load gridPath*"/freqList.jld2" freqRed_map freqList freqList_min parents ops nFermi nBose shift base offset

println("Expanding Vertex for nFermi=",nFermi,", nBose=",nBose,", shift=",shift)
base2 = base*base
freq_red, TwoPartGF_upup_red, TwoPartGF_updo_red = read_vert_chi(dataPath*"/2_part_gf_red") 

gImp = read_gm_wim(2*(nBose+nFermi+1), dataPath*"/gm_wim", storedInverse=false)
g0 = 1.0 ./ read_gm_wim(2*(nBose+nFermi+1), dataPath*"/g0mand", storedInverse=false)
χ0_full = SparseVertex.computeχ0(-nBose:nBose, -(nFermi+shift*nBose):(nFermi+shift*nBose)-1, gImp, β)
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

TwoPartGF_upup, TwoPartGF_updo = expand(TwoPartGF_upup_red, TwoPartGF_updo_red, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);

SparseVertex.subtract_ω0!(freqList, TwoPartGF_upup, gImp, β)
SparseVertex.subtract_ω0!(freqList, TwoPartGF_updo, gImp, β)
χDMFTch = TwoPartGF_upup .+ TwoPartGF_updo
χDMFTsp = TwoPartGF_upup .- TwoPartGF_updo
Γch = -1.0 .* SparseVertex.computeχ(freqList, χDMFTch, χ0_full,nBose,nFermi)
Γsp = -1.0 .* SparseVertex.computeχ(freqList, χDMFTsp, χ0_full,nBose,nFermi)

#SparseVertex.write_fort_dir("gamma", freqList, Γch, Γsp, dataPath*"/gamma_dir", nBose, nFermi)
#SparseVertex.write_fort_dir("chi", freqList, TwoPartGF_upup, TwoPartGF_updo, dataPath*"/chi_dir", nBose, nFermi)
#SparseVertex.write_fort_dir("chi", freqList, χDMFTch, χDMFTsp, dataPath*"/chi_dir", nBose, nFermi)

χDMFTch = permutedims(reshape(χDMFTch, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTsp = permutedims(reshape(χDMFTsp, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])

ϵₖ, Vₖ, μ    = read_anderson_parameters(dataPath * "/hubb.andpar");
U, β, nden   = read_hubb_dat(dataPath * "/hubb.dat");
iνₙ, GImp_ED = readGImp(dataPath*"/gm_wim", only_positive=true)
E_kin_DMFT, E_pot_DMFT  = calc_E_ED(iνₙ[1:length(GImp_ED)], ϵₖ, Vₖ, GImp_ED, nden, U, β, μ)

@save dataPath*"/ED_out.jld2" Γch Γsp χDMFTch χDMFTsp gImp g0 ϵₖ Vₖ μ U β nden E_kin_DMFT E_pot_DMFT
