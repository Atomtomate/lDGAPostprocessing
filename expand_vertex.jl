using Pkg
Pkg.activate(@__DIR__)
using lDGAPostprocessing

using EquivalenceClassesConstructor
using Test
using DataStructures
using DelimitedFiles, JLD2
include("./vertexIntTriple.jl")

gridPath = ARGS[1]  # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
dataPath = ARGS[2] # "/scratch/usr/hhpstobb/lDGA/tests/non_qubic/ed_vertex"
β = parse(Float64, ARGS[3])

println("opening: ", gridPath*"/freqList.jld2")
jldopen(gridPath*"/freqList.jld2", "r") do f
    for k in ["freqList", "nFermi", "nBose", "shift"] #keys(f)
        s=Symbol(k)
        @eval (($s) = ($(f[k])))
    end
end

#@load gridPath*"/freqList.jld2" freqRed_map freqList freqList_min parents ops nFermi nBose shift base offset
println("Expanding Vertex for nFermi=",nFermi,", nBose=",nBose,", shift=",shift)
flush(stderr)
flush(stdout)

function expand_TwoPartGF(freqFile, dataPath)

    jldopen(gridPath*"/freqList.jld2", "r") do f
        for k in ["freqRed_map", "parents", "ops"]
            s=Symbol(k)
            @eval (($s) = ($(f[k])))
        end
    end

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
    freq_red, TwoPartGF_upup_red, TwoPartGF_updo_red = read_vert_chi(dataPath*"/2_part_gf_red") 
    TwoPartGF_upup, TwoPartGF_updo = expand(TwoPartGF_upup_red, TwoPartGF_updo_red, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);
    return TwoPartGF_upup, TwoPartGF_updo
end

TwoPartGF_upup, TwoPartGF_updo = expand_TwoPartGF(gridPath*"/freqList.jld2", dataPath)
println("Done expanding!")
flush(stderr)
flush(stdout)

gImp    = read_gm_wim(2*(nBose+nFermi+1), dataPath*"/gm_wim", storedInverse=false)
g0      = 1.0 ./ read_gm_wim(2*(nBose+nFermi+1), dataPath*"/g0mand", storedInverse=false)
χ0_full = lDGAPostprocessing.computeχ0(-nBose:nBose, -(nFermi+shift*nBose):(nFermi+shift*nBose)-1, gImp, β)
#F_den, F_mag = lDGAPostprocessing.computeF(freq_full, sv_up_test_full, sv_do_test_full, χ0_full)
#ind = indices(sv_up_test_full) 


χDMFTch = TwoPartGF_upup .+ TwoPartGF_updo
χDMFTsp = TwoPartGF_upup .- TwoPartGF_updo
Γch     = -1.0 .* lDGAPostprocessing.computeΓ(freqList, χDMFTch, χ0_full,nBose,nFermi)
Γsp     = -1.0 .* lDGAPostprocessing.computeΓ(freqList, χDMFTsp, χ0_full,nBose,nFermi)
lDGAPostprocessing.add_χ₀!(freqList, TwoPartGF_upup, gImp, β)
lDGAPostprocessing.add_χ₀!(freqList, TwoPartGF_updo, gImp, β)

println("Done calculating vertex!")
flush(stderr)
flush(stdout)
# TODO: activate this via write_fortran flag
lDGAPostprocessing.write_vert_chi(freqList, TwoPartGF_upup, TwoPartGF_updo, dataPath, 2*nBose+1, 2*nFermi)
lDGAPostprocessing.write_fort_dir("gamma", freqList, -Γch, -Γsp, dataPath*"/gamma_dir", 2*nBose+1, 2*nFermi)
lDGAPostprocessing.write_fort_dir("chi", freqList, TwoPartGF_upup, TwoPartGF_updo, dataPath*"/chi_dir", 2*nBose+1, 2*nFermi)
# TODO: find a way to keep memory consumption low
TwoPartGF_updo = nothing
TwoPartGF_upup = nothing

χDMFTch = permutedims(reshape(χDMFTch, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTsp = permutedims(reshape(χDMFTsp, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])

ϵₖ, Vₖ, μ    = read_anderson_parameters(dataPath * "/hubb.andpar");
U, β, nden   = read_hubb_dat(dataPath * "/hubb.dat");
iνₙ, GImp_ED = readGImp(dataPath*"/gm_wim", only_positive=true)
E_kin_DMFT, E_pot_DMFT  = calc_E_ED(iνₙ[1:length(GImp_ED)], ϵₖ, Vₖ, GImp_ED, nden, U, β, μ)
res = isfile(dataPath * "/chi_asympt") ? read_chi_asympt(dataPath * "/chi_asympt") : ([], [], [])
χ_ch_asympt, χ_sp_asympt, χ_pp_asympt = res

gImp = read_gm_wim(10*(nBose+nFermi+1), dataPath*"/gm_wim", storedInverse=false)
g0 = 1.0 ./ read_gm_wim(10*(nBose+nFermi+1), dataPath*"/g0mand", storedInverse=false)

jldopen(dataPath*"/ED_out.jld2", "w") do f
    f["Γch"] = Γch
    f["Γsp"] = Γsp
    f["χDMFTch"] = χDMFTch
    f["χDMFTsp"] = χDMFTsp
    f["χ_ch_asympt"] = χ_ch_asympt
    f["χ_sp_asympt"] = χ_sp_asympt
    f["χ_pp_asympt"] = χ_pp_asympt
    f["gImp"] = gImp
    f["g0"] = g0
    f["ϵₖ"] = ϵₖ
    f["Vₖ"] = Vₖ
    f["μ"] = μ
    f["U"] = U
    f["β"] = β
    f["nden"] = nden
    f["E_kin_DMFT"] = E_kin_DMFT
    f["E_pot_DMFT"] = E_pot_DMFT
    f["FUpDo"] = FUpDo_from_χDMFT(0.5 .* (χDMFTch .- χDMFTsp), gImp, gridPath*"/freqList.jld2", β)
    f["grid_shift"] = shift
    f["grid_nBose"] = nBose
    f["grid_nFermi"] = nFermi
end
