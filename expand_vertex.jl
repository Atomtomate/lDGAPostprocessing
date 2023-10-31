using Pkg
Pkg.activate(@__DIR__)
using lDGAPostprocessing, jED

using EquivalenceClassesConstructor
using Test
using DataStructures
using DelimitedFiles, JLD2
using TOML
include("./vertexIntTriple.jl")

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
    χ_upup, χ_updo = expand(TwoPartGF_upup_red, TwoPartGF_updo_red, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);
    return χ_upup, χ_updo
end

gridPath = ARGS[1]  # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
dataPath = ARGS[2] # "/scratch/usr/hhpstobb/lDGA/tests/non_qubic/ed_vertex"

println("opening: ", gridPath*"/freqList.jld2")
jldopen(gridPath*"/freqList.jld2", "r") do f
    for k in ["freqList", "nFermi", "nBose", "shift"] #keys(f)
        s=Symbol(k)
        @eval (($s) = ($(f[k])))
    end
end
cfg         = TOML.parsefile(joinpath(dataPath, "config.toml"))
U           = cfg["parameters"]["U"]
β           = cfg["parameters"]["beta"]
ϵₖ, Vₖ, μ   = read_anderson_parameters(joinpath(dataPath, "hubb.andpar"))
include("gen_GF.jl")

# ======================================== Unpack 2-Part-GF ========================================

println("Expanding Vertex for nFermi=",nFermi,", nBose=",nBose,", shift=",shift)
flush(stderr)
flush(stdout)


χ_upup, χ_updo = expand_TwoPartGF(gridPath*"/freqList.jld2", dataPath)
println("Done expanding!")
flush(stderr)
flush(stdout)


# This segment computes quntities in the ph channel
# We first subtract the unconnected part of the susceptibility
χ0_full = lDGAPostprocessing.computeχ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp.parent, β)
lDGAPostprocessing.add_χ₀_ω₀!(freqList, χ_upup, GImp.parent, β)
lDGAPostprocessing.add_χ₀_ω₀!(freqList, χ_updo, GImp.parent, β)
Γsp, Γch = -1.0 .* computeΓ_ph(freqList, χ_upup .- χ_updo,  χ_upup .+ χ_updo, χ0_full,nBose,nFermi)



println("Done calculating vertex!")
flush(stderr)
flush(stdout)

E_kin_DMFT, E_pot_DMFT  = calc_E_ED(νnGrid[0:last(axes(GImp,1))], ϵₖ, Vₖ, GImp.parent, nden, U, β, μ)
res = isfile(dataPath * "/chi_asympt") ? read_chi_asympt(dataPath * "/chi_asympt") : error("chi_asympt not found!")
χ_ch_asympt, χ_sp_asympt, χ_pp_asympt = res

jldopen(dataPath*"/ED_out.jld2", "w") do f
    f["Γch"] = permutedims(Γch, [3,1,2])
    f["Γsp"] = permutedims(Γsp, [3,1,2])
    f["χDMFTch"] = permutedims(reshape(χ_upup .+ χ_updo, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
    f["χDMFTsp"] = permutedims(reshape(χ_upup .- χ_updo, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
    f["χ_ch_asympt"] = χ_ch_asympt
    f["χ_sp_asympt"] = χ_sp_asympt
    f["χ_pp_asympt"] = χ_pp_asympt
    f["gImp"] = GImp.parent
    f["g0"] = G0W.parent
    f["ϵₖ"] = ϵₖ
    f["Vₖ"] = Vₖ
    f["μ"] = μ
    f["U"] = U
    f["β"] = β
    f["nden"] = nden
    f["E_kin_DMFT"] = E_kin_DMFT
    f["E_pot_DMFT"] = E_pot_DMFT
    f["grid_shift"] = shift
    f["grid_nBose"] = nBose
    f["grid_nFermi"] = nFermi
end
