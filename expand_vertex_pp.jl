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

gridPath = ARGS[1]    # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
dataPath = ARGS[2]    # "/scratch/usr/hhpstobb/lDGA/tests/non_qubic/ed_vertex"

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


χph_upup, χph_updo = expand_TwoPartGF(gridPath*"/freqList.jld2", dataPath)
lDGAPostprocessing.add_χ₀_ω₀!(freqList, χph_upup, GImp.parent, β)
lDGAPostprocessing.add_χ₀_ω₀!(freqList, χph_updo, GImp.parent, β)

lDGAPostprocessing.write_vert_chi(freqList, χph_upup, χph_updo, ".", nBose, nFermi)
println("Done expanding!")
flush(stderr)
flush(stdout)


# This segment computes the reducible vertex in the pp channel
# - Compute χ_s
# - Compute F_s
# - Compute Γ_s
# - Compute ϕ_s
χ0_pp_full   = computeχ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp.parent, β; mode=:pp)

χpp_s, χpp_t = χph_to_χpp(freqList, χph_upup, χph_updo, χ0_pp_full, shift, nBose, nFermi)
#Fs, Ft       = computeF_pp(freqList, χpp_s, χpp_s, χ0_pp_full)
Γs, Γt       = computeΓ_pp(freqList, χpp_s, χpp_t, χ0_pp_full, shift, nBose, nFermi)
# Φs = reshape(Fs,2*nFermi,2*nFermi,2*nBose+1)[26:75,26:75,51] .- Γs
# Φt = reshape(Ft,2*nFermi,2*nFermi,2*nBose+1)[26:75,26:75,51] .- Γt


# This segment computes quntities in the ph channel
# We first subtract the unconnected part of the susceptibility
#
χ0_full = lDGAPostprocessing.computeχ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp.parent, β)


χm_gen = χph_upup .- χph_updo
χd_gen = χph_upup .+ χph_updo
Γm, Γd = -1.0 .* computeΓ_ph(freqList, χm_gen, χd_gen, χ0_full,nBose,nFermi)
#Γm2, Γd2 = -1.0 .* lDGAPostprocessing.computeΓ_ph2(freqList, χph_upup .+ χph_updo, χph_upup .- χph_updo, χ0_full,nBose,nFermi)
Fm, Fd = computeF_ph(freqList, χph_upup, χph_updo, χ0_full)
Φm = reshape(Fm,2*nFermi,2*nFermi,2*nBose+1)[26:75,26:75,51] .- Γm[26:75,26:75,51]
Φd = reshape(Fd,2*nFermi,2*nFermi,2*nBose+1)[26:75,26:75,51] .- Γd[26:75,26:75,51]

res = isfile(dataPath * "/chi_asympt") ? read_chi_asympt(dataPath * "/chi_asympt") : error("chi_asympt not found!")
χ_d_asympt, χ_m_asympt, χ_pp_asympt = res
