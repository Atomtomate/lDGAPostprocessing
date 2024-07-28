using Pkg
Pkg.activate(@__DIR__)
using lDGAPostprocessing, jED

using EquivalenceClassesConstructor
using Test
using DataStructures
using DelimitedFiles, JLD2
using TOML
include("./vertexIntTriple.jl")

function expand_TwoPartGF(freqFile, dataPath, read_vert_chi_flag)
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
    fname = read_vert_chi_flag ? "vert_chi" : "2_part_gf_red"
    freq_red, TwoPartGF_upup_red, TwoPartGF_updo_red = read_vert_chi(joinpath(dataPath,fname)) 
    χ_upup, χ_updo = expand(TwoPartGF_upup_red, TwoPartGF_updo_red, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);
    if !read_vert_chi_flag
        lDGAPostprocessing.add_χ₀_ω₀!(freqList, χ_upup, GImp.parent, β)
        lDGAPostprocessing.add_χ₀_ω₀!(freqList, χ_updo, GImp.parent, β)
    end
    return χ_upup, χ_updo
end

gridPath = ARGS[1]    # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
dataPath = ARGS[2]    # "/scratch/usr/hhpstobb/lDGA/tests/non_qubic/ed_vertex"
read_vert_chi_flag  = length(ARGS) == 3 ? parse(Bool,lowercase(ARGS[3])) : false

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


χph_upup, χph_updo = expand_TwoPartGF(gridPath*"/freqList.jld2", dataPath, read_vert_chi_flag)
lDGAPostprocessing.add_χ₀_ω₀!(freqList, χph_upup, GImp.parent, β)
lDGAPostprocessing.add_χ₀_ω₀!(freqList, χph_updo, GImp.parent, β)

# lDGAPostprocessing.write_vert_chi(freqList, χph_upup, χph_updo, ".", nBose, nFermi)
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
Fs, Ft       = computeF_pp(freqList, χpp_s[:], χpp_t[:], χ0_pp_full, nBose, nFermi)
Γs, Γt       = computeΓ_pp(freqList, χpp_s, χpp_t, χ0_pp_full, shift, nBose, nFermi)
Φs = reshape(Fs,2*nFermi,2*nFermi,2*nBose+1) .- Γs
Φt = reshape(Ft,2*nFermi,2*nFermi,2*nBose+1) .- Γt
println("Done with pp channel!")

χDMFTch = permutedims(reshape((χph_upup .+ χph_updo), 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTsp = permutedims(reshape((χph_upup .- χph_updo), 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTuu = permutedims(reshape(χph_upup, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTud = permutedims(reshape(χph_updo, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
Γsp, Γch = computeΓ_ph(χDMFTsp,  χDMFTch, GImp, β,nBose,nFermi,shift)

# This segment computes quntities in the ph channel
# We first subtract the unconnected part of the susceptibility
#
χ0_full = lDGAPostprocessing.computeχ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp.parent, β)
Fm, Fd = computeF_ph(freqList, χph_upup, χph_updo, χ0_full, nBose, nFermi)
println("Done with ph channel!")

res = isfile(dataPath * "/chi_asympt") ? read_chi_asympt(dataPath * "/chi_asympt") : (println("WARNING! chi_asympt not found!"); [nothing, nothing, nothing])
χ_d_asympt, χ_m_asympt, χ_pp_asympt = res

E_kin_DMFT, E_pot_DMFT  = calc_E_ED(νnGrid[0:last(axes(GImp,1))], ϵₖ, Vₖ, GImp.parent, nden, U, β, μ)

jldopen(dataPath*"/DMFT_out.jld2", "w") do f
    f["Γch"] = -1.0 .* Γch
    f["Γsp"] = -1.0 .* Γsp
    f["Φpp_s"] = permutedims(Φs, [3,1,2])
    f["Φpp_t"] = permutedims(Φt, [3,1,2])
    f["χDMFTch"] = χDMFTch
    f["χDMFTsp"] = χDMFTsp
    f["χ_ch_asympt"] = χ_d_asympt
    f["χ_sp_asympt"] = χ_m_asympt
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
    f["kGrid"] = cfg["parameters"]["lattice"]
end
