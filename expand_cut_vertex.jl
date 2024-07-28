using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/lDGAPostprocessing")
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
    TwoPartGF_upup, TwoPartGF_updo = expand(TwoPartGF_upup_red, TwoPartGF_updo_red, transform_f, freqRed_map, freqList, parents, ops, nBose, nFermi);
    return TwoPartGF_upup, TwoPartGF_updo
end

function cut_margin(freqList, arr, Ncut_ω::Int, Ncut_ν::Int, Nν_max::Int, Nω_max::Int, shift::Int)
    freqList_cut = similar(freqList, 0)
    arr_cut      = similar(arr, 0)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        si = shift*trunc(Int,ω/2)
        if abs(ω) <= Nω_max-Ncut_ω && ν >=  -(Nν_max-Ncut_ν)-si && ν < (Nν_max-Ncut_ν)-si &&  νp >=  -(Nν_max-Ncut_ν)-si && νp < (Nν_max-Ncut_ν)-si
            push!(freqList_cut, freqList[i])
            push!(arr_cut, arr[i])
        end
    end
    return freqList_cut, arr_cut
end

gridPath = ARGS[1] # "/scratch/usr/hhpstobb/grids/b10_f20_s1"
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

# ======================================== Unpack 2-Part-GF ========================================
TwoPartGF_upup, TwoPartGF_updo = expand_TwoPartGF(gridPath*"/freqList.jld2", dataPath)
println("Done expanding!")
flush(stderr)
flush(stdout)

for Ncut in 0:10:130
    nFermi_N = nFermi-Ncut
    nBose_N  = nBose -Ncut


    freqList_N, TwoPartGF_updo_N = cut_margin(freqList, TwoPartGF_updo, Ncut, Ncut, nFermi, nBose, shift)
    freqList_N, TwoPartGF_upup_N = cut_margin(freqList, TwoPartGF_upup, Ncut, Ncut, nFermi, nBose, shift)

    println("cut from length ", length(freqList), " to ", length(freqList_N))
    kG = jED.gen_kGrid(cfg["parameters"]["lattice"], 100)
     
    p  = AIMParams(ϵₖ, Vₖ)
    Nν = 4*(nBose_N+nFermi_N+1)

    basis  = jED.Basis(length(ϵₖ) + 1);
    νnGrid = jED.OffsetVector([1im * (2*n+1)*π/β for n in 0:Nν-1], 0:Nν-1)
    G0W    = GWeiss(νnGrid, μ, p)

    model  = AIM(ϵₖ, Vₖ, μ, U)
    es     = Eigenspace(model, basis);
    GImp, nden = calc_GF_1(basis, es, νnGrid, β, with_density=true)
    ΣImp   = Σ_from_GImp(G0W, GImp)
    gLoc   = GLoc(ΣImp, μ, νnGrid, kG)

    println("Expanding Vertex for nFermi=",nFermi_N,", nBose=",nBose_N,", shift=",shift)
    flush(stderr)
    flush(stdout)

    χ0_full = lDGAPostprocessing.computeχ0(-nBose_N:nBose_N, -(nFermi_N+2*nBose_N):(nFermi_N+2*nBose_N)-1, GImp.parent, β)

    lDGAPostprocessing.add_χ₀_ω₀!(freqList_N, TwoPartGF_upup_N, GImp.parent, β)
    lDGAPostprocessing.add_χ₀_ω₀!(freqList_N, TwoPartGF_updo_N, GImp.parent, β)
    Γch = -1.0 .* lDGAPostprocessing.computeΓ(freqList_N, TwoPartGF_upup_N .+ TwoPartGF_updo_N, χ0_full,nBose_N,nFermi_N)
    Γsp = -1.0 .* lDGAPostprocessing.computeΓ(freqList_N, TwoPartGF_upup_N .- TwoPartGF_updo_N, χ0_full,nBose_N,nFermi_N)
    println("Done calculating vertex!")
    flush(stderr)
    flush(stdout)

    E_kin_DMFT, E_pot_DMFT  = calc_E_ED(νnGrid[0:last(axes(GImp,1))], ϵₖ, Vₖ, GImp.parent, nden, U, β, μ)
    res = isfile(dataPath * "/chi_asympt") ? read_chi_asympt(dataPath * "/chi_asympt") : error("chi_asympt not found!")
    χ_ch_asympt, χ_sp_asympt, χ_pp_asympt = res

    jldopen(dataPath*"/ED_out_Ncut$Ncut.jld2", "w") do f
        f["Γch"] = Γch
        f["Γsp"] = Γsp
        f["χDMFTch"] = permutedims(reshape(TwoPartGF_upup_N .+ TwoPartGF_updo_N, 2*nFermi_N, 2*nFermi_N, 2*nBose_N+1),[3,2,1])
        f["χDMFTsp"] = permutedims(reshape(TwoPartGF_upup_N .- TwoPartGF_updo_N, 2*nFermi_N, 2*nFermi_N, 2*nBose_N+1),[3,2,1])
        f["χ_ch_asympt"] = χ_ch_asympt
        f["χ_sp_asympt"] = χ_sp_asympt
        f["χ_pp_asympt"] = χ_pp_asympt
        f["gImp"] = GImp.parent
        f["g0"] = gLoc.parent
        f["ϵₖ"] = ϵₖ
        f["Vₖ"] = Vₖ
        f["μ"] = μ
        f["U"] = U
        f["β"] = β
        f["nden"] = nden
        f["E_kin_DMFT"] = E_kin_DMFT
        f["E_pot_DMFT"] = E_pot_DMFT
        f["grid_shift"] = shift
        f["grid_nBose"] = nBose_N
        f["grid_nFermi"] = nFermi_N
        f["kGrid"] = cfg["parameters"]["lattice"]
    end
end
