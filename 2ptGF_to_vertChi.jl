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

    #transform_f helps construct the full two-partilce Green's funciton based which symmetry was used (as indicated by ops)
    #Note that here F refers to G2 and Fup, Fdo to G2upup, G2updown respectively. In this script the full two-particle vertex F never appears explicitly.
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

function write_GAMMA_DM_FULLRANGE_dbg(Γm::Array{ComplexF64,3}, Γd::Array{ComplexF64,3}, fname::String, β::Float64, nBose::Int, nFermi::Int, shift::Union{Int,Bool})
    open(fname, "w") do f
        for (ωi,ωn) in enumerate(-nBose:nBose)
            νgrid =  [i - trunc(Int64,shift*ωn/2) for i in (-nFermi:nFermi-1)]
            for (νi,νn) in enumerate(νgrid)
                for (νpi, νpn) in enumerate(νgrid)
                    lDGAPostprocessing.@printf(f, "%17.10f%17.10f%17.10f%17.10f%17.10f%17.10f%17.10f\n", 
                                               ωn,  νn,  νpn,  
                        real(Γd[ωi,νi,νpi]), imag(Γd[ωi,νi,νpi]),
                        real(Γm[ωi,νi,νpi]), imag(Γm[ωi,νi,νpi]))
                end
            end
        end
    end
end

function write_vert_chi_dbg(ver::Array{ComplexF64,3}, verdo::Array{ComplexF64,3}, fname::String, β::Float64, nBose::Int, nFermi::Int, shift::Union{Int,Bool})
    open(fname, "w") do f
        for (ωi,ωn) in enumerate(-nBose:nBose)
            νgrid =  [i - trunc(Int64,shift*ωn/2) for i in (-nFermi:nFermi-1)]
            for (νi,νn) in enumerate(νgrid)
                for (νpi, νpn) in enumerate(νgrid)
                    lDGAPostprocessing.@printf(f, "  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f\n", 
                                               2*ωn*π/β,  (2*νn+1)*π/β,  (2*νpn+1)*π/β,  
                        real(ver[ωi,νi,νpi]), imag(ver[ωi,νi,νpi]),
                        real(verdo[ωi,νi,νpi]), imag(verdo[ωi,νi,νpi]))
                end
            end
        end
    end
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

χDMFTch = permutedims(reshape((χ_upup .+ χ_updo), 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTsp = permutedims(reshape((χ_upup .- χ_updo), 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTuu = permutedims(reshape(χ_upup, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
χDMFTud = permutedims(reshape(χ_updo, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
Γsp, Γch = computeΓ_ph(χDMFTsp,  χDMFTch, GImp, β,nBose,nFermi,shift)

println("Writing vert_chi!")
write_vert_chi_dbg(χDMFTuu, χDMFTud, joinpath(dataPath,"vert_chi"), β, nBose, nFermi, shift)
println("Writing GAMMA_DM_FULLRANGE")
write_GAMMA_DM_FULLRANGE_dbg(Γsp, Γch, joinpath(dataPath,"GAMMA_DM_FULLRANGE"), β, nBose, nFermi, shift)
