@inline get_symm_f(f::Array{Complex{Float64},1}, i::Int64) = (i < 0) ? conj(f[-i]) : f[i+1]


function identityMap(nB_min, nB_max, nF_min, nF_max, shift, offset, base)
    freqList = [(i,j,k) for i in (nB_min:nB_max) for j in (nF_min:nF_max) .- trunc(Int64,shift*i/2) for k in (nF_min:nF_max) .- trunc(Int64,shift*i/2)]
    fullMap = Dict(tripleToInt(freqList[i][1],freqList[i][2],freqList[i][3],offset,base,base*base) => UInt32(i) for i in 1:length(freqList))
    return fullMap
end


# ==================== GF Stuff ====================
function FUpDo_from_χDMFT(χupdo, GImp, freqFile, β)
    FUpDo = similar(χupdo)
    jldopen(freqFile, "r") do freqFile
        freqList = freqFile["freqList"]
        n_iω = freqFile["nBose"]
        n_iν = freqFile["nFermi"]
        shift = freqFile["shift"]
        for f in freqList
            i = f[1] + n_iω+1
            j = f[2] + n_iν+1 + trunc(Int, shift*f[1]/2)
            k = f[3] + n_iν+1 + trunc(Int, shift*f[1]/2)
            FUpDo[i,j,k] = χupdo[i,j,k]/(β^2 * get_symm_f(GImp,f[2]) * get_symm_f(GImp,f[1]+f[2])
                               * get_symm_f(GImp,f[3]) * get_symm_f(GImp,f[1]+f[3]))
        end
    end
    return FUpDo
end

function add_χ₀!(freqList::Array, arr::Array{T}, gImp::Array{Complex{Float64}, 1}, β::Float64) where T <: Number
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        if ω == 0
            arr[i] -= β*get_symm_f(gImp, ν)*get_symm_f(gImp, νp)
        end
    end
end


function computeχ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, gImp::Array{Complex{Float64}, 1}, β::Float64)
    χ0 = Dict{Tuple{Int,Int},Complex{Float64}}()
    for ω in ω_range, ν in ν_range
        χ0[(ω,ν)] = -β*get_symm_f(gImp, ν)*get_symm_f(gImp, ν+ω)
    end
    return χ0
end

function computeF(freqList::AbstractArray{Int,2}, F_up::SVertex{T}, F_do::SVertex{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}) where T
    F_den = similar(F_up)
    F_mag = similar(F_up)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i,:]
        sub = ν == νp ? χ0[(ω,ν)] : 0.0
        F_den[ω,ν,νp] = (-1.0/χ0[(ω,ν)])*(F_up[ω,ν,νp]+F_do[ω,ν,νp]-sub)*(1.0/χ0[(ω,νp)])
        F_mag[ω,ν,νp] = (-1.0/χ0[(ω,ν)])*(F_up[ω,ν,νp]-F_do[ω,ν,νp]-sub)*(1.0/χ0[(ω,νp)])
    end
    return F_den, F_mag
end

function computeΓ(freqList::Array, χ::Array{T,1}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int64, nFermi::Int64) where T
    res = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    for (ωn,ω) in enumerate(-nBose:nBose)
        freqSegment = (ωn-1)*(2*nFermi)*(2*nFermi)+1:(ωn+0)*(2*nFermi)*(2*nFermi)
        res[ωn,:,:] = inv(reshape(χ[freqSegment],2*nFermi,2*nFermi))
        fermi_grid = freqList[freqSegment]
        for (νn,fg) in enumerate(fermi_grid[1:2*nFermi])
            res[ωn,νn,νn] -= 1.0/χ0[(ω,fg[3])]
        end
    end
    return res
end

function calc_E_ED(iνₙ, ϵₖ, Vₖ, GImp, n, U, β, μ; full=false)
    E_kin = 0.0
    E_pot = 0.0
    vk = sum(Vₖ .^ 2)
    Σ_hartree = n * U/2
    E_pot_tail = (U^2)/2 * n * (1-n/2) - Σ_hartree*(Σ_hartree-μ)
    E_kin_tail = vk

    for n in 1:length(GImp)
        Δ_n = sum((Vₖ .^ 2) ./ (iνₙ[n] .- ϵₖ))
        Σ_n = iνₙ[n] .- Δ_n .- 1.0 ./ GImp[n] .+ μ
        E_kin += 2*real(GImp[n] * Δ_n - E_kin_tail/(iνₙ[n]^2))
        E_pot += 2*real(GImp[n] * Σ_n - E_pot_tail/(iνₙ[n]^2))
    end
    E_kin = E_kin .* (2/β) - (β/2) .* E_kin_tail
    E_pot = E_pot .* (1/β) .+ 0.5*Σ_hartree .- (β/4) .* E_pot_tail
    return E_kin, E_pot
end
