@inline get_symm_f(f::Array{Complex{Float64},1}, i::Int64) = (i < 0) ? conj(f[-i]) : f[i+1]


function identityMap(nB_min, nB_max, nF_min, nF_max, shift, offset, base)
    freqList = [(i,j,k) for i in (nB_min:nB_max) for j in (nF_min:nF_max) .- trunc(Int64,shift*i/2) for k in (nF_min:nF_max) .- trunc(Int64,shift*i/2)]
    fullMap = Dict(tripleToInt(freqList[i][1],freqList[i][2],freqList[i][3],offset,base,base*base) => UInt32(i) for i in 1:length(freqList))
    return fullMap
end

function cut_margin(freqList, arr, Ncut_ω::Int, Ncut_ν::Int, Nν_max::Int, Nω_max::Int, shift::Int)
    freqList_cut = similar(freqList, 0)
    arr_cut      = similar(arr, 0)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        # i < 10 && println(abs(ω), " ?>? ", Nω_max-Ncut_ω,", ", abs(ν+shift*ω/2), " ?>? ", Nν_max-Ncut_ν, ", ", abs(νp+shift*ω/2), " ?>? ", Nν_max-Ncut_ν)
        if abs(ω) <= Nω_max-Ncut_ω && ν+shift*ω/2 >=  -(Nν_max-Ncut_ν) && ν+shift*ω/2 < (Nν_max-Ncut_ν) &&  νp+shift*ω/2 >=  -(Nν_max-Ncut_ν) && νp+shift*ω/2 < (Nν_max-Ncut_ν) 
            push!(freqList_cut, freqList[i])
            push!(arr_cut, arr[i])
            # i < 100 && print("pushing: ",freqList[i])
        end
    end
    return freqList_cut, arr_cut
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

function add_χ₀_ω₀!(freqList::Array, arr::Array{T}, gImp::Array{Complex{Float64}, 1}, β::Float64) where T <: Number
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        if ω == 0
            arr[i] -= β*get_symm_f(gImp, ν)*get_symm_f(gImp, νp)
        end
    end
end


function computeχ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, gImp::Array{Complex{Float64}, 1}, β::Float64; mode=:ph)
    !(mode in [:ph, :pp]) && error("unkown mode")
    χ0 = Dict{Tuple{Int,Int},Complex{Float64}}()
    for ω in ω_range, ν in ν_range
        χ0[(ω,ν)] = (mode == :ph) ? -β*get_symm_f(gImp, ν)*get_symm_f(gImp, ν+ω) : -β*get_symm_f(gImp, ν)*get_symm_f(gImp, ν-ω)
    end
    return χ0
end

function computeF_ph(freqList::AbstractArray{Int,2}, χ_upup::SVertex{T}, χ_updo::SVertex{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}) where T
    F_den = similar(χ_upup)
    F_mag = similar(χ_upup)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i,:]
        sub = ν == νp ? χ0[(ω,ν)] : 0.0
        F_den[ω,ν,νp] = (-1.0/χ0[(ω,ν)])*(χ_upup[ω,ν,νp]+χ_updo[ω,ν,νp]-sub)*(1.0/χ0[(ω,νp)])
        F_mag[ω,ν,νp] = (-1.0/χ0[(ω,ν)])*(χ_upup[ω,ν,νp]-χ_updo[ω,ν,νp]-sub)*(1.0/χ0[(ω,νp)])
    end
    return F_den, F_mag
end

function computeF_pp(freqList::AbstractArray{Int,2}, χ_s::SVertex{T}, χ_t::SVertex{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}) where T
    F_s = similar(χ_s)
    F_t = similar(χ_t)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i,:]
        sub = ν == νp ? χ0[(ω,ν)] : 0.0
        F_s[ω,ν,νp] = (-1.0/χ0[(ω,ν)])*(χ_s[ω,ν,νp]-sub)*(1.0/χ0[(ω,νp)])
        F_t[ω,ν,νp] = (-1.0/χ0[(ω,ν)])*(χ_t[ω,ν,νp]-sub)*(1.0/χ0[(ω,νp)])
    end
    return F_s, F_t
end

function computeΓ_ph(freqList::Array, χm::Array{T,1}, χd::Array{T,1}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int64, nFermi::Int64) where T
    Γm = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    Γd = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    for (ωn,ω) in enumerate(-nBose:nBose)
        freqSegment = findall(fl -> fl[1] == ω, freqList)
        Γm[ωn,:,:] = inv(transpose(reshape(χm[freqSegment],2*nFermi,2*nFermi)))
        Γd[ωn,:,:] = inv(transpose(reshape(χd[freqSegment],2*nFermi,2*nFermi)))
        fermi_grid = freqList[freqSegment]
        for (νn,fg) in enumerate(fermi_grid[1:2*nFermi])
            Γm[ωn,νn,νn] -= 1.0/χ0[(ω,fg[3])]
            Γd[ωn,νn,νn] -= 1.0/χ0[(ω,fg[3])]
        end
    end
    return Γm, Γd
end

function computeΓ_pp(freqList::Array, χs::Array{T,1}, χt::Array{T,1}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int64, nFermi::Int64) where T
    Γs = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    Γt = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    for (ωn,ω) in enumerate(-nBose:nBose)
        freqSegment = findall(fl -> fl[1] == ω, freqList)
        fermi_grid = freqList[freqSegment]
        Γs[ωn,:,:] = transpose(reshape(χs[freqSegment],2*nFermi,2*nFermi))#
        Γt[ωn,:,:] = transpose(reshape(χt[freqSegment],2*nFermi,2*nFermi))#
        for (νn,fg) in enumerate(fermi_grid[1:2*nFermi])
            Γs[ωn,νn,νn] -= χ0[(ω,fg[3])]
            Γt[ωn,νn,νn] += χ0[(ω,fg[3])]
        end
        Γs[ωn,:,:] = 4 .* inv(Γs[ωn,:,:])
        Γt[ωn,:,:] = 4 .* inv(Γt[ωn,:,:])
        for (νn,fg) in enumerate(fermi_grid[1:2*nFermi])
            Γs[ωn,νn,νn] += 2.0/χ0[(ω,fg[3])]
            Γt[ωn,νn,νn] -= 2.0/χ0[(ω,fg[3])]
        end
    end
    return Γs, Γt
end

function χph_to_χpp(freqList::AbstractArray{Int,2}, χph_upup::SVertex{T}, χph_updo::SVertex{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}) where T
    χpp_s = deepcopy(χph_upup) 
    χpp_t = deepcopy(χph_upup) 
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i,:]
        χpp_s = - χ0[(ω,ν)] - χph_upup[ω-ν-νp,ν,νp] + 2*χph_updo[ω-ν-νp,ν,νp]
        χpp_t = + χ0[(ω,ν)] + χph_upup[ω-ν-νp,ν,νp]
    end
    # χpp = similar(χph)
    # chi_inv_singlet(l,s)=-chi_0(i,j,k)-chi_up_up(i-j-k-1,j,k)+2.0d0*chi_up_down(i-j-k-1,j,k)          !if we have chi_ph, i.e. we have to make the frequency shift
    # chi_inv_triplet(l,s)= chi_0(i,j,k)+chi_up_up(i-j-k-1,j,k)                                                            ! to get chi in the pp-notation
    χpp_s, χpp_t
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
