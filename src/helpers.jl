@inline get_symm_f(f::Array{Complex{Float64},1}, i::Int64) = (i < 0) ? conj(f[-i]) : f[i+1]


function identityMap(nB_min, nB_max, nF_min, nF_max, shift, offset, base)
    freqList = [(i,j,k) for i in (nB_min:nB_max) for j in (nF_min:nF_max) .- trunc(Int64,shift*i/2) for k in (nF_min:nF_max) .- trunc(Int64,shift*i/2)]
    fullMap = Dict(tripleToInt(freqList[i][1],freqList[i][2],freqList[i][3],offset,base,base*base) => UInt32(i) for i in 1:length(freqList))
    return fullMap
end


# ==================== GF Stuff ====================
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

function computeχ!(freq_list::AbstractArray{Int,2}, F_sp::Array{Complex{Float64}}, χ0::Dict{Tuple{Int,Int},Complex{Float64}})
end

function computeχ(F::SVertex{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}) where T
    res = inv(F)
    for ωn in 1:size(F,1)
        for νn in 1:size(F,3)
            res[ωn,νn,νn] -= 1.0/χ0[(ωn,νn)]
        end
    end
    #res[ωn,νn,νn] -= 1.0/χ[ωn,νn]
    return res
end

function computeχ(F::Array{T,1}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int64, nFermi::Int64) where T
    res = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    for (ωn,ω) in enumerate(-nBose:nBose)
        res[ωn,:,:] = inv(reshape(F[(ωn-1)*(2*nFermi)*(2*nFermi)+1:(ωn+0)*(2*nFermi)*(2*nFermi)],2*nFermi,2*nFermi))
        for (νn,ν) in enumerate(-2:2-1)
            res[ωn,νn,νn] -= 1.0/χ0[(ω,ν)]
        end
    end
    #res[ωn,νn,νn] -= 1.0/χ[ωn,νn]
    return res
end

