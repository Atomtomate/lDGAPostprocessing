@inline get_symm_f(f::Array{Complex{Float64},1}, i::Int64) = (i < 0) ? conj(f[-i]) : f[i+1]


function identityMap(nB_min, nB_max, nF_min, nF_max, shift, offset, base)
    freqList = [(i,j,k) for i in (nB_min:nB_max) for j in (nF_min:nF_max) .- trunc(Int64,shift*i/2) for k in (nF_min:nF_max) .- trunc(Int64,shift*i/2)]
    fullMap = Dict(tripleToInt(freqList[i][1],freqList[i][2],freqList[i][3],offset,base,base*base) => UInt32(i) for i in 1:length(freqList))
    return fullMap
end


# ==================== GF Stuff ====================
function subtract_ω0!(freqList::Array, arr::Array{T}, gImp::Array{Complex{Float64}, 1}, β::Float64) where T <: Number
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

function computeχ!(freq_list::AbstractArray{Int,2}, F_sp::Array{Complex{Float64}}, χ0::Dict{Tuple{Int,Int},Complex{Float64}})
end

function computeχ(freqList::Array, F::Array{T,1}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int64, nFermi::Int64) where T
    res = Array{T}(undef,2*nBose+1, 2*nFermi, 2*nFermi)
    for (ωn,ω) in enumerate(-nBose:nBose)
        freqSegment = (ωn-1)*(2*nFermi)*(2*nFermi)+1:(ωn+0)*(2*nFermi)*(2*nFermi)
        res[ωn,:,:] = -1.0 * inv(reshape(F[freqSegment],2*nFermi,2*nFermi))
        fermi_grid = freqList[freqSegment]
        for (νn,fg) in enumerate(fermi_grid[1:2*nFermi])
            res[ωn,νn,νn] += 1.0/χ0[(ω,fg[3])]
        end
    end
    return res
end


function write_fort_dir(prefix::String, freqList::Array, arr_ch::Array{Complex{Float64},3}, arr_sp::Array{Complex{Float64},3}, dirname::String, nBose::Int, nFermi::Int)
    if isdir(dirname)
        println("ERROR: Directory already exists. Skipping output")
        return
    else
        mkdir(dirname)
    end
    for ωn in 1:size(arr_ch,1)
        freqSegment = (ωn-1)*(2*nFermi)*(2*nFermi)+1:(ωn+0)*(2*nFermi)*(2*nFermi)
        freq_sub_grid = freqList[freqSegment]
        open(dirname * "/" * prefix * lpad(ωn-1, 3, "0"), "w") do f
            for i in 1:2*nFermi
                for j in 1:2*nFermi
                    @printf(f, "  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f\n", 
                        float(freq_sub_grid[i][1]), float(freq_sub_grid[(i-1)*(2*nFermi)+j][2]), 
                        float(freq_sub_grid[(i-1)*(2*nFermi)+j][3]),
                        real(arr_ch[ωn, i, j]), imag(arr_ch[ωn, i, j]),
                        real(arr_sp[ωn, i, j]), imag(arr_sp[ωn, i, j]))
                end
            end
        end
    end
end

function write_fort_dir(prefix::String, freqList::Array, arr_ch::Array{Complex{Float64},1}, arr_sp::Array{Complex{Float64},1}, dirname::String, nBose::Int, nFermi::Int)
    if isdir(dirname)
        println("ERROR: Directory already exists. Skipping output")
        return
    else
        mkdir(dirname)
    end
    nF2 = (2*nFermi*2*nFermi)
    for ωn in 1:(2*nBose+1)
        freqSegment = (ωn-1)*nF2+1:(ωn+0)*nF2
        freq_sub_grid = freqList[freqSegment]
        open(dirname * "/" * prefix * lpad(ωn-1, 3, "0"), "w") do f
            for i in 1:nF2
                ind = (ωn-1)*nF2+i
                @printf(f, "  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f\n", 
                        float(freq_sub_grid[i][1]), float(freq_sub_grid[i][2]), float(freq_sub_grid[i][3]),
                        real(arr_ch[ind]), imag(arr_ch[ind]),
                        real(arr_sp[ind]), imag(arr_sp[ind]))
            end
        end
    end
end
