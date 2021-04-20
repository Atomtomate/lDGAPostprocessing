function read_vert_chi(fname)
    inp = readdlm(fname, Float64)
    freq = Float64.(inp[:,1:3])
    Fup = inp[:,4] + inp[:,5]*1im
    Fdo = inp[:,6] + inp[:,7]*1im
    return freq, Fup, Fdo
end

function read_gm_wim(nFreq, filename; storedInverse, storeFull=false)
    GFString = open(filename, "r") do f
        readlines(f)
    end

    if size(GFString, 1)*(1 + 1*storeFull) < nFreq
        throw(BoundsError("nFermFreq in simulation parameters too large!"))
    end
    
    tmp = parse.(Float64,hcat(split.(GFString)...)[2:end,:]) # Construct a 2xN array of floats (re,im as 1st index)
    tmp = tmp[1,:] .+ tmp[2,:].*1im

    if storedInverse
        tmp = 1 ./ tmp
    end
    
    GF = Array{Complex{Float64}}(undef, nFreq)
    if storeFull
        NH = Int(nFreq/2)
        GF[1:(NH-1)] = reverse(conj.(tmp[2:NH]))
        GF[NH:nFreq] = tmp[1:NH]
    else
        GF = tmp[1:nFreq]
    end
    return GF
end
