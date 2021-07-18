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
        throw(ArgumentError("nFermFreq in simulation parameters too large! Got $(size(GFString, 1)) lines but need data for $(nFreq) frequencies."))
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

function read_hubb_dat(file)
    content = open(file) do f
        readlines(f)
    end
    U = parse(Float64,split(replace(content[2], "d0" => ""), r"(\s|,)+")[1])
    β = parse(Float64,split(replace(content[4], "d0" => ""), r"(\s|,)+")[1])
    n = parse(Float64,split(replace(content[8], "d0" => ""), r"(\s|,)+")[2])
    return U, β, n
end


function read_anderson_parameters(file)
    content = open(file) do f
        readlines(f)
    end
    
    in_epsk = false
    in_tpar = false
    ϵₖ = []
    Vₖ = []
    μ = 0
    for line in content
        if "Eps(k)" == strip(line)
            in_epsk = true
            continue
        elseif "tpar(k)" == strip(line)
            in_epsk = false
            in_tpar = true
            continue
        end
        
        if in_epsk
            push!(ϵₖ, parse(Float64, line))
        elseif in_tpar
            # skip last line, which is mu
            if length(Vₖ) < length(ϵₖ)
                push!(Vₖ, parse(Float64, line))
            else
                if occursin("#", line)
                    μ = parse(Float64, line[1:(findfirst("#", line))[1] - 1])
                else
                    μ = parse(Float64, line)
                end
            end
        end
    end
    return convert(Array{Float64,1}, ϵₖ), convert(Array{Float64,1}, Vₖ), μ
end
