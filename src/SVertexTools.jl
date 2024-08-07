indices(sv::SVertex) = sort([intToTriple(Int, el, sv.offset, sv.nB) for el in keys(sv.map)])

function full(sv::SVertex{T}, indices) where T
    res = Array{T}(undef, sv.nBose, sv.nFermi, sv.nFermi)
    i = 1
    for nb in 1:sv.nBose
        for nf in 1:sv.nFermi
            for nfp in 1:sv.nFermi
                res[nb,nf,nfp] = sv[indices[i]...]
                i += 1
            end
        end
    end
    return res
end

function full(sv::SVertex{T}) where T
    sv.nBose*sv.nFermi*sv.nFermi != length(sv.map) && throw(ArgumentError("The total number of bosonic and fermionic frequencies must match the length of the map."))
    ind = indice(sv)
    return ind, full(sv, ind)
end


function full!(res::Array{T,3}, sv::SVertex{T}, indices) where T
    i = 1
    for nb in 1:sv.nBose
        for nf in 1:sv.nFermi
            for nfp in 1:sv.nFermi
                res[nb,nf,nfp] = sv[indices[i]...]
                i += 1
            end
        end
    end
end

function expand_test(TwoPartGF_upup, TwoPartGF_updo, freqList_map, freqList, parents, ops, nBose, nFermi)
    off(f) = f[1]+nBose+1,f[2]+nFermi+1,f[3]+nFermi+1
    Fup_full = -1 .* ones(Complex{Float64}, length(freqList))
    Fdo_full = -1 .* ones(Complex{Float64}, length(freqList))
    for (k,v) in freqList_map
        Fup_full[k] = TwoPartGF_upup[v]
        Fdo_full[k] = TwoPartGF_updo[v]
    end
    return Fup_full, Fdo_full
end

function expand(TwoPartGF_upup, TwoPartGF_updo, transform!, freqList_map, freqList, parents, ops, nBose, nFermi)
    #off(f) = f[1]+nBose+1,f[2]+nFermi+1,f[3]+nFermi+1
    Fup_full = Array{eltype(TwoPartGF_upup)}(undef, length(freqList))
    Fdo_full = Array{eltype(TwoPartGF_updo)}(undef, length(freqList))
    done = falses(length(freqList))
    #There are a few frequency tuples in the full grid which directly map to the reduced grid. Here we expand these first.
    for (k,v) in freqList_map
        #println("--- setting $k in full, $v in red")
        Fup_full[k] = TwoPartGF_upup[v]
        Fdo_full[k] = TwoPartGF_updo[v]
        done[k] = true
    end
    open = Stack{eltype(parents)}()
    #Now we expand the frequencies which map onto frequencies which are not directly on the reduced grid. We follow a chain of mappings on the full grid recursevly until we find a frequency tuple on the reduced grid: i.e: freq1_full -> freq2_full -> freq3_full -> freq1_reduced.
    for i in 1:length(freqList)
        next = i
        while !done[next]
            push!(open, next)
            next = parents[next]
        end
        while length(open) > 0 
            prev = next
            next = pop!(open)
            done[next] = true
            #Which symmetry is considered (complex conjugate, *(-1), ...) is given by the transform function, with "ops" indicating which symmetry is used if multiple are applicable.
            transform!(Fup_full, Fdo_full, prev, next, ops)
        end
    end
    return Fup_full, Fdo_full
end
