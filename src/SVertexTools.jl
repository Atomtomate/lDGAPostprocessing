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

function expand_test(Fup_red, Fdo_red, freqList_map, freqList, parents, ops, nBose, nFermi)
    off(f) = f[1]+nBose+1,f[2]+nFermi+1,f[3]+nFermi+1
    Fup_full = -1 .* ones(Complex{Float64}, length(freqList))
    Fdo_full = -1 .* ones(Complex{Float64}, length(freqList))
    for (k,v) in freqList_map
        Fup_full[k] = Fup_red[v]
        Fdo_full[k] = Fdo_red[v]
    end
    return Fup_full, Fdo_full
end
function expand(Fup_red, Fdo_red, transform, freqList_map, freqList, parents, ops, nBose, nFermi)
    off(f) = f[1]+nBose+1,f[2]+nFermi+1,f[3]+nFermi+1
    Fup_full = Array{eltype(Fup_red)}(undef, length(freqList))
    Fdo_full = Array{eltype(Fdo_red)}(undef, length(freqList))
    Fup_done = falses(2*nBose,length(freqList))
    Fdo_done = falses(2*nBose,length(freqList))
    for (k,v) in freqList_map
        println("--- setting $k in full, $v in red")
        Fup_full[k] = Fup_red[v]
        Fdo_full[k] = Fdo_red[v]
        Fup_done[k] = true
        Fdo_done[k] = true
    end
    open = Stack{eltype(parents)}()
    for i in 1:length(freqList)
        next = i
        ωn,νn,νpn  = off(freqList[next])
        println("=== in $i, parents[$i]=$next. status $(Fup_done[i])")
        while !Fup_done[next]
            println(" not done, pushing parent $(parents[next])")
            push!(open, next)
            next = parents[next]
        end
        while length(open) > 0 
            prev = next
            next = pop!(open)
            println(" copying $prev to $next, val = $(Fup_full[prev]), transformed to $(transform(Fup_full[prev], ops[prev])), using op = $(ops[prev]), test: $(ops[next])")
            Fup_done[next] = true
            Fup_full[next] = transform(Fup_full[prev], ops[next])
            Fdo_full[next] = transform(Fdo_full[prev], ops[next])
        end
    end
    return Fup_full, Fdo_full
end
