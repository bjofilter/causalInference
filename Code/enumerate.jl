using Graphs
using LinkedLists

mutable struct Helper
    sets::Vector{LinkedList{Int64}}
    pointers::Vector{ListNode{Int64}}
    size::Vector{Int64}
    mcsorder::Vector{Int64}
    lastpush::Vector{Vector{Int64}}
    succ::Vector{Vector{Int64}}
    uptodate::Vector{BigInt}
    maxcard::Int64
    i::Int64
    count::BigInt
end

function setv(G, H, v)
    H.size[v] = -H.size[v]

    H.mcsorder[H.i] = v

    for w in inneighbors(G, v)
        if H.size[w] >= 1
            deleteat!(H.sets[H.size[w]], H.pointers[w])
            H.size[w] += 1
            H.pointers[w] = push!(H.sets[H.size[w]], w)
            push!(H.lastpush[w], v)
        end
    end

    H.maxcard += 1
    while H.maxcard >= 1 && isempty(H.sets[H.maxcard])
        H.maxcard -= 1
    end

    H.i += 1
end

function resetv(G, H, v)
    H.size[v] = -H.size[v]
    H.pointers[v] = push!(H.sets[H.size[v]], v)
    
    for w in inneighbors(G, v)
        if H.size[w] >= 1
            deleteat!(H.sets[H.size[w]], H.pointers[w])
            H.size[w] -= 1
            H.pointers[w] = push!(H.sets[H.size[w]], w)
            pop!(H.lastpush[w])
        end
    end

    H.maxcard = H.size[v]
    H.i -= 1
end

function dfs(H, u, reachable)
    H.uptodate[u] != H.count && return
    for v in H.succ[u]
        if v != u && H.size[u] == H.size[v]
            push!(reachable, v)
            dfs(H, v, reachable)
        end
    end
end

function recenum(G, H)
    if H.i > nv(G)
        # println(H.mcsorder)
        H.count += 1
        return
    end
    
    v = pop!(H.sets[H.maxcard])
    setv(G, H, v)
    recenum(G, H)
    resetv(G, H, v)

    reachable = Vector{Int64}()
    dfs(H, v, reachable)
    
    lastx = v
    for x in reachable
        if H.size[x] == H.size[v]
            deleteat!(H.sets[H.maxcard], H.pointers[x])
            setv(G, H, x)
            recenum(G, H)
            resetv(G, H, x)
        end
        lastx = x
    end

    if length(H.lastpush[lastx]) > 0
        prec = last(H.lastpush[lastx])
        if H.uptodate[prec] != H.count
            H.uptodate[prec] = H.count
            empty!(H.succ[prec])
        end
        push!(H.succ[last(H.lastpush[lastx])], lastx)
    end

    return
end

function enumerate(G)
    n = nv(G)

    sets = [LinkedLists.LinkedList{Int64}() for i = 1:n+1]
    pointers = Vector(undef,n)
    size = Vector{Int64}(undef, n)
    mcsorder = Vector{Int64}(undef, n)
    lastpush = [Vector{Int64}() for i = 1:n]
    succ = [Vector{Int64}() for i = 1:n]
    updtodate = [BigInt(0) for i = 1:n]
    
    for v in vertices(G)
        size[v] = 1
        pointers[v] = push!(sets[1], v)
    end
    
    maxcard = 1
    i = 1
    count = BigInt(0)

    H = Helper(sets, pointers, size, mcsorder, lastpush, succ, updtodate, maxcard, i, count)


    recenum(G, H)

    return H.count
end
