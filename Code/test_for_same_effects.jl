using GraphPlot

include("genData.jl")
include("causal_functions.jl")
include("plotting_fun.jl")
include("ida.jl")

# Random.seed!(6)
# s = 10000
# reps = 1000
# size_tree = 300
# num_nodes = 9
# sub_tree_size = 100

G = SimpleDiGraph(12)
add_edge!(G, 1, 2)
add_edge!(G, 1, 3)
add_edge!(G, 1, 4)
add_edge!(G, 1, 7)
add_edge!(G, 1, 8)
add_edge!(G, 2, 3)
add_edge!(G, 2, 4)
add_edge!(G, 2, 5)
add_edge!(G, 2, 6)
add_edge!(G, 2, 8)
add_edge!(G, 3, 4)
add_edge!(G, 3, 5)
add_edge!(G, 3, 6)
add_edge!(G, 3, 9)
add_edge!(G, 4, 12)
add_edge!(G, 5, 6)
add_edge!(G, 5, 10)
add_edge!(G, 6, 11)
x = 1
y = 12

# G = SimpleDiGraph(4)
# add_edge!(G, 1, 2)
# add_edge!(G, 1, 4)
# add_edge!(G, 2, 3)
# add_edge!(G, 2, 4)
# wg = assignweights(G)
# ts = collect(1:nv(G))
# x = 1
# y = nv(G)

# make random cordal graph
# G = randCordGraph(size_tree, num_nodes, sub_tree_size)
# ts = randperm(num_nodes)
# makeDAG!(G, ts)
# wg = assignweights(G)
# num_nodes = 6
# en = 5
# G = randgraph(num_nodes, en, "er")
ts = collect(1:nv(G))
wg = assignweights(G)
dt = sampledata(wg, 90, ts)
df = DataFrame(dt, :auto)
true_cov = get_cov(G, ts, wg)
# est_cov = cov(Matrix(df))
cpDAG = cpdag(G)
# nodelabel = 1:nv(G)
# x = 1
# y = num_nodes


function get_effects_adj_sets(df, x, y, G)

    parents = inneighbors(G, x)
    cpDAG = cpdag(G)

    est_effects = DataFrame(Parents = String[], Effect=Float64[], Multip = Float64[], IDA=String[])
    true_effect = 0
        
    effects, multips = multiplicities(cpDAG, x, y, df)
    effects_loc = loc_ida(cpDAG, x, y, df)

    for k in keys(effects)
        if rep == 1
            push!(est_effects,(string(k), get(effects, k, 0), multips[k], "global IDA"))
            push!(est_effects,(string(k), get(effects, k, 0), 1, "local IDA"))
        else
            est_effects[est_effects.Parents .== string(k), :Effect] .+= get(effects, k, 0)
        end

    end

    resp = Term(Meta.parse("x" * string(y)))
    pr = ["x" * string(x)]
    for z in parents
        push!(pr, "x" * string(z))
    end
    pr = map(x -> Meta.parse(x), pr)
    pred = Tuple(Term.(pr))
    true_reg = coef(lm(FormulaTerm(resp, pred), df))[2]
    true_effect += true_reg

    return est_effects, true_effect

end


function get_true_diff_parentSets(multips, cp, x, y)

    parentSets = collect(keys(multips))
    multips_new =copy(multips)
    pSets  = DataFrame(parentSet = Vector{Int64}[], equalSets = String[], multip_glob = Float64[], multip_loc = Float64[])

    for i in 1:length(parentSets)
        if parentSets[i] ∉ keys(multips_new)
            continue
        else
            for j in (i+1):length(parentSets)
                if parentSets[j] ∉ keys(multips_new)
                    continue
                elseif is_adj_set(parentSets[i], parentSets[j], cp, x, y)
                    if parentSets[i] in pSets.parentSet
                        pSets[pSets.parentSet .== [parentSets[i]], :equalSets] .= string(pSets[pSets.parentSet .== [parentSets[i]], :].equalSets[1], ", ", parentSets[j])
                        pSets[pSets.parentSet .== [parentSets[i]], :multip_glob] .= pSets[pSets.parentSet .== [parentSets[i]], :multip_glob][1] + multips[parentSets[j]]
                        pSets[pSets.parentSet .== [parentSets[i]], :multip_loc] .= pSets[pSets.parentSet .== [parentSets[i]], :multip_loc][1] + 1
                    else
                        push!(pSets,(parentSets[i], string(parentSets[j]), multips[parentSets[i]] + multips[parentSets[j]], 2))
                    end
                    multips_new[parentSets[i]] += multips[parentSets[j]]
                    delete!(multips_new, parentSets[j])
                end
            end
            if !(parentSets[i] in pSets.parentSet)
                push!(pSets,(parentSets[i], "none", multips[parentSets[i]], 1))
                #?
                delete!(multips_new, parentSets[i])
            end
        end
    end

    return pSets
end


function is_adj_set(a, b, cp, x, y)

    if y in a && y in b
        return true
    end

    B = adjustGraph(cp, x, b)

    if check_first(B, a, x, y)
        if check_second(B, a, x, y)
            return true
        else
            return false
        end
    else
       return false
    end
end


function adjustGraph(cp, x, parents)
    G = copy(cp)
    for i in all_neighbors(G, x)
        if i in parents
            rem_edge!(G, x, i)
        else
            rem_edge!(G, i, x)
        end
    end

    meek!(G)

    return G
end


function check_first(G, Z, x, y)
    A = copy(G)

    for i in outneighbors(A, x)
        if has_path(A, i, y)
            for z in Z
                if has_path(A, i, z)
                    return false
                end
            end
        end
    end
    return true
end


function check_second(G, Z, x, y)

    function visit(f, j)
        if !marked[[f, j]] && j != x
            enqueue!(q, [f, j])
            marked[[f, j]] = true
        end
    end

    n = nv(G)
    marked = Dict()
    q = Queue{Vector{Int64}}()

    g = copy(G)
    for v in outneighbors(G, x)
        if !(v in inneighbors(G, x)) && has_path(G, v, y)
            rem_edge!(g, x, v)
        end
    end

    for node in 1:n
        for neighbor in inneighbors(g, node)
            marked[[node, neighbor]] = false
        end
        for neighbor in outneighbors(g, node)
            marked[[node, neighbor]] = false
        end
    end
    

    for v in outneighbors(g, x)
        visit(x, v)
    end
    for v in inneighbors(g, x)
        visit(x, v)
    end

    while length(q) > 0
        f, k = dequeue!(q)

        if k == y
            return false
        end


        # f - k
        if !(k in Z) && f in inneighbors(g, k) && f in outneighbors(g, k)
            for node in outneighbors(g, k)
                visit(k, node)
            end
        end

        # f -> k
        if k in Z && f in inneighbors(g, k) && !(f in outneighbors(g, k))
            for node in inneighbors(g, k)
                visit(k, node)
            end
        end

        if !(k in Z) && f in inneighbors(g, k) && !(f in outneighbors(g, k))
            for node in outneighbors(g, k)
                if !(node in inneighbors(g, k))
                    visit(k, node)
                end
            end
        end

        # f <- k
        if !(k in Z) && !(f in inneighbors(g, k)) && f in outneighbors(g, k)
            for node in inneighbors(g, k)
                visit(k, node)
            end
            for node in outneighbors(g, k)
                visit(k, node)
            end
        end
    end

    return true
end


function effects_multip_for_parents(pSets, effects)
    est_effects = DataFrame(Parents = Vector{Int64}[], Effect=Float64[], Multip = Float64[], IDA=String[])

    for parents in eachrow(pSets)
        println(get(effects, parents.parentSet, 0))
        push!(est_effects,(parents.parentSet, get(effects, parents.parentSet, 0), parents.multip_glob, "global IDA"))
        push!(est_effects,(parents.parentSet, get(effects, parents.parentSet, 0), parents.multip_loc, "local IDA"))
    end

    return est_effects
end


multips =  multiplicities_no_effects(cpDAG, x)
pSets = get_true_diff_parentSets(multips, cpDAG, x, y)
# effects = get_effect_parents(pSets, cpDAG, x, y, wg, 70, ts)
effects, multips = multiplicities(cpDAG, x, y, df)

est_effects, true_effect = get_effects_adj_sets(df, x, y, G)
est_effects_parents = effects_multip_for_parents(pSets, effects)

# nodelabel = 1:nv(cpDAG)
# gplot(cpDAG, nodelabel=nodelabel)


# est_effects = compute_effects_one_rep(cpDAG, x, y, 1, true_cov, df, 0.1, 7)