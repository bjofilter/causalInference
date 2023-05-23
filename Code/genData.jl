using Random
using Graphs
using Distributions
using CausalInference
using DataFrames
using StatsModels
using GLM
using CSV
using DataStructures

include("activelearning.jl")

function randgraph(n, d, dist)
    if dist == "er"
        G = SimpleDiGraph(n)
        D = erdos_renyi(n, d)
        for a = 1:n, b in neighbors(D, a)
            add_edge!(G, a, b)
            add_edge!(G, b, a)
        end
    elseif dist == "ba"
        G = SimpleDiGraph(n)
        D = Graphs.barabasi_albert!(Graphs.random_regular_digraph(Int64(round(n*d)), Int64(round(n*d*d))), n, Int64(round((n-1)*d)))
        for a = 1:n, b in neighbors(D, a)
            add_edge!(G, a, b)
            add_edge!(G, b, a)
        end
    end
    return G
end

function random_Tree(size_tree)
    G = SimpleDiGraph(size_tree)
    for i in 2:size_tree
        r = rand(1:i-1)
        add_edge!(G, i, r)
        add_edge!(G, r, i)
    end
    return G
end

function gen_random_subtrees(G, num_subtrees, sub_tree_size)
    num_nodes = nv(G)
    subtrees = Array{Set}(undef, num_subtrees)
    for i in 1:num_subtrees
        start = rand(1:num_nodes)
        subtree = Set(start)
        border = Set(start)

        for j in 2:sub_tree_size
            subtree, border = pick_random(G, subtree, border)
        end
        subtrees[i] = subtree

    end
    return subtrees
end

function pick_random(G, subtree, border)
    border_node = rand(1:length(border))
    border_node = collect(border)[border_node]
    neighs = all_neighbors(G, border_node)
    neighs = neighs[(!in).(neighs,Ref(collect(subtree)))]
    node = rand(1:length(neighs))
    node = neighs[node]
    push!(subtree, node)
    
    if length(neighs) == 1
        pop!(border, border_node)
    end

    new_neighs = all_neighbors(G, node)
    if  length(new_neighs[(!in).(new_neighs,Ref(collect(subtree)))]) > 0
        push!(border, node)
    end

    return subtree, border

end

function randCordGraph(size_tree, num_nodes, sub_tree_size)
    U = random_Tree(size_tree)

    subtrees = gen_random_subtrees(U, num_nodes, sub_tree_size)

    G = SimpleDiGraph(num_nodes)
    for i in 1:num_nodes
        for j in 1:i-1
            if !isempty(intersect(subtrees[i], subtrees[j]))
                add_edge!(G, i, j)
                add_edge!(G, j, i)
            end
        end
    end
    return G
end

# make random cordal graph
function make_cord_graph(num_nodes, density, range)
    if density  == 1
        sub_tree_size = 40
    else
        sub_tree_size = density * 20 + 2
    end
    size_tree = 100
    cpDAG = randCordGraph(size_tree, num_nodes, sub_tree_size)
    while ((ne(cpDAG) / (nv(cpDAG)*(nv(cpDAG)-1))) > density + range || (ne(cpDAG) / (nv(cpDAG)*(nv(cpDAG)-1))) < density - range)
        cpDAG = randCordGraph(size_tree, num_nodes, sub_tree_size)
    end
    G = sampleDAG(cpDAG)
    wg = assignweights(G)
    ts = ts_from_DAG(G)
    return G, cpDAG, wg, ts
end

function randDAG!(G)
    n = nv(G)
    ts = randperm(n)
    makeDAG!(G, ts)
end

function makeDAG!(G, ts)
    n = nv(G)
    for a in 1:n, b in inneighbors(G, a)
        if ts[b] < ts[a]
            rem_edge!(G, a, b)
        end
    end
end

function assignweights(G)
    n = nv(G)
    wg = [Dict{Integer, AbstractFloat}() for i = 1:n]
    ud = Uniform(1, 2)
    for a in vertices(G), b in outneighbors(G, a)
        wg[a][b] = rand(ud)
    end
    return wg
end

function sampledata(wg, s, ts)
    n = size(wg, 1)
    nd = Normal()
    dt = rand(nd, s, n)
    # for i in 1:n
    #     while abs(1 - var(dt[1:s, i])) > 0.00001
    #         if 1 - var(dt[1:s, i]) > 0
    #             dt[1:s, i] .*= 1 + abs(1 - var(dt[1:s, i]))/10
    #         else
    #             dt[1:s, i] .*= 1 - abs(1 - var(dt[1:s, i]))/10
    #         end
    #     end
    # end

    its = zeros(Int64, n)
    for i in 1:n
        its[ts[i]] = i
    end
    for a in its, b in wg[a]
        dt[1:s, b.first] += dt[1:s, a] * b.second
    end
    return dt
end

function make_graph_data(graph, num_nodes, density, s)
    if graph == "Gil"
        G = randgraph(num_nodes, density, "er")
        ts = randperm(num_nodes)
        makeDAG!(G, ts)
        wg = assignweights(G)
        cpDAG = cpdag(G)

    elseif graph == "BA"
        G = randgraph(num_nodes, density, "ba")
        ts = randperm(num_nodes)
        makeDAG!(G, ts)
        wg = assignweights(G)
        cpDAG = cpdag(G)

    else graph == "Cho"
        range_d = 0.02
        G, cpDAG, wg, ts = make_cord_graph(num_nodes, density, range_d)
        while (!has_path(cpDAG, 1, num_nodes))
            G, cpDAG, wg, ts = make_cord_graph(num_nodes, density, range_d)
        end
    end

    dt = sampledata(wg, s, ts)
    df = DataFrame(dt, :auto)
    true_cov = get_cov(G, ts, wg)

    return cpDAG, df, true_cov
end

function generateidadata(n, d, s)
    G = randgraph(n, d, "er")
    ts = randperm(n)
    makeDAG!(G, ts)
    wg = assignweights(G)
    return sampledata(wg, s, ts)
end

function ts_from_DAG(G)
    ts = zeros(nv(G))
    del = zeros(nv(G))
    D = copy(G)
    num = 1
    while 0 in ts
        for vert in vertices(D)
            if (inneighbors(D, vert) == Int64[] && !(vert in del))
                while !(outneighbors(D, vert) == Int64[])
                    rem_edge!(D, vert, outneighbors(D, vert)[1])
                end
                ts[vert] = num
                del[num] = vert
                num += 1
                break
            end
        end
    end
    ts = convert(Array{Int64}, ts)
    return ts
end