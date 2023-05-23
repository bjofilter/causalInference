# Apply Meek-Rules 1-4 to a CPDAG
function meek!(graph)
    n = nv(graph)
    done = false
    while !done
        done = true
        # rule 1:
        # a -> b - c => a -> b -> c
        for a = 1:n, b in outneighbors(graph, a)
            if !has_edge(graph, b, a)
                for c in outneighbors(graph, b)
                    if has_edge(graph, c, b) && !has_edge(graph, a, c) && !has_edge(graph, c, a)
                        rem_edge!(graph, c, b)
                        done = false
                    end
                end
            end
        end

        # rule 2:
        # a -> b -> c and a - c => a -> c
        for a = 1:n, b in outneighbors(graph, a)
            if !has_edge(graph, b, a)
                for c in outneighbors(graph, b)
                    if !has_edge(graph, c, b) && has_edge(graph, a, c) && has_edge(graph, c, a)
                        rem_edge!(graph, c, a)
                        done = false
                    end
                end
            end
        end

        # rule 3:
        # a - b  and a - c -> b and a - d -> b => a -> b
        for a = 1:n, b in outneighbors(graph, a)
            if has_edge(graph, b, a)
                for c in inneighbors(graph, b), d in inneighbors(graph, b)
                    if !has_edge(graph, b, c) && !has_edge(graph, b, d) &&
                        has_edge(graph, a, c) && has_edge(graph, c, a) &&
                        has_edge(graph, a, d) && has_edge(graph, d, a) &&
                        rem_edge!(graph, b, a)
                        done = false
                    end
                end
            end
        end

        # rule 4:
        # a -> b -> c and a - d - c and b - d  => d -> c
        for a = 1:n, b in outneighbors(graph, a)
            if !has_edge(graph, b, a)
                for c in outneighbors(graph, b)
                    if !has_edge(graph, c, b)
                        for d in outneighbors(graph, a)
                            if has_edge(graph, d, a) &&
                                has_edge(graph, b, d) && has_edge(graph, d, b) &&
                                has_edge(graph, d, c) && has_edge(graph, c, d)
                                rem_edge!(graph, c, d)
                                done = false
                            end
                        end
                    end
                end
            end
        end
    end
end

# Given a CPDAG, recurse over each undirected edge, direct it both ways, appl meek rules, pick next edge and so on, until the set 
# of all possible dags is computed
function orient_cpdag(graph)
    meek!(graph)
    n = nv(graph)
    done = false
    called_anew = false
    while !done
        done = true
        for a = 1:n, b in outneighbors(graph, a)
            if has_edge(graph, b, a)
                graph_b = copy(graph)
                rem_edge!(graph, b, a)
                rem_edge!(graph_b, a, b)
                # println(a, "->", b)
                orient_cpdag(graph)
                # println(a, "<-", b)
                orient_cpdag(graph_b)
                done = false
                called_anew = true
            end        
        end
        if done && !called_anew
            push!(dags, graph)
        end
    end
    return dags
end


function get_parent_sets(dags, x)
    parent_sets = Dict()
    for i in 1:length(dags)
        parents = inneighbors(dags[i], x)
        if parents in keys(parent_sets)
            parent_sets[parents] += 1
        else
            merge!(parent_sets, Dict(parents => 1))
        end
    end
    return parent_sets
end