using GraphPlot
using Statistics
using LinearSolve
using BenchmarkTools

include("genData.jl")
include("ida.jl")
include("plotting_fun.jl")
include("causal_functions.jl")

#____________________________________________________________________________________

function compute_effects_one_rep(cpDAG, df, x, y, true_cov)

    est_effects = DataFrame(Effect=Float64[], Multip=Float64[], IDA=String[])
        
    effects, multips, true_eff = multiplicities(cpDAG, x, y, true_cov, df)

    for k in keys(effects)
        push!(est_effects,(get(effects, k, 0), 1, "local IDA"))
        push!(est_effects,(get(effects, k, 0), multips[k], "global IDA"))
        push!(est_effects,(get(true_eff, k, 0), multips[k], "true IDA"))
    end

    return est_effects

end

function prep_for_dense(effects)
    effects_for_dense = DataFrame(Effect=Float64[], IDA=String[])
    effects_g = effects[effects.IDA .== "global IDA", :]
    min = minimum(effects_g.Multip)
    while min > 100000
        effects_g.Multip = effects_g.Multip ./ sqrt(min)
        min = minimum(effects_g.Multip)
    end
    for j in 1:size(effects_g)[1]
        push!(effects_for_dense,(effects_g[j, :Effect], "local IDA"))
        for k in 1:round(effects_g[j, :Multip])
            push!(effects_for_dense,(effects_g[j, :Effect], "global IDA"))
        end
    end
    return effects_for_dense
end

function prep_for_hist(effects)
    effects_for_dense = DataFrame(Effect=Float64[], IDA=String[])
    effects_g = effects[effects.IDA .== "global IDA", :]
    min = minimum(effects_g.multip)
    while min > 1000
        effects_g.multip = effects_g.multip ./ sqrt(min)
        min = minimum(effects_g.multip)
    end
    for j in 1:size(effects_g)[1]
        for k in 1:round(effects_g[j, :multip])
            push!(effects_for_dense,(effects_g[j, :Effect], "global IDA"))
        end
    end
    len = size(effects_for_dense)[1]
    num_effects = size(effects)[1] / 2
    num_loc = len / num_effects
    for j in 1:size(effects_g)[1]
        for k in 1:round(num_loc)
            push!(effects_for_dense,(effects_g[j, :Effect], "local IDA"))
        end
    end
    return effects_for_dense
end

function check(graph)
    res = DataFrame(Num_Effects=Float64[], Density=Float64[], Method=String[], Size=Int64[])
    x = 1
    s = 2000
    densities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    num_nodes = [10, 25, 50, 75, 100, 1000]
    for n in num_nodes
        for d in densities
            for j in 1:1000
                println("n=", n, "; d=", d, "; rep=", j)

                cpDAG, df, true_cov = make_graph_data(graph, n, d, s)
                num_effs = length(multiplicities_no_effects(cpDAG, x))
                push!(res,(num_effs, d, graph, n))

            end
        end
    end
    return res
end

function count_undirec_edges(graph)
    num = [10, 25, 50, 100, 1000]
    res = DataFrame(Undirec_Edges=Float64[], Density=Float64[], Size = Int64[], Method=String[])
    for num_nodes in num
        for d in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
            for j in 1:1000
                println("d=", d, "; rep=", j)

                cpDAG, df, true_cov = make_graph_data(graph, n, d, s)
                undirec_edges = 0
                for k in 1:n
                    undirec_edges += size(setdiff(inneighbors(cpDAG, k), setdiff(inneighbors(cpDAG, k), outneighbors(cpDAG, k))))[1]
                end
                undirec_edges /= ne(cpDAG)
                push!(res,(undirec_edges, d, n, graph))
            end
        end
    end
    return res
end

function plot_num_effects(res, graph)
    num_er = res[res.Method .== graph, :]
    direc_er = undirec_edges[undirec_edges.Method .== graph, :]
    num_er_sig = num_er[num_er.Num_Effects .> 1, :]

    perc = DataFrame(Perc=Float64[], Density=Float64[], Size=Int64[])

    for s in unique(num_er.Size)
        n = num_er[num_er.Size .== s, :]
        for i in unique(num_er.Density)
            a = n[n.Density .== i, :]
            b = a[a.Num_Effects .!= 1, :]
            per = size(b)[1] / size(a)[1]
            push!(perc, (per, i, s))
        end
    end

    @df perc StatsPlots.plot(:Density, :Perc, 
                    group=:Size,
                    xlabel = "Density",
                    ylabel = "Amount With Multiple Possible Effects",
                    grid=false,
                    legend=:topleft,
                    # palette = :Spectral_5,
                    linewidth=2,
                    tick_direction=:out,
                    framestyle = :box)


    direc_e = DataFrame(Perc=Float64[], Density=Float64[], Size=Int64[])

    for s in unique(undirec.Size)
        n = undirec[undirec.Size .== s, :]
        for i in unique(undirec.Density)
            a = n[n.Density .== i, :]
            a = filter(row -> ! isnan(row.Undirec_Edges), a)
            av = mean(a.Undirec_Edges)
            push!(direc_e, (av, i, s))
        end
    end

    @df direc_e StatsPlots.plot(:Density, :Perc, 
                    group=:Size,
                    xlabel = "Density",
                    ylabel = "Average Amount of Undirected Edges\nin the CPDAG",
                    grid=false,
                    legend=:topright,
                    linewidth=2,
                    tick_direction=:out,
                    framestyle = :box)

    StatsPlots.boxplot(num_er.Density, num_er.Num_Effects, color="white",
                                            xlabel = "Density",
                                            ylabel = "Number of {{X}}-MECs",
                                            grid=false,
                                            yaxis=:log,
                                            linecolor = "black",
                                            legend=false,
                                            outliers=false,
                                            tick_direction=:out,
                                            framestyle = :box,
                                            yticks = [1, 10^4, 10^8, 10^12, 10^16, Float64(10)^20,  Float64(10)^24, Float64(10)^28, Float64(10)^32],
                                            xticks=(1:9, ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]),
                                            whisker_width=:half,
                                            dpi=400)
    savefig("Num_er.pdf")  


    @df direc_er StatsPlots.boxplot(:Density, :Undirec_Edges, 
                                            color="white",
                                            linecolor = "black",
                                            xlabel = "Density",
                                            ylabel = "Amount of Undirected Edges",
                                            grid=false,
                                            legend=false,
                                            outliers=true,
                                            tick_direction=:out,
                                            framestyle = :box,
                                            xticks=xticks = (range(10, step=10, length=9), range(0.1, step=0.1, length=9)),
                                            whisker_width=:half,
                                            dpi=400)
    savefig("Undirec_er.pdf")  
end





#____________________________________________________________________________________

# res = check("Gil")
# undirec_edges = count_undirec_edges()



function check_num_effects()
    res = DataFrame(Num_Effects=Float64[], Size=Int64[], Rep = Int64[])
    nodes = [20, 40, 60, 80, 100]
    x = 1
    d = 0.1
    for n in nodes
        for i in 1:1000
            println(n, " : ", i, "/1000")
            G, cpDAG, wg, ts = make_cord_graph(n, d, 0.01)
            num_effs = Float64(2)^Float64(size(neighbors(cpDAG, x))[1])

            push!(res,(num_effs, n, i))
        end
    end
    return res
end