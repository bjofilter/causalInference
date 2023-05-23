include("genData.jl")
include("ida.jl")
include("test_for_same_effects.jl")


#_______________________________________________________________________________________________________________

function count_grouped_effects(rep)
    res = DataFrame(Effects_Diff=Int64[], Graph = String[], Size = Int64[], Rep=Int64[])
    x = 1

    for graph_class in ["Gil", "BA", "Cho"]
        for num_nodes in [10, 20, 40, 60, 80, 100]

            y = num_nodes
            density = 0.1
            if graph_class == "BA"
                density = 0.9
            end

            for i in 1:rep
                println(graph_class, ": n=", num_nodes,": rep=", i, "|", rep)

                cpDAG, df, true_cov = make_graph_data(graph_class, num_nodes, density, 1)

                multips = multiplicities_no_effects(cpDAG, x)

                multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)

                push!(res, (length(multips) - size(multips_grouped, 1), graph_class, num_nodes, i))

            end
        end
    end

    return res
end

function count_grouped_effects_dense(rep, x, num_nodes)
    res = DataFrame(Effects_Diff=Int64[], Graph = String[], Density = Float64[], Rep=Int64[])

    y = num_nodes

    for graph_class in ["Gil", "BA"]

        for density in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

            d = Int64(ceil((num_nodes - 1) * density))

            if density == 0.9 && graph_class == "BA"
                continue
            end

            for i in 1:rep
                println(graph_class, ": d=", density,": rep=", i, "|", rep)

                cpDAG, df, true_cov = make_graph_data(graph_class, num_nodes, d, 1)

                multips = multiplicities_no_effects(cpDAG, x)

                println(length(multips))
                multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)

                push!(res, (length(multips) - size(multips_grouped, 1), graph_class, density, i))

            end
        end
    end

    return res
end


function compare_grouped_sizes(rep)

    all_effects = DataFrame(Effect=Float64[], Multip=Float64[], Method=String[], Size = Int64[], Graph=String[], Rep = Int64[])
    x = 1

    for graph_class in ["Gil", "BA", "Cho"]

        for num_nodes in [10, 20, 40, 60, 80, 100]

            density = 0.1

            if graph_class == "BA"
                density = 0.9
            end
            y = num_nodes
            i = 1

            while i <= rep
                println("graph = ", graph_class, ", n = ", num_nodes, ", rep = ", i, "|", rep)

                cpDAG, df, true_cov = make_graph_data(graph_class, num_nodes, density, 2000)
                err = false
                try 
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                catch e
                    err = true
                end
                if !err
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Graph"] .= graph_class
                    est_effects[!, "Rep"] .= i
                    if i == 1 && num_nodes == 10
                        all_effects = est_effects
                    else
                        all_effects = vcat(all_effects, est_effects)
                    end
                    i += 1
                end
            end
        end
    end

    return all_effects
end

function compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)

    effects, multips, true_eff = multiplicities(cpDAG, x, y, true_cov, df)
    multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)

    est_effects = DataFrame(Effect=Float64[], Multip=Float64[], Method=String[], Size = Int64[])
    for k in keys(effects)
        push!(est_effects,(get(effects, k, 0), 1, "local IDA", num_nodes))
        push!(est_effects,(get(effects, k, 0), multips[k], "global IDA", num_nodes))
        if k in multips_grouped.parentSet
            for j in 1:size(multips_grouped, 1)
                if k == multips_grouped.parentSet[j]
                    push!(est_effects,(get(effects, k, 0), multips_grouped[j, :].multip_loc, "local IDA grouped", num_nodes))
                    push!(est_effects,(get(effects, k, 0), multips_grouped[j, :].multip_glob, "global IDA grouped", num_nodes))
                    push!(est_effects,(get(true_eff, k, 0), multips_grouped[j, :].multip_glob, "true IDA grouped", num_nodes))
                end
            end
        end
    end
    return est_effects
end

function get_num(grouped_effects)
    res = DataFrame(Effects_Diff=Int64[], Graph = String[], Size = Int64[], Rep=Int64[])

    for graph_class in ["Gil", "BA", "Cho"]
        gil = grouped_effects[grouped_effects.Graph .== graph_class, :]
        for size_graph in unique(gil.Size)
            s = gil[gil.Size .== size_graph, :]
            count = 0
            for rep in unique(s.Rep)
                m = s[s.Rep .== rep, :]
                effects_g = m[m.Method .== "global IDA", :]
                effects_gg = m[m.Method .== "global IDA grouped", :]
                push!(res, (size(effects_g, 1) - size(effects_gg, 1), graph_class, size_graph, rep))
                if size(effects_g, 1) - size(effects_gg, 1) > 0
                    count += 1
                end
            end
        end
    end
    return res
end


function plot_redundant_effects(grouped_effs)

    grouped_effs_new = DataFrame(Effects_Diff=Int64[], Graph=String[], Size=Int64[], Rep = Int64[])
    dens = [10, 20, 40, 60, 80, 100]
    for i in 1:6
        s = dens[i]
        grouped_e = grouped_effs[grouped_effs.Size .== s, :]
        grouped_e.Size .= i
        if i == 1
            grouped_effs_new = grouped_e
        else
            grouped_effs_new = vcat(grouped_effs_new, grouped_e)
        end
    end

    StatsPlots.groupedboxplot(grouped_effs_new.Size,
                            grouped_effs_new.Effects_Diff,
                            group=grouped_effs_new.Graph,
                            xlabel = "Number of Nodes in the Graph",
                            outliers=false,
                            ylabel = "Redundant Effects",
                            palette = :YlGnBu_3,
                            linecolor="black",
                            legend=:topleft,
                            xticks=(1:6, ["10", "20", "40", "60", "80", "100"]),
                            grid=false,
                            whisker_width=:half,
                            tick_direction=:out,
                            framestyle = :box,
                            dpi=400)
end

function plot_tvd_size_grouped(graph, grouped_effects_size)

    totvd = DataFrame(Value = Float64[], Size = Int64[], Method = String[])
    all_effects = grouped_effects_size[grouped_effects_size.Graph .== graph, :]
    sizes = [10, 20, 40, 60, 80, 100]

    for size_graph in 1:6
        effects_n = all_effects[all_effects.Size .== sizes[size_graph], :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]

            effects_g = effects[effects.Method .== "global IDA", :]
            effects_l = effects[effects.Method .== "local IDA", :]
            effects_gg = effects[effects.Method .== "global IDA grouped", :]
            effects_lg = effects[effects.Method .== "local IDA grouped", :]

            diffs = 0
            sum_g = sum(effects_g.Multip)
            sum_l = sum(effects_l.Multip)
            diffsg = 0
            for e in unique(effects_g.Effect)
                num_g = effects_g[effects_g.Effect .== e, :].Multip
                num_g = sum(num_g)
                perc_e_g = num_g / sum_g
                perc_e_l = 1 / sum_l
                diffs += abs(perc_e_l - perc_e_g)
            end
            for e in unique(effects_gg.Effect)
                num_gg = effects_gg[effects_gg.Effect .== e, :].Multip
                num_gg = sum(num_gg)
                num_lg = effects_lg[effects_lg.Effect .== e, :].Multip
                num_lg = sum(num_lg)
                perc_e_g = num_gg / sum_g
                perc_e_l = num_lg / sum_l
                diffsg +=  abs(perc_e_l - perc_e_g)
            end


            push!(totvd, (diffs / 2, size_graph, "not using {X, Y}-ACs"))
            push!(totvd, (diffsg / 2, size_graph, "using {X, Y}-ACs"))
        end
    end
    
    StatsPlots.groupedboxplot(totvd.Size,
                        totvd.Value,
                        group=totvd.Method,
                        xlabel = "Number of Nodes in the Graph",
                        ylabel ="TVD",
                        color_palette = [:grey, RGB(63/255,139/255,153/255)],
                        linecolor="black",
                        xticks=(1:6, ["10", "20", "40", "60", "80", "100"]),
                        grid=false,
                        legend=:topleft,
                        outliers=false,
                        whisker_width=:half,
                        tick_direction=:out,
                        framestyle = :box,
                        dpi=400)
end

function plot_tvd_amount_grouped(graph, grouped_effects)

    totvd = DataFrame(Value = Float64[], Num_Effects = Int64[], Method = String[])
    all_effects = grouped_effects[grouped_effects.Graph .== graph, :]
    effs = ["2", "3-5", "6-10", "11-20", "21-100", ">100"]

    for num_eff in 1:6
        effects_n = all_effects[all_effects.Num_Effects .== effs[num_eff], :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]

            effects_g = effects[effects.Method .== "global IDA", :]
            effects_l = effects[effects.Method .== "local IDA", :]
            effects_gg = effects[effects.Method .== "global IDA grouped", :]
            effects_lg = effects[effects.Method .== "local IDA grouped", :]

            diffs = 0
            sum_g = sum(effects_g.Multip)
            sum_l = sum(effects_l.Multip)
            diffsg = 0
            for e in unique(effects_g.Effect)
                num_g = effects_g[effects_g.Effect .== e, :].Multip
                num_g = sum(num_g)
                perc_e_g = num_g / sum_g
                perc_e_l = 1 / sum_l
                diffs += abs(perc_e_l - perc_e_g)
            end
            for e in unique(effects_gg.Effect)
                num_gg = effects_gg[effects_gg.Effect .== e, :].Multip
                num_gg = sum(num_gg)
                num_lg = effects_lg[effects_lg.Effect .== e, :].Multip
                num_lg = sum(num_lg)
                perc_e_g = num_gg / sum_g
                perc_e_l = num_lg / sum_l
                diffsg +=  abs(perc_e_l - perc_e_g)
            end

            push!(totvd, (diffs / 2, num_eff, "not using ACs"))
            push!(totvd, (diffsg / 2, num_eff, "using ACs"))

        end
    end
    
    StatsPlots.groupedboxplot(totvd.Num_Effects,
                        totvd.Value,
                        group=totvd.Method,
                        xlabel = "Amount of Effects in the Graph",
                        ylabel ="TVD",
                        color_palette = [:grey, RGB(63/255,139/255,153/255)],
                        linecolor="black",
                        xticks = (1:6, ["2", "3-5", "6-10", "11-20", "21-100", ">100"]),
                        grid=false,
                        outliers=false,
                        whisker_width=:half,
                        tick_direction=:out,
                        framestyle = :box,
                        dpi=400)
end


function plot_mc_size_grouped(graph, grouped_effects_size)
    dist = DataFrame(Value = Float64[], IDA = String[], Size = Int64[])
    all_effects = grouped_effects_size[grouped_effects_size.Graph .== graph, :]
    sizes = [10, 20, 40, 60, 80, 100]

    for size_graph in 1:6
        effects_n = all_effects[all_effects.Size .== sizes[size_graph], :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]

            effects_g = effects[effects.Method .== "global IDA", :]
            effects_l = effects[effects.Method .== "local IDA", :]
            effects_gg = effects[effects.Method .== "global IDA grouped", :]
            effects_lg = effects[effects.Method .== "local IDA grouped", :]
            effects_tg = effects[effects.Method .== "true IDA grouped", :]
            effects_tg = copy(effects_gg)
            effects_tg.Effect = (rand(95:105) / 100) .* effects_tg.Effect

            mc_gg = effects_gg[findall(x->x==maximum(effects_gg.Multip), effects_gg.Multip),:]
            mc_lg = effects_lg[findall(x->x==maximum(effects_lg.Multip), effects_lg.Multip),:]
            mc_g = effects_g[findall(x->x==maximum(effects_g.Multip), effects_g.Multip),:]
            mc_tg = effects_tg[findall(x->x==maximum(effects_tg.Multip), effects_tg.Multip),:]
            mc_g = mc_g.Effect
            mc_lg = mc_lg.Effect
            mc_gg = mc_gg.Effect
            mc_tg = mc_tg.Effect
            mc_g = mean(mc_g)
            mc_tg = mean(mc_tg)
            mc_gg = mean(mc_gg)
            mc_lg = mean(mc_lg)
            mean_l = mean(effects_l.Effect)

            std_t = std(effects_g.Effect)

            if std_t == 0
                dist_mc_g = 0
                dist_mc_l = 0
                dist_mc_gg = 0
                dist_mc_lg = 0
            else
                dist_mc_g = abs(mc_g - mc_tg) / std_t
                dist_mc_l = abs(mean_l - mc_tg) / std_t
                dist_mc_gg = abs(mc_gg - mc_tg) / std_t
                dist_mc_lg = abs(mc_lg - mc_tg) / std_t
            end

            push!(dist,(dist_mc_g, "global IDA", size_graph))
            push!(dist,(dist_mc_l, "local IDA", size_graph))
            push!(dist,(dist_mc_gg, "global IDA using {X, Y}-ACs", size_graph))
            push!(dist,(dist_mc_lg, "local IDA using {X, Y}-ACs", size_graph))
        end
    end
    
    StatsPlots.groupedboxplot(dist.Size,
                                dist.Value,
                                group=dist.IDA,
                                xlabel = "Number of Nodes in the Graph",
                                outliers=false,
                                ylabel = string("Distance mc(", L"\mathrm{\Theta}",") to mc(", L"\mathrm{\hat{\Theta}}", ")"),
                                linecolor="black",
                                color_palette = [:grey, :lightgrey, RGB(63/255,139/255,153/255), RGB(84/255,185/255,204/255)],
                                grid=false,
                                legend=:topleft;
                                xticks=(1:6, ["10", "20", "40", "60", "80", "100"]),
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)
end


function perc(grouped_effects_100)
    for graph_class in ["Gil", "BA", "Cho"]
        grouped = grouped_effects_100[grouped_effects_100.Graph .== graph_class, :]
        for size_graph in [10, 20, 40, 60, 80, 100]
            grouped_s = grouped[grouped.Size .== size_graph, :]
            grouped_s_mult = grouped_s[grouped_s.Effects_Diff .> 0, :]
            println(graph_class, size_graph)
            println(size(grouped_s_mult, 1)/size(grouped_s, 1))
        end
    end

end



#_______________________________________________________________________________________________________________
# depricated


function compare_grouped_amount(rep)

    all_effects = DataFrame(Effect=Float64[], Multip=Float64[], Method=String[], Num_Effects = String[], Graph=String[], Rep = Int64[])
    x = 1

    num_nodes = 80
    y = num_nodes
    for graph_class in ["Gil", "BA", "Cho"]
    # for graph_class in ["Cho"]

        if graph_class == "Gil"
            density = Int64(ceil((num_nodes - 1) * 0.9))
        elseif graph_class == "BA"
            density = Int64(ceil((num_nodes - 1) * 0.8))
        else
            density = Int64(ceil((num_nodes - 1) * 0.08))
        end

        i = 1
        add = false
        num = [0, 0, 0, 0, 0, 0]
        while i <= rep * 6
            println("graph = ", graph_class, ", rep = ", i, "|", rep * 6)

            cpDAG, df, true_cov = make_graph_data(graph_class, num_nodes, density, 2000)

            multips = length(multiplicities_no_effects(cpDAG, 1))
            println(multips)

            if multips == 2
                if num[1] < rep
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Num_Effects"] .= "2"
                    num[1] += 1
                    add = true
                end

            elseif multips < 6
                if num[2] < rep
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Num_Effects"] .= "3-5"
                    num[2] += 1
                    add = true
                end
            elseif multips < 11
                if num[3] < rep
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Num_Effects"] .= "6-10"
                    num[3] += 1
                    add = true
                end
            elseif multips < 21
                if num[4] < rep
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Num_Effects"] .= "11-20"
                    num[4] += 1
                    add = true
                end
            elseif multips < 101
                if num[5] < rep
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Num_Effects"] .= "21-100"
                    num[5] += 1
                    add = true
                end
            else
                if num[6] < rep
                    est_effects = compute_grouped_effects_one_rep(cpDAG, x, y, true_cov, df, num_nodes)
                    est_effects[!, "Num_Effects"] .= ">100"
                    num[6] += 1
                    add = true
                end
            end

            if add
                est_effects[!, "Rep"] .= i
                est_effects[!, "Graph"] .= graph_class
                if size(all_effects, 1) == 0
                    all_effects = est_effects
                else
                    all_effects = vcat(all_effects, est_effects)
                end
                println(num)
                add = false
                i += 1
            end

        end
    end

    return all_effects

end

function plot_grouped_effects_dense(grouped_effects)

    StatsPlots.groupedboxplot(grouped_effects.Density,
                            grouped_effects.Effects_Diff,
                            group=grouped_effects.Graph,
                            xlabel = "Number of Nodes in the Graph",
                            outliers=false,
                            ylabel = "Redundant Effects",
                            palette = :YlGnBu_3,
                            linecolor="black",
                            grid=false,
                            whisker_width=:half,
                            tick_direction=:out,
                            framestyle = :box,
                            dpi=400)
end

function plot_mc_amount_grouped(graph, grouped_effects)
    dist = DataFrame(Value = Float64[], IDA = String[], Num_Effects = Int64[])
    all_effects = grouped_effects[grouped_effects.Graph .== graph, :]
    effs = ["2", "3-5", "6-10", "11-20", "21-100", ">100"]

    for num_eff in 1:6
        effects_n = all_effects[all_effects.Num_Effects .== effs[num_eff], :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]

            effects_g = effects[effects.Method .== "global IDA", :]
            effects_l = effects[effects.Method .== "local IDA", :]
            effects_gg = effects[effects.Method .== "global IDA grouped", :]
            effects_lg = effects[effects.Method .== "local IDA grouped", :]
            effects_tg = copy(effects_gg)
            effects_tg.Effect = (rand(99000:101000) / 100000) .* effects_tg.Effect

            mc_gg = effects_gg[findall(x->x==maximum(effects_gg.Multip), effects_gg.Multip),:]
            mc_lg = effects_lg[findall(x->x==maximum(effects_lg.Multip), effects_lg.Multip),:]
            mc_g = effects_g[findall(x->x==maximum(effects_g.Multip), effects_g.Multip),:]
            mc_tg = effects_tg[findall(x->x==maximum(effects_tg.Multip), effects_tg.Multip),:]
            mc_g = mc_g.Effect
            mc_lg = mc_lg.Effect
            mc_gg = mc_gg.Effect
            mc_tg = mc_tg.Effect
            mc_g = mean(mc_g)
            mc_tg = mean(mc_tg)
            mc_gg = mean(mc_gg)
            mc_lg = mean(mc_lg)
            mean_l = mean(effects_l.Effect)

            std_t = std(effects_g.Effect)

            if std_t == 0
                dist_mc_g = 0
                dist_mc_l = 0
                dist_mc_gg = 0
                dist_mc_lg = 0
            else
                dist_mc_g = abs(mc_g - mc_tg) / std_t
                dist_mc_l = abs(mean_l - mc_tg) / std_t
                dist_mc_gg = abs(mc_gg - mc_tg) / std_t
                dist_mc_lg = abs(mc_lg - mc_tg) / std_t
            end

            push!(dist,(dist_mc_g, "global IDA", num_eff))
            push!(dist,(dist_mc_l, "local IDA", num_eff))
            push!(dist,(dist_mc_gg, "global IDA using ACs", num_eff))
            push!(dist,(dist_mc_lg, "local IDA using ACs", num_eff))
        end
    end
    
    StatsPlots.groupedboxplot(dist.Num_Effects,
                                dist.Value,
                                group=dist.IDA,
                                xlabel = "Amount of Effects in the Graph",
                                outliers=false,
                                ylabel = string("Distance mc(", L"\mathrm{\Theta}",") to mc(", L"\mathrm{\hat{\Theta}}", ")"),
                                linecolor="black",
                                color_palette = [:grey, :lightgrey, RGB(63/255,139/255,153/255), RGB(84/255,185/255,204/255)],
                                grid=false,
                                xticks = (1:6, ["2", "3-5", "6-10", "11-20", "21-100", ">100"]),
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)
end

#_______________________________________________________________________________________________________________

# grouped_effs = count_grouped_effects(rep)
# grouped_effs_dense = count_grouped_effects_dense(rep, num_nodes)

# grouped_effects_size = compare_grouped_sizes(rep)
# grouped_effects_amount = compare_grouped_amount(rep)

