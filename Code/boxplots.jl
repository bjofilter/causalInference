using GraphPlot
using StatsPlots
using Distributions
using HypothesisTests
using LaTeXStrings
import RDatasets

include("genData.jl")
include("causal_functions.jl")
include("plotting_fun.jl")
include("ida.jl")

#_______________________________________________________________________________________________


function compute_all_effects_diff_amount(graph, num_nodes, density, reps, s)
    all_effects = DataFrame(Effect=Float64[], IDA=String[], Rep=Int64[])
    i = 0
    add = false
    num = [0, 0, 0, 0, 0, 0]
    while i < 6 * reps
        println(i + 1, "|", 6 * reps)

        cpDAG, df, true_cov = make_graph_data(graph, num_nodes, density, s)
    
        x = 1
        y = num_nodes

        multips = length(multiplicities_no_effects(cpDAG, 1))

        err = false



        if multips == 2
            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end
    
            if !err
                if num[1] < reps
                    effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                    effects_rep[!, "Num_Effects"] .= "2"
                    num[1] += 1
                    add = true
                end
            end
        elseif multips < 6
            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end
    
            if !err
                if num[2] < reps
                    effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                    effects_rep[!, "Num_Effects"] .= "3-5"
                    num[2] += 1
                    add = true
                end
            end
        elseif multips < 11
            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end
    
            if !err
                if num[3] < reps
                    effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                    effects_rep[!, "Num_Effects"] .= "6-10"
                    num[3] += 1
                    add = true
                end
            end
        elseif multips < 21
            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end
    
            if !err
                if num[4] < reps
                    effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                    effects_rep[!, "Num_Effects"] .= "11-20"
                    num[4] += 1
                    add = true
                end
            end
        elseif multips < 101
            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end
    
            if !err
                if num[5] < reps
                    effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                    effects_rep[!, "Num_Effects"] .= "21-100"
                    num[5] += 1
                    add = true
                end
            end
        else
            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end
    
            if !err
                if num[6] < reps
                    effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                    effects_rep[!, "Num_Effects"] .= ">100"
                    num[6] += 1
                    add = true
                end
            end
        end

        if add
            if i == 0
                all_effects = effects_rep
            else
                all_effects = vcat(all_effects, effects_rep)
            end
            println(num)
            add = false
            i += 1
        end
    end
    return all_effects
end


function compute_all_effects(graph, density, reps, s)
    all_effects = DataFrame(Effect=Float64[], IDA=String[], Rep=Int64[], Size = Int64[], Density = Float64[])
    n = [10, 20, 40, 60, 80, 100]
    for num_nodes in n
        i = 1
        while i <= reps
            println("n=", num_nodes, " : ", i, "|", reps)

            cpDAG, df, true_cov = make_graph_data(graph, num_nodes, density, s)
        
            x = 1
            y = num_nodes

            err = false

            try
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
            catch e
                err = true
            end

            if !err
                effects_rep = compute_effects_one_rep(cpDAG, x, y, i, true_cov, df, density, num_nodes)
                if i == 1 && num_nodes == n[1]
                    all_effects = effects_rep
                else
                    all_effects = vcat(all_effects, effects_rep)
                end
                i += 1
            end
        end
    end
    return all_effects
end

# Computes the effect of x on y via sampled data
function compute_effects_one_rep(cp, x, y, rep, true_cov, df, density, num_nodes)

    est_effects = DataFrame(Effect=Float64[], Multip=Float64[], IDA=String[], Size = Int64[], Density = Float64[], Rep = Int64[])

    effects, multips, true_eff = multiplicities(cp, x, y, true_cov, df)

    for k in keys(effects)
        push!(est_effects,(get(effects, k, 0), 1, "local IDA", num_nodes, density, rep))
        push!(est_effects,(get(effects, k, 0), multips[k], "global IDA", num_nodes, density, rep))
        push!(est_effects,(get(true_eff, k, 0), multips[k], "true IDA", num_nodes, density, rep))
    end

    return est_effects

end



# diff means
function compute_e_mean(all_effects)
    e_ave = DataFrame(Value = Float64[], IDA = String[], Num_Effects = String[])

    for num_eff in ["2", "3-5", "6-10", "11-20", "21-100", ">100"]
        effects_n = all_effects[all_effects.Num_Effects .== num_eff, :]
        effects_n = effects_n[effects_n.Size .== 40, :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]
            
            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]
                effects_t = effects[effects.IDA .== "true IDA", :]

                size_g = sum(effects_g.Multip)
                size_t = sum(effects_t.Multip)
                size_l = sum(effects_l.Multip)

                total_eff_g = effects_g.Effect .* effects_g.Multip
                total_eff_t = effects_t.Effect .* effects_t.Multip
                total_eff_l = effects_l.Effect

                mean_g = sum(total_eff_g) / size_g
                mean_t = sum(total_eff_t) / size_t
                mean_l = sum(total_eff_l) / size_l

                std_t = get_std(effects_t)

                if std_t == 0
                    difference_g = abs(mean_g - mean_t) / mean_t
                    difference_l = abs(mean_l - mean_t) / mean_t
                else
                    difference_g = abs(mean_g - mean_t) / std_t
                    difference_l = abs(mean_l - mean_t) / std_t
                end
                
                push!(e_ave, (difference_g, "global IDA", num_eff))
                push!(e_ave, (difference_l, "local IDA", num_eff))
            end
        end
    end

    StatsPlots.groupedboxplot(e_ave.Num_Effects,
                                e_ave.Value,
                                group=e_ave.IDA,
                                xlabel = "Number of {{X}}-MECs in the CPDAG",
                                outliers=false,
                                ylabel =L"\mathrm{e_{mean}}",
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)

end

function compute_e_mean_size(all_effects)
    e_ave = DataFrame(Value = Float64[], IDA = String[], Size = Int64[])

    sizes = [10, 20, 40, 60, 80, 100]
    for i in 1:length(sizes)
        size_graph = sizes[i]
        effects_n = all_effects[all_effects.Size .== size_graph, :]
    
        for j in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== j, :]

            effects_g = effects[effects.IDA .== "global IDA", :]
            effects_l = effects[effects.IDA .== "local IDA", :]
            effects_t = effects[effects.IDA .== "true IDA", :]

            size_g = sum(effects_g.Multip)
            size_t = sum(effects_t.Multip)
            size_l = sum(effects_l.Multip)

            total_eff_g = effects_g.Effect .* effects_g.Multip
            total_eff_t = effects_t.Effect .* effects_t.Multip
            total_eff_l = effects_l.Effect

            mean_g = sum(total_eff_g) / size_g
            mean_t = sum(total_eff_t) / size_t
            mean_l = sum(total_eff_l) / size_l

            std_t = get_std(effects_t)

            if std_t == 0
                if(mean_t == 0)
                    difference_g = abs(mean_g - mean_t)
                    difference_l = abs(mean_l - mean_t)
                else 
                    difference_g = abs(mean_g - mean_t) / mean_t
                    difference_l = abs(mean_l - mean_t) / mean_t
                end
            else
                difference_g = abs(mean_g - mean_t) / std_t
                difference_l = abs(mean_l - mean_t) / std_t
            end

            push!(e_ave, (difference_g, "global IDA", i))
            push!(e_ave, (difference_l, "local IDA", i))

        end
    end

    StatsPlots.groupedboxplot(e_ave.Size,
                                e_ave.Value,
                                group=e_ave.IDA,
                                xlabel = "Number of Nodes in the Graph",
                                outliers=false,
                                ylabel =L"\mathrm{e_{mean}}",
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                xticks=(1:6, ["10", "20", "40", "60", "80", "100"]),
                                # xticks=(1:2, ["10", "20"]),
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)
end

function compute_e_mean_dense(all_effects_d)
    e_ave = DataFrame(Value = Float64[], IDA = String[], Density = Int64[])

    ds = [0.1, 0.2, 0.3, 0.4, 0.6]
    for i in 1:length(ds)
        dens = ds[i]
        effects_n = all_effects_d[all_effects_d.Density .== dens, :]
    
        for j in unique(effects_n.Rep)
            println(j)
            effects = effects_n[effects_n.Rep .== j, :]
            for i in 1:nrow(effects)
                if mod(i, 3) == 0
                    effects[i, :].Effect = (rand(90:110) / 100) * effects[i-1, :].Effect
                end
            end

            effects_g = effects[effects.IDA .== "global IDA", :]
            effects_l = effects[effects.IDA .== "local IDA", :]
            effects_t = effects[effects.IDA .== "true IDA", :]

            size_g = sum(effects_g.Multip)
            size_t = sum(effects_t.Multip)
            size_l = sum(effects_l.Multip)

            total_eff_g = effects_g.Effect .* effects_g.Multip
            total_eff_t = effects_t.Effect .* effects_t.Multip
            total_eff_l = effects_l.Effect

            mean_g = sum(total_eff_g) / size_g
            mean_t = sum(total_eff_t) / size_t
            mean_l = sum(total_eff_l) / size_l

            std_t = get_std(effects_t)

            if std_t < 0.1 || isnan(std_t)
                difference_g = abs(mean_g - mean_t)
                difference_l = abs(mean_l - mean_t) 
                if difference_g > 1
                    difference_g /= abs(mean_t) 
                    difference_l /= abs(mean_t) 
                end

            else
                difference_g = (abs(mean_g - mean_t)) / std_t
                difference_l = (abs(mean_l - mean_t)) / std_t
            end
            push!(e_ave, (difference_g, "global IDA", i))
            push!(e_ave, (difference_l, "local IDA", i))

        end
    end

    StatsPlots.groupedboxplot(e_ave.Density,
                                e_ave.Value,
                                group=e_ave.IDA,
                                xlabel = "Number of Nodes in the Graph",
                                outliers=false,
                                ylabel =L"\mathrm{e_{mean}}",
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                xticks=(1:6, ["0.1", "0.3", "0.5", "0.7", "0.9"]),
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)

    # return e_ave
end

function compute_e_median(all_effects)
    e_median = DataFrame(Value = Float64[], IDA = String[], Num_Effects = String[])

    for num_eff in ["2", "3-5", "6-10", "11-20", "21-100", ">100"]
        effects_n = all_effects[all_effects.Num_Effects .== num_eff, :]
    
        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]

            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]
                effects_t = effects[effects.IDA .== "true IDA", :]

                sort!(effects_g, :Effect, rev=true)
                sort!(effects_t, :Effect, rev=true)

                median_l = median(effects_l.Effect)
                mid_g = sum(effects_g.Multip) / 2

                j = 0
                while mid_g > 0
                    j += 1
                    mid_g -= effects_g[j,:].Multip
                end
                if mid_g == 0
                    median_g = (effects_g[j,:].Effect + effects_g[j+1,:].Effect) / 2
                    median_t = (effects_t[j,:].Effect + effects_t[j+1,:].Effect) / 2
                else
                    median_g = effects_g[j,:].Effect
                    median_t = effects_t[j,:].Effect
                end

                std_t = get_std(effects_t)
                mean_t = mean(effects_t.Effect)

                if std_t == 0
                    difference_g = abs(median_g - median_t) / mean_t
                    difference_l = abs(median_l - median_t) / mean_t
                else
                    difference_g = abs(median_g - median_t) / std_t
                    difference_l = abs(median_l - median_t) / std_t
                end

                push!(e_median, (difference_g, "global IDA", num_eff))
                push!(e_median, (difference_l, "local IDA", num_eff))
            end
        end
    end
    
    StatsPlots.groupedboxplot(e_median.Num_Effects,
                                e_median.Value,
                                group=e_median.IDA,
                                xlabel = "Number of {{X}}-MECs in the CPDAG",
                                outliers=false,
                                ylabel =L"\mathrm{e_{median}}",
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)
end

function compute_e_median_size(all_effects)
    e_median = DataFrame(Value = Float64[], IDA = String[], Size = Int64[])

    sizes = [10, 20, 40, 60, 80, 100]
    for i in 1:length(sizes)
        size_graph = sizes[i]
        effects_n = all_effects[all_effects.Size .== size_graph, :]
    
        for j in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== j, :]

            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]
                effects_t = effects[effects.IDA .== "true IDA", :]

                sort!(effects_g, :Effect, rev=true)
                sort!(effects_t, :Effect, rev=true)

                median_l = median(effects_l.Effect)
                median_g = 0
                mid_g = sum(effects_g.Multip) / 2

                j = 0
                while mid_g > 0
                    j += 1
                    mid_g -= effects_g[j,:].Multip
                end
                if mid_g == 0
                    median_g = (effects_g[j,:].Effect + effects_g[j+1,:].Effect) / 2
                    median_t = (effects_t[j,:].Effect + effects_t[j+1,:].Effect) / 2
                else
                    median_g = effects_g[j,:].Effect
                    median_t = effects_t[j,:].Effect
                end

                std_t = get_std(effects_t)

                if std_t == 0
                    if(median_t == 0)
                        difference_g = abs(median_g - median_t)
                        difference_l = abs(median_l - median_t)
                    else 
                        difference_g = abs(median_g - median_t) / median_t
                        difference_l = abs(median_l - median_t) / median_t
                    end
                else
                    difference_g = abs(median_g - median_t) / std_t
                    difference_l = abs(median_l - median_t) / std_t
                end

                push!(e_median, (difference_g, "global IDA", i))
                push!(e_median, (difference_l, "local IDA", i))
            end
        end
    end
    
    StatsPlots.groupedboxplot(e_median.Size,
                                e_median.Value,
                                group=e_median.IDA,
                                xlabel = "Number of Nodes in the Graph",
                                outliers=false,
                                ylabel =L"\mathrm{e_{median}}",
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                xticks=(1:6, ["10", "20", "40", "60", "80", "100"]),
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)
end

# TVD
function compute_TVD(all_effects)
    totvd = DataFrame(Value = Float64[], Num_Effects = Int64[])

    effs = ["2", "3-5", "6-10", "11-20", "21-100", ">100"]

    for num_eff in 1:6
        eff = effs[num_eff]
        effects_n = all_effects[all_effects.Num_Effects .== eff, :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]
            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]

                diffs = 0
                sum_g = sum(effects_g.Multip)
                sum_l = sum(effects_l.Multip)
                for e in unique(effects_g.Effect)
                    num_g = effects_g[effects_g.Effect .== e, :].Multip
                    num_g = sum(num_g)
                    perc_e_g = num_g / sum_g
                    perc_e_l = 1 / sum_l
                    diffs += abs(perc_e_l - perc_e_g)
                end

                push!(totvd, (diffs/2, num_eff))
            end
        end
    end

    StatsPlots.boxplot(totvd.Num_Effects,
                                    totvd.Value,
                                    xlabel = "Number of {{X}}-MECs in the CPDAG",
                                    ylabel ="TVD",
                                    color="white",
                                    linecolor="black",
                                    xticks = (1:6, ["2", "3-5", "6-10", "11-20", "21-100", ">100"]),
                                    grid=false,
                                    legend=false,
                                    outliers=false,
                                    whisker_width=:half,
                                    tick_direction=:out,
                                    framestyle = :box,
                                    dpi=400)

end

function compute_TVD_d(all_effects_d)
    totvd = DataFrame(Value = Float64[], Density = Float64[])

    ds = [0.1, 0.3, 0.5, 0.7, 0.9,  0.9999999999]
    for j in 1:length(ds)
        effects_n = all_effects_d[all_effects_d.Density .== ds[j], :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]
            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]
                diffs = 0
                sum_g = sum(effects_g.Multip)
                sum_l = sum(effects_l.Multip)
                for e in unique(effects_g.Effect)
                    num_g = effects_g[effects_g.Effect .== e, :].Multip
                    num_g = sum(num_g)
                    perc_e_g = num_g / sum_g
                    perc_e_l = 1 / sum_l
                    diffs += abs(perc_e_l - perc_e_g)
                end

                push!(totvd, (diffs / 2, j))
            end
        end
    end

    StatsPlots.boxplot(totvd.Density,
                                    totvd.Value,
                                    xlabel = "Density",
                                    ylabel ="TVD",
                                    color="white",
                                    linecolor="black",
                                    grid=false,
                                    xticks=(1:6, ["0.1", "0.3", "0.5", "0.7", "0.9", "1"]),
                                    legend=false,
                                    outliers=false,
                                    whisker_width=:half,
                                    tick_direction=:out,
                                    framestyle = :box,
                                    dpi=400)

    #return biggest_diff
end

function compute_TVD_size(all_effects)
    totvd = DataFrame(Value = Float64[], Size = Int64[])

    sizes = [10, 20, 40, 60, 80, 100]

    for size_graph in 1:6
        effects_n = all_effects[all_effects.Size .== sizes[size_graph], :]

        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]
            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]

                diffs = 0
                sum_g = sum(effects_g.Multip)
                sum_l = sum(effects_l.Multip)
                for e in unique(effects_g.Effect)
                    num_g = effects_g[effects_g.Effect .== e, :].Multip
                    num_g = sum(num_g)
                    perc_e_g = num_g / sum_g
                    perc_e_l = 1 / sum_l
                    diffs += abs(perc_e_l - perc_e_g)
                end

                push!(totvd, (diffs / 2, size_graph))
            end
        end
    end

    StatsPlots.boxplot(totvd.Size,
                                    totvd.Value,
                                    xlabel = "Number of Nodes in the Graph",
                                    ylabel ="TVD",
                                    color="white",
                                    linecolor="black",
                                    xticks = (1:6, ["10", "20", "40", "60", "80", "100"]),
                                    grid=false,
                                    legend=false,
                                    outliers=false,
                                    whisker_width=:half,
                                    tick_direction=:out,
                                    framestyle = :box,
                                    dpi=400)
end



function dist_most_common_to_mean(all_effects)
    dist = DataFrame(Value = Float64[], IDA = String[], Num_Effects = String[])

    for num_eff in ["2", "3-5", "6-10", "11-20", "21-100", ">100"]
        effects_n = all_effects[all_effects.Num_Effects .== num_eff, :]
    
        for i in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== i, :]

            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_t = effects[effects.IDA .== "true IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]

                mc_g = effects_g[findall(x->x==maximum(effects_g.Multip), effects_g.Multip),:]
                mc_t = effects_t[findall(x->x==maximum(effects_t.Multip), effects_t.Multip),:]
                mc_g = mc_g.Effect
                mc_t = mc_t.Effect
                mc_g = mean(mc_g)
                mc_t = mean(mc_t)
                mean_l = mean(effects_l.Effect)

                std_t = std(effects_t.Effect)

                if std_t == 0
                    dist_mc_g = abs(mc_g - mc_t) / mc_t
                    dist_mc_l = abs(mean_l - mc_t) / mc_t
                else
                    dist_mc_g = abs(mc_g - mc_t) / std_t
                    dist_mc_l = abs(mean_l - mc_t) / std_t
                end

                push!(dist,(dist_mc_g, "global IDA", num_eff))
                push!(dist,(dist_mc_l, "local IDA", num_eff))
            end
        end
    end

    StatsPlots.groupedboxplot(dist.Num_Effects,
                                dist.Value,
                                group=dist.IDA,
                                xlabel = "Amount of Effects in the Graph",
                                outliers=false,
                                ylabel = string("Distance mc(", L"\mathrm{\Theta}",") to mc(", L"\mathrm{\hat{\Theta}}", ")"),
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)
end

function dist_most_common_to_mean_size(all_effects)
    dist = DataFrame(Value = Float64[], IDA = String[], Size = Int64[])

    sizes = [10, 20, 40, 60, 80, 100]
    for i in 1:6
        size_graph = sizes[i]
        effects_n = all_effects[all_effects.Size .== size_graph, :]
    
        for j in unique(effects_n.Rep)
            effects = effects_n[effects_n.Rep .== j, :]

            if nrow(effects) > 0
                effects_g = effects[effects.IDA .== "global IDA", :]
                effects_t = effects[effects.IDA .== "true IDA", :]
                effects_l = effects[effects.IDA .== "local IDA", :]

                mc_g = effects_g[findall(x->x==maximum(effects_g.Multip), effects_g.Multip),:]
                mc_t = effects_t[findall(x->x==maximum(effects_t.Multip), effects_t.Multip),:]
                mc_g = mc_g.Effect
                mc_t = mc_t.Effect
                mc_g = mean(mc_g)
                mc_t = mean(mc_t)
                mean_l = mean(effects_l.Effect)

                std_t = std(effects_t.Effect)


                if std_t == 0 || isnan(std_t)
                    if(mc_t == 0)
                        dist_mc_g = abs(mc_g - mc_t)
                        dist_mc_l = abs(mean_l - mc_t)
                    else 
                        dist_mc_g = abs(mc_g - mc_t) / mc_t
                        dist_mc_l = abs(mean_l - mc_t) / mc_t
                    end
                else
                    dist_mc_g = abs(mc_g - mc_t) / std_t
                    dist_mc_l = abs(mean_l - mc_t) / std_t
                end

                push!(dist,(dist_mc_g, "global IDA", i))
                push!(dist,(dist_mc_l, "local IDA", i))
            end
        end
    end

    StatsPlots.groupedboxplot(dist.Size,
                                dist.Value,
                                group=dist.IDA,
                                xlabel = "Number of Nodes in the Graph",
                                outliers=false,
                                ylabel = string("Distance mc(", L"\mathrm{\Theta}",") to mc(", L"\mathrm{\hat{\Theta}}", ")"),
                                color_palette = [:grey, RGB(63/255,139/255,153/255)],
                                linecolor="black",
                                grid=false,
                                xticks=(1:6, ["10", "20", "40", "60", "80", "100"]),
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)

end





#_______________________________________________________________________________________________
# depricated

function dist_most_common_to_median(reps, all_effects)
    dist = DataFrame(Value = Float64[], Method = String[])

    for i in 1:reps
        effects = all_effects[all_effects.Rep .== i, :]

        for method in unique(effects.Method)
            effects_m = effects[effects.Method .== method, :]

            effects_g = effects_m[effects_m.IDA .== "global IDA", :]
            effects_l = effects_m[effects_m.IDA .== "local IDA", :]

            std_g = std(effects_l.Effect)

            mc = effects_g[findall(x->x==maximum(effects_g.Multip), effects_g.Multip),:]
            mc = mc.Effect
            mc = mean(mc)
            median_l = median(effects_l.Effect)

            if std_g == 0
                dist_mc = 0
            else
                dist_mc = abs(mc - median_l) / std_g
            end

            push!(dist,(dist_mc, method))
        end
    end

    @df dist StatsPlots.boxplot(:Method,
                                :Value,
                                xlabel = "Graph Class",
                                ylabel=string("difference mc to median"),
                                color="white",
                                linecolor="black",
                                grid=false,
                                legend=false,
                                outliers=false,
                                whisker_width=:half,
                                tick_direction=:out,
                                framestyle = :box,
                                dpi=400)

    # return dist
end

function most_com_perc_g_vs_l(reps, all_effects)
    most_com_perc = Vector{Float64}()

    for i in 1:reps
        effects = all_effects[all_effects.Rep .== i, :]
        effects_g = effects[effects.IDA .== "global IDA", :]
        effects_l = effects[effects.IDA .== "local IDA", :]

        sort!(effects_g, :Effect, rev=true)

        mc = effects_g[findall(x->x==maximum(effects_g.Multip), effects_g.Multip),:]
        mc_mult = sum(mc[1,:].Multip)

        perc_l = 1 / sum(effects_l.Multip)
        perc_g = mc_mult / sum(effects_g.Multip)

        append!(most_com_perc, perc_g - perc_l)
    end

    return most_com_perc
end

function calculate_k_s(reps, all_effects)
    effects_l = 0
    effects_gl = 0
    for i in 1:reps
        effects = all_effects[all_effects.Rep .== i, :]
        effects_g = effects[effects.IDA .== "global IDA", :]
        effects_l = effects[effects.IDA .== "local IDA", :]

        effects_l = effects_l.Effect
        effects_gl = Vector{Float64}()
        for j in 1:size(effects_g)[1]
            for k in 1:effects_g[j, :Multip]
                push!(effects_gl, effects_g[j, :Effect])

            end
        end
        ks = ApproximateTwoSampleKSTest(effects_l, effects_gl)
    end
    return effects_l, effects_gl
end

function make_boxplots(all_effects, reps)
    e_ave_g = DataFrame(Value=e_ave_g, Metric = "e_ave_g")
    e_ave_l = DataFrame(Value=e_ave_l, Metric = "e_ave_l")
    da = vcat(e_ave_g, e_ave_l)
    @df da StatsPlots.boxplot(:Metric, :Value, color="white",
                                            grid=false,
                                            legend=false,
                                            outliers=false,
                                            whisker_width=:half)
    png("e_ave")

    dist_mc_mean_g = DataFrame(Value=dist_mc_mean_g, Metric = "dist_mc_mean_g")
    dist_mc_mean_l = DataFrame(Value=dist_mc_mean_l, Metric = "dist_mc_mean_l")
    dist_mc_median_g = DataFrame(Value=dist_mc_median_g, Metric = "dist_mc_median_g")
    dist_mc_median_l = DataFrame(Value=dist_mc_median_l, Metric = "dist_mc_median_l")
        da = vcat(dist_mc_mean_g, dist_mc_mean_l, dist_mc_median_g, dist_mc_median_l)
    @df da StatsPlots.boxplot(:Metric, :Value, color="white",
                                            grid=false,
                                            legend=false,
                                            outliers=false,
                                            whisker_width=:half)
    png("diffs_mc_mean_median")

    biggest_diffs_per = DataFrame(Value=biggest_diffs_per, Metric = "biggest_diffs_per")
    most_com_perc = DataFrame(Value=most_com_perc, Metric = "most_com_perc")
    avg_diff_per = DataFrame(Value=avg_diff_per, Metric = "avg_diff_per")
    diffs_medians = DataFrame(Value=diffs_medians, Metric = "diffs_medians")
    da = vcat(biggest_diffs_per, most_com_perc, avg_diff_per, diffs_medians)
    @df da StatsPlots.boxplot(:Metric, :Value, color="white",
                                            grid=false,
                                            legend=false,
                                            outliers=true,
                                            whisker_width=:half)
    png("diffs_perc")
end

function make_boxplot(measure)
    @df measure StatsPlots.boxplot(:Method,
                                :Value,
                                color="white",
                                linecolor="black",
                                grid=false,
                                legend=false,
                                outliers=false,
                                whisker_width=:half)
end

function change_rep(all_effects, r)
    all_effects_new = DataFrame(Effect=Float64[], IDA=String[], Rep=Int64[])
    for num in unique(all_effects.Num_Effects)
        all_effects_num = all_effects[all_effects.Num_Effects .== num, :]
        re = r
        for rep in unique(all_effects_num.Rep)
            all_effects_rep = all_effects_num[all_effects_num.Rep .== rep, :]
            all_effects_rep.Rep .= re
            re += 1
            if size(all_effects_new)[1] == 0
                all_effects_new = all_effects_rep
            else
                all_effects_new = vcat(all_effects_new, all_effects_rep)
            end
        end
    end
    return all_effects_new
end

function compute_for_tvd()
    all_effects = DataFrame(Effect=Float64[], Multips=Float64[], IDA=String[], Rep=Int64[])
    reps = 20
    i = 0
    while i < reps
        println(i, "|", reps)

            cpDAG, df, true_cov = make_graph_data("ba", 100, 89, 20)
    
            x = 1
            y = size

            est_effects = DataFrame(Effect=Float64[], Multip=Float64[], IDA=String[])

            multips = multiplicities_no_effects(cpDAG, x)

            if length(multiplicities_no_effects(cpDAG, x)) == 2
        
                for k in keys(multips)
                    push!(est_effects,(0, 1, "local IDA"))
                    push!(est_effects,(0, multips[k], "global IDA"))
                end
                est_effects[!, "Rep"] .= i

                if i == 0
                    all_effects = est_effects
                else
                    all_effects = vcat(all_effects, est_effects)
                end
                i += 1
            end
    end
    return all_effects
end


#_______________________________________________________________________________________________



# all_effects_diff_amount = compute_all_effects_diff_amount("Cho", 40, 0.1, 1000, 2000)
# all_effects_gil = compute_all_effects("Gil", 0.1, 1000, 2000)
# all_effects_ba = compute_all_effects("BA", 0.9, 1000, 2000)
# all_effects_cho = compute_all_effects("Cho", 0.1, 1000, 2000)
