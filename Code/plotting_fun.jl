using Plots
using Cairo, Fontconfig
using StatsPlots
# using CairoMakie
using LaTeXStrings
using Formatting


function plot_true_path_vs_reg(effects_for_dense)

    # max_val = maximum(est_effects[!, "Effect"])
    # width = max_val/250

    @df effects_for_dense StatsPlots.density(:Effect,
        # weights = :multip,
        group = :IDA,
        color_palette = [:grey, RGB(63/255,139/255,153/255)],
        linewidth=3,
        xlabel = "Estimated Coefficient",
        ylabel = "Density",
        grid=false,
        tick_direction=:out,
        framestyle = :box,
        dpi=400)
        # savefig("dense_effects.pdf")  

end


function bar_plot_effects(est_effects)

    global_df = est_effects[est_effects.IDA .== "global IDA", :]
    # global_df.IDA .= "alternative global IDA"
    local_df = est_effects[est_effects.IDA .== "local IDA", :]
    # local_df.IDA .= "alternative local IDA"
    sum_global = sum(global_df.Multip)
    sum_local = sum(local_df.Multip)
    global_df.Multip  = global_df.Multip ./ sum_global
    local_df.Multip  = local_df.Multip ./ sum_local

    est_effects = vcat(global_df, local_df)

    max_val = maximum(est_effects[!, "Effect"])
    width = max_val/25

    StatsPlots.groupedbar(est_effects.Effect,
            est_effects.Multip,
            group = est_effects.IDA,
            color_palette = [:grey, RGB(63/255,139/255,153/255)],
            tick_direction=:out,
            framestyle = :box,
            legend=:topright,
            bar_width = width, 
            xlabel = "Estimated Coefficient",
            ylabel = "Probability",
            grid=false)
    savefig("bar.pdf")
end

function plot_ec(est_effects)
    effects_g = est_effects[est_effects.IDA .== "global IDA", :]
    effects_l = est_effects[est_effects.IDA .== "local IDA", :]

    effects_l = effects_l.Effect
    effects_gl = Vector{Float64}()
    min = minimum(effects_g.multip)
    while min > 1000
        effects_g.multip = effects_g.multip ./ sqrt(min)
        min = minimum(effects_g.multip)
    end
    for j in 1:size(effects_g)[1]
        for k in 1:effects_g[j, :multip]
            push!(effects_gl, effects_g[j, :Effect])

        end
    end

    gcdf_l = StatsPlots.ecdf(effects_l)
    gcdf_g = StatsPlots.ecdf(effects_gl)

    Plots.plot(x -> gcdf_g(x),
                0,
                color = :black,
                maximum(effects_l)*1.1,
                xlabel = "ECDF of Estimated Coefficient",
                ylabel = "Probability",
                tick_direction=:out,
                framestyle = :box,
                linewidth = 3,
                grid=false,
                legend = :bottomright,
                labels="global IDA")
    Plots.plot!(x -> gcdf_l(x),
                0,
                color = RGB(63/255,139/255,153/255),
                maximum(effects_l) *1.1,
                xlabel = "ECDF of Estimated Coefficient",
                ylabel = "Probability",
                tick_direction=:out,
                framestyle = :box,
                linewidth = 3,
                grid=false,
                legend = :bottomright,
                labels="local IDA")
    # savefig("ecdf.pdf")
end

function bar_plot_diff_efffects(effects_all_alphas)

    max_val = maximum(effects_all_alphas[!, "Effect"])
    width = .00001

    StatsPlots.groupedbar(effects_all_alphas.Effect,
    effects_all_alphas.Multip,
    group = effects_all_alphas.Alpha,
    # color_palette = [:grey, RGB(63/255,139/255,153/255)],
    tick_direction=:out,
    framestyle = :box,
    legend=:topright,
    #bar_width = width, 
    xlabel = "Estimated Coefficient",
    ylabel = "Probability",
    grid=false)

    # savefig("bar.pdf")
end