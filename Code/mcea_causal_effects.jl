using CSV
using DataFrames
using CausalInference
using GraphPlot

include("ida.jl")
include("plotting_fun.jl")


alphas = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
x = 1
y = 8




function add_all_effects!(effects_all_alphas, est_effects, multips, alpha)

    effs = DataFrame(Alpha = String[], Effect=Float64[], Multip=Float64[])

    for k in keys(est_effects)
        effect = get(est_effects, k, 0)
        if effect in effs[:, :Effect]
            effs[!, :Multip][effs[!, :Effect] .== effect] .= [get(multips, k, 0)/sum(values(multips))] + effs[!, :Multip][effs[!, :Effect] .== effect]
        else
            push!(effs, (string(alpha), get(est_effects, k, 0), get(multips, k, 0)/sum(values(multips))))
        end
    end

    for i in 1:size(effs)[1]
        push!(effects_all_alphas, effs[i, :])
    end

    return effects_all_alphas
end

function get_effects_all_alphas(df, alphas, x, y, )
    effects_all_alphas = DataFrame(Alpha = String[], Effect=Float64[], Multip=Float64[])
    for alpha in alphas
        est_CPDAG = pcalg(df, alpha, gausscitest)
        est_effects, multips = multiplicities(est_CPDAG, x, y, df)
        effects_all_alphas = add_all_effects!(effects_all_alphas, est_effects, multips, alpha)
    
    end
    return effects_all_alphas
end




cd("/Users/bjornfilter/Desktop/MCEA_Data")
df = DataFrame(CSV.File("data-prepared-julia.csv", header=1, delim=";", decimal=','))

effects_all_alphas = get_effects_all_alphas(df, alphas, x, y, )
bar_plot_diff_efffects(effects_all_alphas)








# nodelabel = 1:nv(est_CPDAG)
# gplot(est_CPDAG, nodelabel=nodelabel)