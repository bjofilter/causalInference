
using CSV
using DataFrames
using CausalInference
using GraphPlot

include("ida.jl")


alphas = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
x = 1
y = 100


cd("/Users/bjornfilter/Desktop/GAIA/csv")
syntrea_data = CSV.read("syndata.csv", DataFrame)
syntrea_net = CSV.read("synnet.csv", DataFrame)

net = SimpleDiGraph(Matrix(syntrea_net))


function get_effects_all_alphas(df, alphas, x, y)
    effects_all_alphas = DataFrame(Alpha = String[], Effect=Float64[], Multip=Float64[])
    for alpha in alphas
        est_CPDAG = pcalg(df, alpha, gausscitest)
        est_effects, multips = multiplicities(est_CPDAG, x, y, df)
        effects_all_alphas = add_all_effects!(effects_all_alphas, est_effects, multips, alpha)
    
    end
    return effects_all_alphas
end



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


effects_all_alphas = get_effects_all_alphas(syntrea_data, alphas, x, y)


# nodelabel = 1:nv(net)
# gplot(net, nodelabel=nodelabel)