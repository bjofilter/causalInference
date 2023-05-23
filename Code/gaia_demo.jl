
using CSV
using DataFrames
using CausalInference
using GraphPlot

include("ida.jl")





# function get_effects_all_alphas(df, alphas, x, y)
#     effects_all_alphas = DataFrame(Alpha = String[], Effect=Float64[], Multip=Float64[])
#     for alpha in alphas
#         est_CPDAG = pcalg(df, alpha, gausscitest)
#         est_effects, multips = multiplicities(est_CPDAG, x, y, df)
#         effects_all_alphas = add_all_effects!(effects_all_alphas, est_effects, multips, alpha)
    
#     end
#     return effects_all_alphas
# end



function add_all_effects!(est_effects, multips, alpha)

    effs = DataFrame(Effect=Float64[], Multip=Float64[])

    for k in keys(est_effects)
        effect = get(est_effects, k, 0)
            push!(effs, (get(est_effects, k, 0), get(multips, k, 0)/sum(values(multips))))
    end

    return effs
end


function gaia_demo(csv, x, y, alpha_val)
    data = DataFrame(csv)
    est_CPDAG = pcalg(data, alpha_val, gausscitest)
    est_effects, multips = multiplicities(est_CPDAG, x, y, data)
    all_effects = add_all_effects!(est_effects, multips, alpha)
    CSV.write("all_effects.csv", df)
    return all_effects
end


# effects_all_alphas = get_effects_all_alphas(syntrea_data, alphas, x, y)


# nodelabel = 1:nv(net)
# gplot(net, nodelabel=nodelabel)