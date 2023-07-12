
using CSV
using DataFrames
using CausalInference
using GraphPlot

include("ida.jl")





function get_effects_all_alphas(df, alphas, x, y)
    effects_all_alphas = DataFrame(Alpha = String[], Effect=Float64[], Multip=Float64[])
    for alpha in alphas
        est_CPDAG = pcalg(df, alpha, gausscitest)
        est_effects, multips = multiplicities(est_CPDAG, x, y, df)
        effects_all_alphas = add_all_effects!(effects_all_alphas, est_effects, multips, alpha)
    
    end
    return effects_all_alphas
end



function add_all_effects_global(est_effects, multips)

    effs = DataFrame(Effect=Float64[], Probability=Float64[])

    for k in keys(est_effects)
        if get(est_effects, k, 0) in effs[:,1]
        n = get(est_effects, k, 0)
        effs[effs.Effect .== n, :Probability] .= effs[effs.Effect .== n, :Probability] .+ get(multips, k, 0)/sum(values(multips))
        else
            push!(effs, (get(est_effects, k, 0), get(multips, k, 0)/sum(values(multips))))
        end
    end

    return effs
end

function add_all_effects_local(est_effects)

    effs = DataFrame(Effect=Float64[], Probability=Float64[])

    for k in keys(est_effects)
        if get(est_effects, k, 0) in effs[:,1]
        n = get(est_effects, k, 0)
        effs[effs.Effect .== n, :Probability] .= effs[effs.Effect .== n, :Probability] .+ 1/length(est_effects)
        else
            push!(effs, (get(est_effects, k, 0), 1/length(est_effects)))
        end
    end

    return effs
end


function gaia_demo_global(csv, x, y, alpha_val)
    data = DataFrame(csv)
    est_CPDAG = pcalg(data, alpha_val, gausscitest)
    est_effects, multips = multiplicities(est_CPDAG, x, y, data)
    all_effects = add_all_effects_global(est_effects, multips)
    CSV.write("all_effects.csv", all_effects)
    return all_effects
end

function gaia_demo_local(csv, x, y, alpha_val)
    data = DataFrame(csv)
    est_CPDAG = pcalg(data, alpha_val, gausscitest)
    est_effects = loc_ida(est_CPDAG, x, y, data)
    all_effects = add_all_effects_local(est_effects)
    CSV.write("all_effects.csv", all_effects)
    return all_effects
end


# effects_all_alphas = get_effects_all_alphas(syntrea_data, alphas, x, y)



csv = CSV.read("/Users/bjornfilter/Desktop/GAIA/csv/syndata.csv", DataFrame)
x = 7
y = 16
alpha_val = 0.02
all_effects = gaia_demo_local(csv, x, y, alpha_val)


# nodelabel = 1:nv(net)
# gplot(net, nodelabel=nodelabel)