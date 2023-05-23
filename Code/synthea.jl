using CSV
using DataFrames
using CausalInference
using GraphPlot

include("ida.jl")
include("plotting_fun.jl")


alphas = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
x = 5
y = 4


cd("/Users/bjornfilter/Desktop/GAIA/csv")
patients = DataFrame(CSV.File("patients.csv", header=1, delim=",", decimal='.'))
observations = DataFrame(CSV.File("observations.csv", header=1, delim=",", decimal='.'))

patients = patients[:, [:Id, :BIRTHDATE, :DEATHDATE, :GENDER, :HEALTHCARE_EXPENSES, :HEALTHCARE_COVERAGE, :INCOME]]
observations = observations[:, [:PATIENT, :ENCOUNTER, :CATEGORY, :DESCRIPTION, :VALUE, :UNITS, :TYPE]]
rename!(observations, :PATIENT => :Id)

desc = unique(observations[:, :DESCRIPTION])
desc = collect(skipmissing(desc))
desc = desc[[1; collect(3:14); collect(45:57); collect(95:102); collect(104:111); collect(115:130); collect(133:137); collect(139:151); collect(172:178); 203; 204; collect(231:240); collect(251:225); 276; collect(278:280)] ]
desc_blood = desc[[collect(5:13); collect(15:22)]]

desc = desc_blood

obs = unique(observations[:, :ENCOUNTER])
obs = collect(skipmissing(obs))
num_obs = size(obs)[1]

df = DataFrame([Vector{Union{String, Missing}}(missing, num_obs) for _ = desc] , desc)
df[!, :Encounter] = obs

for i in 1:size(observations)[1]
    if mod(i, 1000) == 0
        println(i, '/', size(observations)[1])
    end
    enc = observations[i, :].ENCOUNTER
    if !ismissing(enc)
        descript =  observations[i, :].DESCRIPTION
        if descript in desc
            val = observations[i, :].VALUE
            df[df.Encounter .== enc, descript] .= val
        end
    end
end

df = dropmissing(df)
select!(df, Not(:"Encounter"))

function toInt(df)
    df[:, "Tobacco smoking status"] = [item == "Never smoked tobacco (finding)" ? "0" : "1" for item in df[:, "Tobacco smoking status"]]
    mapcols!(change_data, df)
end

function change_data(x)
    x = parse.(Float64, x)
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

function get_effects_all_alphas(df, alphas, x, y)
    effects_all_alphas = DataFrame(Alpha = String[], Effect=Float64[], Multip=Float64[])
    for alpha in alphas
        est_CPDAG = pcalg(df, alpha, gausscitest)
        est_effects, multips = multiplicities(est_CPDAG, x, y, df)
        effects_all_alphas = add_all_effects!(effects_all_alphas, est_effects, multips, alpha)
    
    end
    return effects_all_alphas
end



df = toInt(df)

effects_all_alphas = get_effects_all_alphas(df, alphas, x, y)
# bar_plot_diff_efffects(effects_all_alphas)




# nodelabel = 1:nv(est_CPDAG)
# gplot(est_CPDAG, nodelabel=nodelabel)