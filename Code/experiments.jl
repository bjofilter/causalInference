using Statistics
include("activelearning.jl")
include("ida.jl")
include("genData.jl")
include("causal_functions.jl")

function IDAexperiments()
    num_nodes = [10, 20, 40, 60, 80, 100]
    for n in num_nodes
        ltm = []
        gtm = []
        println("n = ", n, ": ")
        rep = 1
        while length(gtm) < 1000
            println(string(n) * " " * string(rep))
            cpDAG, df, true_cov = make_graph_data("Gil", n, 0.1, 2000)
            cpDAG, df, true_cov = make_graph_data("BA", n, 0.9, 2000)
            cpDAG, df, true_cov = make_graph_data("Cho", n, 0.1, 2000)
            x = 1
            y = n

            err = false
            try 
                effects, mults = multiplicities(cpDAG, x, y, df)
                # multips = multiplicities_no_effects(cpDAG, x)
                # multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)
                # effects, mults = glob_ida_mults(x, y, df, multips_grouped)
            catch e
                err = true
            end

            if !err
                rep += 1
                # don't include init time
                effects, mults = multiplicities(cpDAG, x, y, df)
                effects = loc_ida(cpDAG, x, y, df)
                # multips = multiplicities_no_effects(cpDAG, x)
                # multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)
                # effects, mults = glob_ida_mults(x, y, df, multips_grouped)

                if rand(1:2) == 1
                    tmp = @elapsed begin
                        effects, mults = multiplicities(cpDAG, x, y, df)
                        # multips = multiplicities_no_effects(cpDAG, x)
                        # multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)
                        # effects, mults = glob_ida_mults(x, y, df, multips_grouped)
                    end
                    push!(gtm, tmp)
                    
                    # multips = multiplicities_no_effects(cpDAG, x)
                    tmp = @elapsed begin
                        effects = loc_ida(cpDAG, x, y, df)
                        # multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)
                        # effects, mults = loc_ida_mults(x, y, df, multips_grouped)
                    end
                    push!(ltm, tmp)
                else
                    # multips = multiplicities_no_effects(cpDAG, x)
                    tmp = @elapsed begin
                        effects = loc_ida(cpDAG, x, y, df)
                        # multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)
                        # effects, mults = loc_ida_mults(x, y, df, multips_grouped)
                    end
                    push!(ltm, tmp)
                    
                    tmp = @elapsed begin
                        effects, mults = multiplicities(cpDAG, x, y, df)
                        # multips = multiplicities_no_effects(cpDAG, x)
                        # multips_grouped = get_true_diff_parentSets(multips, cpDAG, x, y)
                        # effects, mults = glob_ida_mults(x, y, df, multips_grouped)
                    end
                    push!(gtm, tmp)
                end
            end
        end
        println("  [global] avg time in s: " * string(mean(gtm)))
        println("  [global] std dev for time in s: " * string(std(gtm)))
        println("  [local] avg time in s: " * string(mean(ltm)))
        println("  [local] std dev for time in s: " * string(std(ltm)))
    end
end



##### IDA experiments:
IDAexperiments()