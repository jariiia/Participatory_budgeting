using JuMP
import HiGHS
import MultiObjectiveAlgorithms as MOA
import Plots

input_file = string("data/",ARGS[1],".txt")
results_file = string("results/PB_", ARGS[1],".txt")
plot_file = string("plots/PB_", ARGS[1],".png")

open(input_file) do f

 	# line_number
  	line = 0 

	# read till end of file
        while ! eof(f)  
        
            # read a new / next line for every iteration           
            s = readline(f)
	    line += 1
	    token = split(s)[3]
	    if line == 1
	        global citsat = parse.(Float64, split(chop(token; head=1, tail=1), ','))
	    elseif line == 2
	    	global alignment = parse.(Float64, split(chop(token; head=1, tail=1), ','))
	    elseif line == 3
	    	global cost = parse.(Float64, split(chop(token; head=1, tail=1), ','))
	    else
		global budget = parse(Float64,token)
	    end
         end
end


N = length(citsat)

model = Model()
@variable(model, x[1:N], Bin)
@constraint(model, sum(cost[i] * x[i] for i in 1:N) <= budget)
@expression(model, citsat_expr, sum(citsat[i] * x[i] for i in 1:N))
@expression(model, alignment_expr, sum(alignment[i] * x[i] for i in 1:N))
@objective(model, Max, [citsat_expr, alignment_expr])

set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_silent(model)

# set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
set_attribute(model, MOA.Algorithm(), MOA.TambyVanderpooten())

optimize!(model)
 # solution_summary(model)

number_of_solutions = result_count(model)


# Compute selected proposals per solution

selected_proposals = Vector{Vector{Int}}(undef, number_of_solutions)
for i in 1:number_of_solutions
    selected_proposals[i] = [j for j in 1:N if value(x[j]; result = i) > 0]
end


# Compute Nash and Kalai-Smordinsky solutions

max_sat = maximum([value(citsat_expr; result = i) for i in 1:number_of_solutions])
max_al = maximum([value(alignment_expr; result = i) for i in 1:number_of_solutions])
min_al = minimum([value(alignment_expr; result = i) for i in 1:number_of_solutions])
den_ks = sqrt(max_sat * max_sat + max_al * max_al)

nash_products = Vector{Float64}(undef, number_of_solutions)
ks_ratios = Vector{Float64}(undef, number_of_solutions)
for i in 1:number_of_solutions
    nash_products[i] = value(citsat_expr; result = i) * value(alignment_expr; result = i)
    ks_ratios[i] = abs(max_al*value(citsat_expr; result = i) - max_sat*value(alignment_expr; result = i))/den_ks
end
nash_solution = argmax(nash_products)
ks_solution = argmin(ks_ratios)



 # Compute Szymkiewicz–Simpson coefficient

nash_vs_classic = size(intersect(selected_proposals[nash_solution],selected_proposals[number_of_solutions]))[1]/min(size(selected_proposals[nash_solution])[1],size(selected_proposals[number_of_solutions])[1])

ks_vs_classic = size(intersect(selected_proposals[ks_solution],selected_proposals[number_of_solutions]))[1]/min(size(selected_proposals[ks_solution])[1],size(selected_proposals[number_of_solutions])[1])


# Log results



open(results_file,"w") do io
     println(io,"# Pareto front points: solution index, citizen satisfaction, alignment")
     for i in 1:number_of_solutions
     	 println(io, i,",", value(citsat_expr; result = i),",",value(alignment_expr; result = i))
     end
     println(io,"# Nash solution index")
     println(io, nash_solution)
     println(io,"# Kalai-Smordinsky solution index")
     println(io, ks_solution)
     println(io,"# Szymkiewicz–Simpson coefficients")
     println(io,"# Nash vs citizen satisfaction")
     println(io,nash_vs_classic)
     println(io,"# Kalai-Smordinsky vs citizen satisfaction")
     println(io,ks_vs_classic)
     println(io,"# solution index , list of selected proposals")
     for i in 1:number_of_solutions
	 println(io, i,",",join(selected_proposals[i],","))
     end
end

# Plot results

plot = Plots.scatter(
    [value(citsat_expr; result = i) for i in 1:number_of_solutions],
    [value(alignment_expr; result = i) for i in 1:number_of_solutions];
    xlabel = "Citizen satisfaction",
    ylabel = "Value alignment",
    title = "Objective space",
    label = "",
#    xlims = (18000, 21500),
)
for i in 1:number_of_solutions
    y = objective_value(model; result = i)
    Plots.annotate!(y[1] - 1, y[2], (i, 10))
end
ideal_point = objective_bound(model)
Plots.scatter!([ideal_point[1]], [ideal_point[2]]; label = "Ideal point")

Plots.savefig(plot,plot_file)