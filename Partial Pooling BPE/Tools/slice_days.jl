"""
Auxiliary function: from a numerical solution, only selected the
relant data points
# Arguments
- step: step size used for the numerical method. Must be a divisor of 1 
- selected days: days to extrac from the numerical approximation
- num_sol: the numerical solution
"""
function slice_days(step, selected_days, num_sol)
    extracted_days = Vector{Float64}(undef, length(selected_days))
    for (index, day) in enumerate(selected_days)
        
    end
end