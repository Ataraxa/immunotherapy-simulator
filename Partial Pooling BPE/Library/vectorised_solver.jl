
"Function that accepts all the parameters as a vector and solve the equation"
function sensitivity_solver(params_vector)
    # Assume vector = [p; vl0]
    p = params_vector[1:end-1] 
    u0 = [0.0084; 9.56; 4.95; params_vector[end]; 0]
    prob1 = remake(problem;p=p, u0=u0)
    sol = solve(prob1;saveat=0.1)
    return(sol[4,end] + sol[5,end])
end