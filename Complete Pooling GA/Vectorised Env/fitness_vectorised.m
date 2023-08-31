function error = fitness_vectorised(params_matrix)
    % along rows: population, 
    % along cols: parameters (check order!)

    load Data/preprocessed_data data % Syntax must be compatible with HPC cluster

    treatments_available = methods(treatment_factory);
    treatment_processed = 1;

    % Iterate over all treatments
    error_per_treatment = [];
    for treatment = treatments_available(1:end-1).' % Need to remove last method bc it's not a real methods (it's a callback)
        treatment_spec = treatment_factory.(treatment{1})();
        sol = immuno_solver(params_matrix, treatment_spec);

        % Calculate MSE between experiment data and solution for each
        % treatment and rake average of average for overall error
        i = 1;
        error_per_day = [];
        for day = data.(treatment{1}).valid_days
            in_vivo = data.(treatment{1}).vector_avg(i);
            [~, sol_index] = min(abs(day .* ones(1, numel(sol.x)) - sol.x));
            total_tumour_vol = sol.y(4,:) + sol.y(5,:);
            in_silico = total_tumour_vol(sol_index);
            error_per_day(end+1) = ((in_silico - in_vivo)^2)/1000;
            i = i + 1;
        end
    error_per_treatment(treatment_processed) = mean(error_per_day);
    treatment_processed = treatment_processed + 1;
    end
error = mean(error_per_treatment);
disp(error)
end

