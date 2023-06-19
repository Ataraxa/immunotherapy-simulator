function error = fitness_function(params_array)
    % 21 parameters
    params_list = ["td", "t_delay", "t_last", "t_delay12", "t_last12", ...
        "k1", "k2", "k3", "k4", "k5", "k6", ...
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", ...
        "s1", "s2"];

    load Data\preprocessed_data.mat data
    % Boring part. Thanks matlab. surely there's a better way to do that
    for i = 1:numel(params_list)
        p.(params_list(i)) = params_array(i);
    end
    error = 0;
    treatments_available = methods(treatment_factory);
    error_per_treatment = 0;
    % Iterate over all treatments
    for treatment = treatments_available(1:end-1).' % Need to remove last method bc it's not a real methods (it's a callback)
        treatment_spec = treatment_factory.(treatment{1})();
        sol = immuno_solver(p, treatment_spec);

        % Calculate MSE between experiment data and solution
        i = 1;
        for day = data.(treatment{1}).valid_days
            in_vivo = data.(treatment{1}).vector_avg(i);
            [~, sol_index] = min(abs(day .* ones(1, numel(sol.x)) - sol.x));
            in_silico = sol.y(sol_index);
%             error = error + (in_silico - in_vivo)^2;
            error_per_treatment(end+1) = (in_silico - in_vivo)^2;
            i = i + 1;
        end
    end
    error = mean(error_per_treatment);
    disp(error)
end

