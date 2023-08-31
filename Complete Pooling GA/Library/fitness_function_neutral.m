function error = fitness_function_neutral(params_array)
    debug_flag = false;

    % 25 parameters
    params_list = ["td", "t_delay", "t_last", "t_delay12", "t_last12", ...
        "k1", "k2", "k3", "k4", "k5", "k6", ...
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", ...
        "s1", "s2", ...
        "ig0", "c0", "p0", "vl0"];

    load Data/preprocessed_data data

    for i = 1:numel(params_list)
        p.(params_list(i)) = params_array(i);
    end
    
    treatments_available = methods(treatment_factory);
    treatment_processed = 1;
    error_per_treatment = [];
    error_vector = [];

    % Iterate over all treatments
    for treatment = treatments_available(1:end-1).' % Need to remove last method bc it's not a real methods (it's a callback)
        treatment_spec = treatment_factory.(treatment{1})();

        if debug_flag 
            plot_info.flag = true;
            plot_info.treatment_name = treatment{1};
            plot_info.subplot_id = treatment_processed;
        else
            plot_info.flag = false;
        end

        sol = immuno_solver(p, treatment_spec, plot_info);
        total_tumour_vol = sol.y(4,:) + sol.y(5,:);

        % Calculate MSE between experiment data and solution for each
        % treatment and take average of average for overall error
        i = 1;

        error_per_day = [];
        for day = data.(treatment{1}).valid_days
            in_vivo = data.(treatment{1}).vector_avg(i);
            [~, sol_index] = min(abs(day .* ones(1, numel(sol.x)) - sol.x));
            if debug_flag
                figure(1); subplot(2, 3, treatment_processed); hold on
                scatter(sol.x(sol_index), total_tumour_vol(sol_index), "b", "HandleVisibility","off")
            end
            in_silico = total_tumour_vol(sol_index);
            error_per_day(end+1) = ((in_silico - in_vivo)^2)/1000;
            error_vector(end+1) = (in_silico - in_vivo)^2;
            i = i + 1;
        end

    error_per_treatment(treatment_processed) = sum(error_per_day);
    treatment_processed = treatment_processed + 1;
    end

% Combine error from all 5 treatments  
% error = mean(error_per_treatment.^2);
% error = error_per_treatment(1);
% error = -sum(error_vector);
error = sum(error_per_treatment(1:3)) - sum(error_per_treatment(4:6));
disp(error)
% disp(error_per_treatment)
end

