function error = fitness_function(params_array)

    % 21 parameters
    params_list = ["td", "t_delay", "t_last", "t_delay12", "t_last12", ...
        "k1", "k2", "k3", "k4", "k5", "k6", ...
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", ...
        "s1", "s2"];

    load Data/preprocessed_data data

    for i = 1:numel(params_list)
        p.(params_list(i)) = params_array(i);
    end
%     disp(p)
    
    treatments_available = methods(treatment_factory);
    treatment_processed = 1;
    error_per_treatment = [];



    % Iterate over all treatments
    for treatment = treatments_available(1:end-1).' % Need to remove last method bc it's not a real methods (it's a callback)
        treatment_spec = treatment_factory.(treatment{1})();
        plot_info.flag = false; plot_info.treatment_name = treatment{1};
        plot_info.subplot_id = treatment_processed;
        sol = immuno_solver(p, treatment_spec, plot_info);
        % Calculate MSE between experiment data and solution for each
        % treatment and rake average of average for overall error
        i = 1;
        error_per_day = [];
        for day = data.(treatment{1}).valid_days
            in_vivo = data.(treatment{1}).vector_avg(i);
            [~, sol_index] = min(abs(day .* ones(1, numel(sol.x)) - sol.x));
            figure(1); subplot(2, 3, treatment_processed); hold on
            total_tumour_vol = sol.y(4,:) + sol.y(5,:);
            scatter(sol.x(sol_index), total_tumour_vol(sol_index), "b", "HandleVisibility","off")
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

