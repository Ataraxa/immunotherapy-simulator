function error = fitness_function(...
    td, t_delay, t_last, t_delay12, t_last12,...
    k1, k2, k3, k4, k5, k6,...
    d1, d2, d3, d4, d5, d6, d7, d8,...
    s1, s2)
    % 21 parameters

    % Boring part. Thanks matlab
    p.td = td; p.t_delay = t_delay; p.t_last = t_last; 
    p.t_delay12 = t_delay12; p.t_last12 = t_last12;
    p.k1 = k1; p.k2 = k2; p.k3 = k3; p.k4 = k4; p.k5 = k5; p.k6 = k6;
    p.d1 = d1; p.d2 = d2; p.d3 = d3; p.d4 = d4; p.d5 = d5; p.d6 = d6;
    p.d7 = d7; p.d8 = d8;
    p.s1 = s1; p.s2 = s2;
    
    % Settings
    treatments_to_fit = {...
        @treatment_factory.cbd_7, @treatment_factory.cbd_9_14, @treatment_factory.cpi_9,...
        @treatment_factory.combo, @treatment_factory.il_12, @treatment_factory.placebo,...
        }; % List of treatments used to perform parameter fitting
    error = 0;
    for treatment = treatments_to_fit
        treatment_spec = treatment{1}();
        sol = immuno_solver(p, treatment_spec);
    end
    error = error + 1;
end

