% DDE solver that can process parameterised functions
function sol = immuno_solver(p, tr, plot_info)
    % [input] p -> parameters for the system of DDEs&
    % [input] tr -> treatment to be simulated
    % [output] sol -> results of the simulation (struct)
    
    % This is so OP, thanks matlab
    arguments
        p {struct}
        tr {struct}
        plot_info {struct}
    end

    % Settings of the dde23 solver
    lag = abs(p.td); % TODO: verify that
    tspan = [0 27];
    v_max = 600;

    sol = dde23(@ddefun, lag, @history, tspan);

    function dydt = ddefun(t, y, Z)
        % d_CBD and d_12
        d_cbd = treatment_doser(t, tr.t_in, p.t_delay, p.t_last, tr.active_cbd);
        d_12 = treatment_doser(t, tr.t_in12,p.t_delay12, p.t_last12, tr.active_il12);

        % d_CPT
        if (tr.t_inCPI) < t && tr.active_cpi
            d_cpi = 1;
        else
            d_cpi = 0;
        end
        
        y_lagged = Z(:,1);
        
        dydt = [
          p.k1 + p.k2*(d_cbd + d_12) - p.d1*y(1); % dg/dt
          p.k3 + p.k4*y_lagged(1) - p.d2*y(2); % dc/dt
          p.k5 - (p.d3+p.d4*y(1))*y(3); % dp/dt
          p.k6*(1-(y(4)+y(5))/v_max)*y(4) - (p.d5 + (p.d6*y(2)/(1+p.s1*y(3)*(1-d_cpi)) + p.d7*y(1))/(1+p.s2*(y(4)+y(5))))*y(4); % dv_l/dt
          (p.d5 + (p.d6*y(2)/(1+p.s1*y(3)*(1-d_cpi)) + p.d7*y(1))/(1+p.s2*(y(4)+y(5))))*y(4) - p.d8*y(5);
        ];
    end
    
    function s = history(t)
        s = [
            p.ig0; % Amount of IFNg
            p.c0; % Amount of CD8+
            p.p0; % Amount of PD-1
            p.vl0; % Volume of injected living tumour
            0; % Initial volume of dead tumour
        ];
    end

    if plot_info.flag
        load Data\preprocessed_data.mat data
        figure(1); hold on; % subplot, one for each treatment

        subplot(2, 3, plot_info.subplot_id); hold on
        plot(sol.x,(sol.y(4,:)+sol.y(5,:)),'-','LineWidth',1) % In-silico
        plot(data.(plot_info.treatment_name).valid_days, data.(plot_info.treatment_name).vector_avg, '-o', 'LineWidth',1) % In-vivo
        xlabel('Day since tumour inoculation', 'Interpreter','latex');
        ylabel('Tumour Volume ($mm^3$)', 'Interpreter','latex');
        ylim([0 inf])
        title(regexprep(sprintf('%s', plot_info.treatment_name), '_', ' '), Interpreter="none")
        legend('Simulated', 'Experimental','Location','SouthEast');
    end
end

