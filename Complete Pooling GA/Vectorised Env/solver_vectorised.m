% DDE solver that can process parameterised functions
function sol = solver_vectorised(p, tr)
    % [input] p -> parameters for the system of DDEs&
    % [input] tr -> treatment to be simulated
    % [output] sol -> results of the simulation (struct)
    
    % This is so OP, thanks matlab
    arguments
        p
        tr {struct}
    end

    % Settings of the dde23 solver
    lag = abs(p(1)); % TODO: verify that
    tspan = [0 27];
    v_max = 600;

    sol = dde23(@ddefun, lag, @history, tspan);

    function dydt = ddefun(t, y, Z)
        % d_CBD and d_12
        d_cbd = treatment_doser(t, tr.t_in, tr.active_cbd, p(2), p(3));
        d_12 = treatment_doser(t, tr.t_in12, tr.active_il12, p(4), p(5));

        % d_CPT
        if (tr.t_inCPI) < t && tr.active_cpi
            d_cpi = 1;
        else
            d_cpi = 0;
        end
        
        y_lagged = Z(:,1);
        
        dydt = [
          p(6) + p(7)*(d_cbd + d_12) - p(12)*y(1); % dg/dt
          p(8) + p(9)*y_lagged(1) - p(13)*y(2); % dc/dt
          p(10) - (p(14)+p(15)*y(1))*y(3); % dp/dt
          p(11)*(1-(y(4)+y(5))/v_max)*y(4) - (p(16) + (p(17)*y(2)/(1+p(20)*y(3)*(1-d_cpi)) + p(18)*y(1))/(1+p(21)*(y(4)+y(5))))*y(4); % dv_l/dt
          (p(16) + (p(17)*y(2)/(1+p(20)*y(3)*(1-d_cpi)) + p(18)*y(1))/(1+p(21)*(y(4)+y(5))))*y(4) - p(19)*y(5);
        ];
    end
end

