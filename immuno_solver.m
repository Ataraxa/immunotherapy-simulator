% DDE solver that can process parameterised functions
function sol = immuno_solver(p, tr)
    % [input] p -> parameters for the system of DDEs (struct)
    % [input] tr -> treatment to be simulated (struct)
    % [output] sol -> results of the simulation (struct)
    
    % Settings of the dde23 solver
    lag = p.td;
    tspan = [0 27];
    v_max = 600;

    sol = dde23(@ddefun, lag, @history, tspan);

    function dydt = ddefun(t, y, Z)
        % d_CBD and d_12
        d_cbd = treatment_doser(t, tr.t_in, p, tr.active_cbd);
        d_12 = treatment_doser(t, tr.t_in12, p, tr.active_il12);

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
end

