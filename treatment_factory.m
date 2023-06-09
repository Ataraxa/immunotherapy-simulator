classdef treatment_factory
    % Class of object that can produce treatments

    methods(Static)
        function treat_struct = placebo()
            treat_struct.t_in = 1;
            treat_struct.t_in12 = 1;
            treat_struct.t_inCPI = 1;
            treat_struct.active_cbd = false;
            treat_struct.active_il12 = false;
            treat_struct.active_cpi = false;
        end
        
        function treat_struct = cbd_7()
            treat_struct.t_in = 7;
            treat_struct.t_in12 = 1;
            treat_struct.t_inCPI = 1;
            treat_struct.active_cbd = true;
            treat_struct.active_il12 = false;
            treat_struct.active_cpi = false;
        end

        function treat_struct = cbd_9_14() % NOT FINISHED YET
            treat_struct.t_in = 9;
            treat_struct.t_in12 = 1;
            treat_struct.t_inCPI = 1;
            treat_struct.active_cbd = true;
            treat_struct.active_il12 = false;
            treat_struct.active_cpi = false;
        end

        function treat_struct = il_12()
            treat_struct.t_in = 1;
            treat_struct.t_in12 = 7;
            treat_struct.t_inCPI = 1;
            treat_struct.active_cbd = false;
            treat_struct.active_il12 = true;
            treat_struct.active_cpi = false;
        end
    end
end



