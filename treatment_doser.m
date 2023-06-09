% Function that tells if a drug is active or not
% Useful to process multi-dose treatments
function is_active = treatment_doser(t, t_in_array, p, is_injected)
    is_active = 0;
    for t_in = t_in_array
        if (t_in + p.t_delay) < t && t < (t_in + p.t_delay + p.t_last) && is_injected
            is_active = 1;
            break
        end
    end
end

