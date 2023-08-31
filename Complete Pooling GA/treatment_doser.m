% Function that tells if a drug is active or not
% Useful to process multi-dose treatments
function is_active = treatment_doser(t, t_in_array, delay, last, is_injected)
    is_active = false;
    for t_in = t_in_array
        if (t_in + delay) < t && t < (t_in + delay + last) && is_injected
            is_active = true;
            break
        end
    end
end

