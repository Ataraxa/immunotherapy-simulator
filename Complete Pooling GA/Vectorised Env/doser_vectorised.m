% Function that tells if a drug is active or not
% Useful to process multi-dose treatments
function is_active = doser_vectorised(t, t_in_array, is_injected, t_delay, t_last)
    is_active = 0;
    for t_in = t_in_array
        if (t_in + t_delay) < t && t < (t_in + t_delay + t_last) && is_injected
            is_active = 1;
            break
        end
    end
end

