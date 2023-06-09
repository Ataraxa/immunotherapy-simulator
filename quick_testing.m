% to test random crap
t = 0;
t_in = [1, 10];
t_delay = 5; t_last = 0.5;
if (t_in + t_delay) < t && t < (t_in + t_delay + t_last) && active_cbd
    d_cbd = 1;
else
    d_cbd = 0;
end