% to test random crap
clc; clear; close all;
treatments = ["treat1", "treat2"];
for tr = treatments
    params.(tr{1}) = 'test';
end