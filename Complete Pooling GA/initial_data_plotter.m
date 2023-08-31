data = readtable("tumour_vol_data.csv");

%% Create a line for each mouse
clc; close 
figure(1); hold on
for mouse_id = 1:9
    fprintf('Mouse %i \n', mouse_id)
    individual_array = zeros(9, 1);
    for day_id = 1:9
        
        row = 100 + mouse_id + 9 * (day_id-1) - 1;
        fprintf('%i: %f \n',row, data{row,3})
        individual_array(day_id) = data{row, 3};
    end
    plot(individual_array, '-o')
    disp('____________________________________________')
end