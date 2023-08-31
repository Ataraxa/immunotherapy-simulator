% Code to plot the experimental data using boxplot
clc; clear; close all
fprintf("starting running the program...\n")
figure(1); hold on; % Figure will be 7x4

% Reading the files
tumour_data = readtable("Data/tumour_vol_data.csv");
ifng_data = readtable("Data/ifng_data.csv");
cbd_d7_data = readtable("Data/CBD-D7_data.csv");

% TODO: reverse row order of table (for consistency)
% tumour_data = tumour_data{:,:};
treatments = ["Placebo", "CBD-IL-12 Day 7", "CBD-IL-12 Day 9, 14", "IL-12 Day 7'", "CPI Day 9, 14", "CBD-IL-12 + CPI Day 9, 14"];
matrix_treat = zeros(30, 30, numel(treatments));
%%
subject_id = 1;
for row = 1:334
    temp = tumour_data{row,2:end};
    treatment_id = find(treatments==tumour_data{row,1}); % Index of the correct plane
    matrix_treat(temp(1), subject_id, treatment_id) = temp(2);
    if row == 1
        continue
    elseif tumour_data{row-1,2} == tumour_data{row, 2}
        subject_id = subject_id + 1;
    else
        subject_id = 1;
    end
end

for treatment = 1:numel(treatments)
    slice = matrix_treat(:,:,treatment);
    row_id = treatment + 3*(treatment - 1);
    slice( ~any(slice,2), : ) = []; % rows
    slice( :, ~any(slice,1) ) = []; % columns
    subplot(6,4,row_id); hold on
    title(treatments(treatment))
    vector_min = min(slice, 1);
    vector_max = max(slice, 1);
    vector_avg = mean(slice, 1);
    errorbar((vector_max+vector_min)/2, (vector_max-vector_min)/2);
end
fprintf("exit code 1 \n")


