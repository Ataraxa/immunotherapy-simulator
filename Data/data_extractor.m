% Script to interface with the experimental data in an automated way
% It simply reads the files and saved them in a more convenient way (mat
% file)

% How does the final strcture work?
% data -> treatment_name -> 
% ind_info = a matrix with the individual tumour volume data. In row
%            direction are indexed patiens, in col direction are days. It is ordered
%            in increasing days, but does not contain the exact timestamp
%
% valid_days = an array with the days where data has been collected
clc; clear; close all;
tumour_data = readtable("Data/tumour_vol_data.csv");
ifng_data = readtable("Data/ifng_data.csv");
cbd_d7_data = readtable("Data/CBD-D7_data.csv");

%% Dealing with the tumour data
% Organise tumour as struct
% Each element of 1st-order struct is treatment
% In each treatment there is array of days, mean and std (all same size),
% and a matrix with all the individual info
clc;

% Settings and Initialisation
max_days = 30; max_nb_patients = 30;
data = struct();

% Process each row individually
for data_point = 1:size(tumour_data, 1)
    % Individual data processing and refactoring
    treatment_name = toValidFieldName(tumour_data{data_point, "treatment"}{1});
    day = tumour_data{data_point, "daysSinceTumourInoculation"};

    % Initialise matrix if first encounter
    if ~(isfield(data, treatment_name))
        data.(treatment_name).ind_info = zeros(max_days, max_nb_patients); % rows are patients and cols are days
        data.(treatment_name).valid_days = day;
    end
    
    % Update the ind_info (individual info) matrix (change first zero-ed cell)
    for patient = 1:max_nb_patients
        if data.(treatment_name).ind_info(patient, day) == 0
            data.(treatment_name).ind_info(patient, day) = tumour_data{data_point, "tumourVolume_mm_3_"};
            break
        end
    end
   
    % Update the valid_days array (only add day if not already in array)
    if (data.(treatment_name).valid_days(:) ~= day)
        data.(treatment_name).valid_days(end+1) = day;
    end
end

% Reprocess each treatment field 
treatment_list = fieldnames(data);
for treatment_name = treatment_list.' % Need to take transpose of array to loop over it
    slice = data.(treatment_name{1}).ind_info;
    data.(treatment_name{1}).ind_info( ~any(slice,2), : ) = []; % removes empty rows
    data.(treatment_name{1}).ind_info( :, ~any(slice,1) ) = []; % removes empty columns
    
    re_slice =  data.(treatment_name{1}).ind_info;
    data.(treatment_name{1}).vector_avg = mean(re_slice, 1);
    data.(treatment_name{1}).vector_std = std(re_slice, 1);
    data.(treatment_name{1}).vector_min = min(re_slice);
    data.(treatment_name{1}).vector_max = max(re_slice);
end

save("preprocessed_data.mat", "data")
fprintf("Processing done!\n")

% Auxiliary functions
function valid_str = toValidFieldName(invalid_str)
    % Function that converts an invalid field name to a valid one by
    % removing spaces and dashes
    valid_str = regexprep(invalid_str, ' ', '_'); % replaces space by underscore
    valid_str = regexprep(valid_str, '-', '_'); % replaces - by _
    valid_str = regexprep(valid_str, ',', ''); % removes commas
    valid_str = regexprep(valid_str, '+', 'and'); % replaces + by 'and'
end
