% File to visualise the CDB-D7 experimental results

clc; clear; close all;

raw_data = readtable("../Data/CBD-D7_data.csv");
vectorized_data = horzcat(table2array(raw_data(:, 2)), table2array(raw_data(:, 4:end)));


%% Graphs Settings
total_mice = 5;
CR_pool = [3]; % Select the mice that achieved CR
collection_days = [7 8 9 11 14 17 20];
data_per_collection_day = [5 4 5 5 5 4 5];
% ... more settings ?

%% Plotting the data
collector = [0, cumsum(data_per_collection_day(1:end-1))];

% Tumour data
tumour_data = zeros(total_mice, numel(collection_days), 3);

for id = 1:total_mice
    tumour_data(id, :, 1) = vectorized_data(collector+id, 5); % CD8+
    tumour_data(id, :, 2) = vectorized_data(collector+id, 13); % PD1
    tumour_data(id, :, 3) = vectorized_data(collector+id, end); % tumour
end

% Create and setup figure
figure(1); hold on;
% subplot(2, 1, 1); hold on
% title("Evolution of tumour volume in mice who achieved CR")
% subplot(2, 1, 2); hold on
% title("Evolution of tumour volume in mice who did not achieved CR")

for id = 1:total_mice
    if any(CR_pool(:) == id)
        subplot(2, 3, 1); hold on
        plot(collection_days, tumour_data(id,:,1), '-o')
        subplot(2, 3, 2); hold on
        plot(collection_days, tumour_data(id,:,2), '-o')
        subplot(2, 3, 3); hold on
        plot(collection_days, tumour_data(id,:,3), '-o')
    else
        subplot(2, 3, 4); hold on
        plot(collection_days, tumour_data(id,:,1), '-o')
        subplot(2, 3, 5); hold on
        plot(collection_days, tumour_data(id,:,2), '-o')
        subplot(2, 3, 6); hold on
        plot(collection_days, tumour_data(id,:,3), '-o')
    end
end
