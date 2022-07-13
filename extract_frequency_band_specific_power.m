close all;
clear;
clc;

%% Loading the raw data

[file, path] = uigetfile;
load(strcat(path, file));

%% Baseline-corrected power

Fs = 5e3;
min_freq = min(frex);
max_freq = max(frex);
num_frex = length(frex);
% frex = linspace(min_freq, max_freq, num_frex);

baselinetime = [ 400 500 ]; % in ms
tx = linspace(0, 1000, size(tf1, 2)); % time vector in ms

% Convert baseline window time to indices
[~, baselineidx(1)] = min(abs(tx - baselinetime(1)));
[~, baselineidx(2)] = min(abs(tx - baselinetime(2)));

final_baselineZ = zeros(size(tf1, 1), size(tf1, 2), size(tf1, 3));


for  i = 1:size(tf1, 3)
       
        temp_tf = squeeze(tf1(:, :,  i));
        baseline_power = mean(temp_tf(:, baselineidx(1):baselineidx(2)), 2);
    
        % Z-transform
        baseline_powerZ = temp_tf(:, baselineidx(1):baselineidx(2));
        baselineZ = (temp_tf - repmat(mean(baseline_powerZ, 2), 1, size(temp_tf, 2)))...
                     ./ repmat(std(baseline_powerZ, [], 2), 1, size(temp_tf, 2));

        final_baselineZ(:,  :, i) = baselineZ;
end

%% Extract frequency band specific power

freq = [30 120];
freq_ind = dsearchn(frex',freq');

time = [500 510];
time_ind = dsearchn(tx',time');
% baselineZ = squeeze(final_baselineZ(1, :, :));
% n_channels = size(final_baselineZ, 1);
n_trials = size(final_baselineZ, 3);
power = zeros(n_trials, 1);
pwr_temp = zeros(freq_ind(2)-freq_ind(1)+1, time_ind(2)-time_ind(1)+1,  n_trials);
% for j = 1:n_channels

    
    for i = 1:length(power)
        temp_data = squeeze(final_baselineZ(:, :, i));
        pwr = temp_data(freq_ind(1):freq_ind(2), time_ind(1):time_ind(2));
        pwr_temp(:, :,  i) = pwr;
        pwr_sum = sum(sum(pwr, 2), 1);
        power(i, 1) = pwr_sum;
    end
    

% end
n = 10; % block size 
block_average = arrayfun(@(i) mean(power(i:i+n-1)),1:n:length(power)-n+1)'; % the averaged vector


%% 
save(strcat(path, file, 'gamma power 10 ms.mat'), 'pwr_temp', 'block_average')

figure;
plot(block_average, 'o-k', 'linew',  2)

%% end of script
