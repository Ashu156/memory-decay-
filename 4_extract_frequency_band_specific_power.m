%%

% This script extracts the user-defined frequency band specific power from the single trial time-frequency decomposition

% Basically it plots the results obtained from sript 01.
% INPUT: Results saved from script 03.
% OUTPUT: User-defined frequency band specific power for all trials and then segregates this into blocks where blocks are represent power averaged over
%         certain number of trials defined the user.

% Written in  MATLAB 2018b.
% Tested in MATLAB 2018b and 2022a.

%%

tic;       % start timer
close all; % close all open tabs in MATLAB
clear;     % clear workspace
clc;       % clear command window

%% Loading the raw data 

[file, path] = uigetfile; % choose the file using GUI
load(strcat(path, file)); % load the file

%% Define values of some variables 

prompt = {'Sampling frequency (in Hz):','Minimum frequency (in Hz):','Maximum frequency (in Hz):', 'Length of frequency vector:'};
dlgtitle = 'Frequency input';
dims = [1 50];
definput = {'5000', '1', '500', '500'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%% Baseline-corrected power

Fs = str2double(answer(1)); % Sampling frequency
min_freq = min(frex);       % minimum of frequency vector
max_freq = max(frex);       % maximum of frequency vector
num_frex = length(frex);    % length of frequency vector
% frex = linspace(min_freq, max_freq, num_frex);

prompt = {'Baseline starts at (in ms):','Baseline ends at (in ms):'};
dlgtitle = 'Baseline';
dims = [1 35];
definput = {'400', '500'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

baselinetime = [str2double(answer(1)) str2double(answer(2))]; % baseline time for normaliation (in ms)
tx = linspace(0, 1*1000, size(tf1, 2)); % time vector in ms

% Convert baseline window time to indices
[~, baselineidx(1)] = min(abs(tx - baselinetime(1)));
[~, baselineidx(2)] = min(abs(tx - baselinetime(2)));

% final_baselineZ = zeros(size(tf1, 1), size(tf1, 2), size(tf1, 3));


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

prompt = {'Starting frequency (in Hz):','Ending frequency (in Hz):'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'30','120'};
freq_band = inputdlg(prompt,dlgtitle,dims,definput);


freq = [str2double(freq_band(1)) str2double(freq_band(2))];
freq_ind = dsearchn(frex',freq');

prompt = {'Starting time (in ms):','Ending time (in ms):'};
dlgtitle = 'Input';
dims = [1 35];
definput = {' ',' '};
time_band = inputdlg(prompt,dlgtitle,dims,definput);

time = [str2double(time_band(1)) str2double(time_band(2))];
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

prompt = {'Number of  trials to  be averaged in a block:'};
dlgtitle = 'Block size';
dims = [1 35];
definput = {'5'};
block_size = inputdlg(prompt,dlgtitle,dims,definput);

n = str2double(block_size(1)); % block size 
block_average = arrayfun(@(i) mean(power(i:i+n-1)),1:n:length(power)-n+1)'; % the averaged vector


%% 
save(strcat(path, file(1:end-25), ' gamma power 10 ms.mat'), 'pwr_temp', 'block_average')

figure;
plot(block_average, 'o-k', 'linew',  2)

%% end of script
