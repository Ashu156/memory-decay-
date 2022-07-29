%% Description:

% This script plots th
% 1. Trial-averaged time-domain signal (averaged over multiple trials) for all the channels (4 in this case)
% 2. Trial-averaged spectrograms (time-frequency decomposition results averaged over multile trials) for all the channels (4 in this case)

% Basically it plots the results obtained from sript 01.
% INPUT: Results saved from script 01.
% OUTPUT: Plots for results obtained from script 01.

% Written in  MATLAB 2018b.
% Tested in MATLAB 2018b and 2022a.

% Requires the Signal Processing Toolbox.

%% 

tic;       % start timer
close all; % close all open tabs in MATLAB
clear;     % clear workspace
clc;       % clear command window

%% Loading the raw data 

[file, path] = uigetfile; % choose the file using GUI
load(strcat(path, file)); % load the file

%%
figure('Color',[1 1 1])

for i = 1:size(lfp_data, 1)
    subplot(2,2,i)                                               % change this if you have <4 or > 4channels, otherwise you will get an error
    plot(tx, mean(squeeze(lfp_data(i,:, :)), 2), 'linew', 1)
    xlabel('Time (ms)', 'FontSize',  14)
    ylabel('Amplitude (mV)',  'FontSize', 14)
    title(strcat('Channel', {' '}, num2str(i) ))
end

%%
figure('Color', [1 1 1]), clf

for i = 1:size(lfp_data, 1) % Loop over all channels
    subplot(2, 2, i)                                                       % change this if you have <4 or > 4channels, otherwise you will get an error       
    pcolor(tx, frex, squeeze(final_baselineZ(i,  :,  :))); shading interp; % plot the sepctrogram
    colorbar; % show colorbar
    xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14); % X- and Y-axis labels
    hold on % hold on to the figure(subplot) handle
%     plot([baselinetime(2) baselinetime(2)], [min_freq max_freq], '--w', 'LineWidth', 1.5) % plot a vertical dashed white line at the time of event onset
%     xlim([baselinetime(2) - 50 baselinetime(2) + 50]) % Confine X-axis to 50 ms bfore and 50 ms after the event
    ylim([min(frex) max(frex)])
    caxis([-5 10]) % Define the range of the colorbar
    title(strcat('Channel', {' '}, num2str(i) )) % title for the subplot
    hold off % end holding on to the figure(subplot) handle
    
end % end of for loop

%% end of script
