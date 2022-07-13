tic;       % start timer
close all; % close all open tabs in MATLAB
clear;     % clear workspace
clc;       % clear command window

%% Loading the raw data 

[file, path] = uigetfile; % choose the file using GUI
load(strcat(path, file)); % load the file

%% Extract the variables of choice form the workspace 

s1 = whos;               % Extract all variables from the workspace
s2 = {s1.name};          % Extract the name of workspace variables
ind = [];                % Initiate an empty array for storing the indices of variables which are struct array

for k = 1:numel(s1) % Loop over all the workspace variables
    c = s1(k).class; % Extract class of  workspace variables
    
    if strcmp (c, 'struct') == 1 % If the variable is a struct array, then
        ind = [ind, k];          % update the empty array with the index of the variable
    end % end of if-else statment (conditional statement)
    
end      % end of for loop



for i = 1:numel(ind)               % Loop over the elements of the indices array to extract and store only those variables that are struct arrays
    temp = evalin('base',  s2{i}); % Extract data from the struct array
    combined{i} = temp;            % Store the data
end                                % end of for loop


selected = []; % Initiate an empty array for stoing the indices of variables of interest

for k = 1:numel(combined) % Loop over all the struct arrays to choose channels of choice 
    
    % Also check if channel names are E:E1 or just E1
    keyword_target = {'E:E1'; 'E:E2'; 'E:E3'; 'E:E4'; 'LIGHT'}; % Channels of choice (first 4 are for raw data retrievaland the last one is for timestamp  retrieval)
    keyword_retrieved = combined{k}.title;                      % channel name of struct array
    
    for j = 1:numel(keyword_target)   % Loop over the channel names of choice
        
        if strcmp(keyword_retrieved,  keyword_target{j}) == 1 % Compare the channel names
            selected = [selected, k];                         % and store the indices of matched channel names
        end
        
    end
        
end

%% Load the selected data for analysis

timestamps = combined{selected(1)}.times;   % Timestamps for the optical stimulus
Ch1_data = combined{selected(2)}.values;    % Data from the 1st channel
Ch2_data = combined{selected(3)}.values;    % Data from the 2nd channel
Ch3_data = combined{selected(4)}.values;    % Data from te 3rd channel
Ch4_data = combined{selected(5)}.values;    % Data  from the 4th channel

%% parameters for downsampling the raw data

downsampling_factor = 2;          % downsample the raw data by this factor
Fs = 1e4/downsampling_factor;     % new sampling rate after downsampling

%% downsampling the raw data

Ch1_ds = decimate(Ch1_data, downsampling_factor); % downsampled data from channel 1
Ch2_ds = decimate(Ch2_data, downsampling_factor); % downsampled data from channel 2
Ch3_ds = decimate(Ch3_data, downsampling_factor); % downsampled data from channel 3
Ch4_ds = decimate(Ch4_data, downsampling_factor); % downsampled data from channel 4

times = linspace(0, length(Ch1_ds), length(Ch1_ds)-1)/Fs; % vector with timevalues (in samples)
ind = dsearchn(times', timestamps);                       % find indices of events in the time vector

%% Extracting windowed data of choice

n_channels = 4; % number of channels recorded from

t_start = timestamps(1, 1);           % 1st event
ind_start = dsearchn(times',t_start); % index of 1st event in the time vector

time = 1;                 % Extract data around the event for this much time (in seconds)
t_pt = Fs*time + 1;       % number of time points
step_size = Fs/2;         % step size

lfp_data = zeros(n_channels, t_pt, length(timestamps)); % initiate an empty vector to store data for all channels 
data = zeros(t_pt, n_channels);                         % initiate temporary vector to store data


for i = 1:length(timestamps)         % Loop over number of events
    t_start = timestamps(i, 1);      % event timestamp
    idx = dsearchn(times', t_start); % index  of event in the time vector
    data(:, 1) = Ch1_ds(idx - step_size:idx + step_size, 1); % 
    data(:, 2) = Ch2_ds(idx - step_size:idx + step_size, 1);
    data(:, 3) = Ch3_ds(idx - step_size:idx + step_size, 1);
    data(:, 4) = Ch4_ds(idx - step_size:idx + step_size, 1);
    lfp_data(1, :, i) = data(:, 1); % data from 1st channel
    lfp_data(2, :, i) = data(:, 2); % data from  2nd channel
    lfp_data(3, :, i) = data(:, 3); % data from 3rd channel
    lfp_data(4, :, i) = data(:, 4); % data from 4th channel
    
end

% clear unnecessary variables from workspace to free up memory
clear temp c combined Ch1_data Ch2_data Ch3_data Ch4_data Ch1_ds Ch2_ds Ch3_ds Ch4_ds data keyword_target keyword_retrieved selected

%% Plot the evoked response

figure('Color',[1 1 1])

for i = 1:4
    subplot(2,2,i)
    plot(mean(squeeze(lfp_data(i,:, :)), 2), 'linew', 1)
    xlabel('Time (samples)', 'FontSize',  14)
    ylabel('Amplitude (mV)',  'FontSize', 14)
    title(strcat('Channel', {' '}, num2str(i) ))
end

%% Wavelet transform (Time-frequency decomposition)

% Wavelet parameters
min_freq = 0.1; % minimum frequency
max_freq = 500; % maximum frequency
num_frex = 500; % number of  frequencies

% Other wavelet parameters
frex = linspace(min_freq, max_freq, num_frex); % frequency vector
time = -0.02:1/Fs:0.02;                        % time vector for wavelet transform
half_wave = (length(time) - 1)/2;              % half of the length of the time vector

range_cycles = [ 2  10 ];                                        % range of cycle paramter
cycles = linspace(range_cycles(1), range_cycles(end), num_frex); % cycle vector
num_cycles = length(range_cycles);                               % length of cycl vector

TF = zeros(n_channels, length(frex), size(lfp_data, 2)); % initiating a zero matrix for storing the TF decomposition data

for k = 1:n_channels                 % Loop over all the channels
    lfp = squeeze(lfp_data(k,:,:));  % data from each channel

    % FFT parameters
    nKern = length(time);                 % length of the kernel
    nData = size(lfp, 1)*size(lfp, 2);    % number of data points
    nConv = nKern + nData - 1;            % number of points in convolution

   % Initialize output time-frequency data

    tf = zeros(length(frex), size(lfp, 1)); % Intiate a temporary zero  matrix for storing TF data
   


     % FFT of total data
     fft_lfp = fft( reshape(lfp, 1, []), nConv); % FFT of single channel LFP 




       for cyclei = 1:length(num_cycles) % Loop over all the cycles
    
            for fi = 1:length(frex)      % Loop over all the frequencies
    
               % Create wavelet and get its FFT
             s        = cycles(fi) /(2*pi*frex(fi)); % width of the wavelet for every cycle and every frequency

             wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2)); % Morlet wavelet function
    
             % Run convolution for each of total, induced, and evoked
   

               % Need separate FFT 
               waveletX = fft(wavelet, nConv);       % FFT of the Morlet wavelet
               waveletX = waveletX./max(waveletX);   % Normalized FFT  
        
        
        %  Notice that the fft_lfp cell changes on each iteration
        as = ifft(waveletX.*fft_lfp, nConv);   % Multiply FFT of wavelet and FFT of LFP and then perform  inverse FFT on the product
        as = as(half_wave + 1:end - half_wave); % Keeping the rsults for positive frequencies and discarding that for negative frequencies
        as = reshape(as, size(lfp, 1), size(lfp, 2)); % reshape the result as it was in the input format (i.e, times X trials)
        
            % Compute power
            tf(fi, :) = mean(abs(as).^2, 2); % TF decomposition result averaged over all the trials
             
      
     % end loop around total, evoked, induced
end % end frequency loop

end % end cycle loop

TF(k,  :, :) = reshape(tf, 1, size(tf, 1), size(tf, 2)); % store the results for each channel

end % end of channel loop

tx = linspace(0, size(lfp, 1), size(lfp, 1))./ Fs*1e3; % time vector in ms

%% Baseline corrected power

baselinetime = [400 500]; % baseline time for normaliation (in ms)

% Convert baseline window time to indices
[~, baselineidx(1)] = min(abs(tx - baselinetime(1))); % Find starting baseline time in the time vector
[~, baselineidx(2)] = min(abs(tx - baselinetime(2))); % Find ending baseline time in the time vector

final_dbconverted = zeros(n_channels, numel(frex), numel(tx)); % Initiate a zero matrix for storing baseline-normaized dB converted TF data
final_pctchange = zeros(n_channels, numel(frex), numel(tx));   % Initiate a zero matrix for storing baseline-normaized percent changed TF data
final_baselinediv = zeros(n_channels, numel(frex), numel(tx)); % Initiate a zero matrix for storing baseline-normaized ratio change TF data
final_baselineZ = zeros(n_channels, numel(frex), numel(tx));   % Initiate a zero matrix for storing baseline-normaized z-scored TF data


for i = 1:n_channels % Loop over all the channels
    
    % dB converted power
    baseline_power = mean(squeeze(TF(i, :, baselineidx(1):baselineidx(2))), 2); % baseline power
    dbconverted = 10*log10( bsxfun(@rdivide, squeeze(TF(i, :, :)), baseline_power)); % baseline-normalized dB converted data
    final_dbconverted(i, :,  :) = dbconverted;    % store results for each  channel

    % Percent change in power w.r.t. baseline
    pctchange = 100 * (squeeze(TF(i, :, :)) - repmat(baseline_power, 1, size(lfp, 1)))./ repmat(baseline_power, 1, size(lfp, 1));
    final_pctchange(i, :,  :) = pctchange; % store results for each channel

    % Baseline division
    baselinediv = squeeze(TF(i, :, :)) ./ repmat(baseline_power, 1, size(lfp, 1));
    final_baselinediv(i, :,  :) = baselinediv; % store results for each channel

    % Z-transform
    baseline_powerZ = squeeze(TF(i,  :, baselineidx(1):baselineidx(2))); % baseline power
    baselineZ = (squeeze(TF(i, :, :)) - repmat(mean(baseline_powerZ, 2), 1, size(TF, 3)))...
                ./ repmat(std(baseline_powerZ, [], 2), 1, size(TF, 3));
    final_baselineZ(i, :, :) = baselineZ; % stor results for each channel

end

%% Plotting the results

figure('Color', [1 1 1]), clf

for i = 1:n_channels % Loop over all channels
    subplot(2, 2, i) % Subplots in the figure
    pcolor(tx, frex, squeeze(final_baselineZ(i,  :,  :))); shading interp; % plot the sepctrogram
    colorbar; % show colorbar
    xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14); % X- and Y-axis labels
    hold on % hold on to the figure(subplot) handle
    plot([500 500], [min_freq max_freq], '--w', 'LineWidth', 1.5) % plot a vertical dashed white line at the time of event onset
    xlim([450 550]) % Confine X-axis to 450 to 550 ms
    ylim([0.1 500])
    caxis([-5 10]) % Define the range of the colorbar
    title(strcat('Channel', {' '}, num2str(i) )) % title for the subplot
    hold off % end holding on to the figure(subplot) handle
    
end % end of for loop

%% saving the results

save(strcat(path, 'R3 before pretest power 5kHz.mat'), 'lfp_data', 'final_baselineZ', 'tx',  'frex'); % save the results
toc; % stop the timer

%% end of script
