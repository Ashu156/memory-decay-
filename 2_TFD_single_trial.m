%% Description:

% This script performs the time-frequency decomposition of a time-domain signal using the Morlet wavelet function. Parameters defined are user-defined through a 
% pop-up dialog box. Finally, it saves the results as a matrix (frequencis X time X trials) averaged over all the trials.

% Written in  MATLAB 2018b.
% Tested in MATLAB 2018b and 2022a.

% Requires the Signal Processing Toolbox.

%% 

tic; % start timing

close all; % close all open tabs in MATLAB
clear;     % clear workspace
clc;       % clear command window

%% Loading the raw data 

[file, path] = uigetfile; % choose the file using GUI
load(strcat(path, file)); % load the file

%% Wavelet transform

clear final_baselineZ % clear this variable as it is not required and occupies a lot of space

prompt = {'Sampling frequency:', 'Minimum frequency:','Maximum frequency:', 'Length of frequency vector:'};
dlgtitle = 'Frequency vector';
dims = [1 35];
definput = {'5000', '0.1',  '500',  '500'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

% Wavelet parameters
Fs = str2double(answer(1));          % To be changed manually as per the downsampling factor 
min_freq = str2double(answer(2));     % minimum of frequency vector
max_freq = str2double(answer(3));     % maximum of  fequency vector
num_frex = str2double(answer(4));     % length of  frequency vector


% Other wavelet parameters
frex = linspace(min_freq, max_freq, num_frex);  % frequency vector
time = -0.02:1/Fs:0.02;                         % time support for Morlet wavelet
half_wave = (length(time) - 1)/2;               % half length of the time support

lfp = lfp_data; 
clear lfp_data 

% FFT parameters
nKern = length(time);                           % length of the kernel

% Initialize output time-frequency data
n_channels = size(lfp, 1);
n_trials = size(lfp, 3);

%% 

tf1 = zeros(length(frex), size(lfp, 2), size(lfp,  3)); % Initiate a zeros matrix to store the results

%%

prompt = {'Which channel do  you want to analyze:'};
dlgtitle = 'Channels';
dims = [1 35];
definput = {'2'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

temp_data1 = squeeze(lfp(str2double(answer(1)),:,:)); % LFP data from channel of choice

nData1 = size(temp_data1, 1)*size(temp_data1, 2);     % number of data points
nConv1 = nKern + nData1 - 1;                          % number of  points in convolution

range_cycles = [ 2  10 ];                             % range of the cycle parameter
cycles = linspace(range_cycles(1), range_cycles(end), num_frex); % cycle vector
num_cycles = length(range_cycles);                               % length of cycle vector

% FFT of total data
fft_lfp1 = fft( reshape(temp_data1, 1, []), nConv1); % FFT of single channel LFP 
    
for cyclei = 1:length(num_cycles) % Loop over all the cycles
        
    for fi = 1:length(frex) % Loop over all the frequencies
            
            % Create wavelet and get its FFT
    s        = cycles(fi) /(2*pi*frex(fi)); % width of the wavelet for every cycle and every frequency

    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2)); % Morlet wavelet function
    
    % Run convolution for each of total, induced, and evoked
   

        % Need separate FFT 
        waveletX = fft(wavelet, nConv1);    % FFT of the Morlet wavelet
        waveletX = waveletX./max(waveletX); % Normalized FFT  
        
        
        %  Notice that the fft_lfp cell changes on each iteration
        as1 = ifft(waveletX.*fft_lfp1, nConv1);   % Multiply FFT of wavelet and FFT of LFP and then perform  inverse FFT on the product
        as1 = as1(half_wave + 1:end - half_wave); % Keeping the rsults for positive frequencies and discarding that for negative frequencies
        as1 = reshape(as1, size(temp_data1, 1), size(temp_data1, 2)); % reshape the result as it was in the input format (i.e, times X trials)
        
            % Compute power
           
            tf1(fi, :, :) = abs(as1).^2;  % store power values for all time points and all trials
            
            
    
             
             
     % end loop around total, evoked, induced
    
    end % end frequency loop

end % end cycle loop


%% saving the results

target = strcat(path, file(1:end-15),' Ch', answer(1), ' trialwise power 5 kHz.mat');
save(target{1}, 'tf1', 'tx', 'frex');
clear tf1


%% end of script
