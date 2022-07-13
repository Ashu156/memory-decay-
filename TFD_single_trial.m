tic; % start timing

close all; % close all open tabs in MATLAB
clear;     % clear workspace
clc;       % clear command window

%% Loading the raw data 

[file, path] = uigetfile;
load(strcat(path, file));

%% Wavelet transform

clear final_baselineZ

% Wavelet parameters
min_freq = 0.1;     % minimum of frequency vector
max_freq = 500;     % maximum of  fequency vector
num_frex = 500;     % length of  frequency vector
Fs = 5000;          % To be changed manually as per the downsampling factor 

% Other wavelet parameters
frex = linspace(min_freq, max_freq, num_frex);
time = -0.02:1/Fs:0.02;
half_wave = (length(time) - 1)/2;

lfp = lfp_data;


% FFT parameters
nKern = length(time);

% Initialize output time-frequency data
n_channels = size(lfp, 1);
n_trials = size(lfp, 3);

clear lfp_data 

%% 

tf1 = zeros(length(frex), size(lfp, 2), size(lfp,  3));

%%

temp_data1 = squeeze(lfp(3,:,:)); % To be changed manually for chnnelof choice,e.g., here it is for Ch2

nData1 = size(temp_data1, 1)*size(temp_data1, 2); % number of data points
nConv1 = nKern + nData1 - 1;                      % number of  points in convolution

range_cycles = [ 2  10 ];                          % range of the cycle parameter
cycles = linspace(range_cycles(1), range_cycles(end), num_frex);
num_cycles = length(range_cycles);

% FFT of total data
fft_lfp1 = fft( reshape(temp_data1, 1, []), nConv1);
    
for cyclei = 1:length(num_cycles)
        
    for fi = 1:length(frex)
            
            % Create wavelet and get its FFT
    s        = cycles(fi) /(2*pi*frex(fi));

    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    
    % Run convolution for each of total, induced, and evoked
   

        % Need separate FFT 
        waveletX = fft(wavelet, nConv1);
        waveletX = waveletX./max(waveletX);
        
        
        %  Notice that the fft_lfp cell changes on each iteration
        as1 = ifft(waveletX.*fft_lfp1, nConv1);
        as1 = as1(half_wave + 1:end - half_wave);
        as1 = reshape(as1, size(temp_data1, 1), size(temp_data1, 2));
        
            % Compute power
           
            tf1(fi, :, :) = abs(as1).^2;  % store powr values for all time points and all trials
            
            
    
             
             
     % end loop around total, evoked, induced
    
    end % end frequency loop

end % end cycle loop

    
tx = linspace(0, 1000, size(lfp, 2)); % time vector in ms

%%

save(strcat(path, 'R3 before pretest_Ch3 trialwise power 5 kHz.mat'), 'tf1', 'tx', 'frex');
clear tf1

%% end of script
