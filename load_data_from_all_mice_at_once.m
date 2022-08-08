%% Description:
% 
% This script loads the previously analyzed and saved data from all the mice in the workspace.
%
% The data needs to be organized in the following format:
% 
% Parent folder  -> one parent folder
% |                     
% ||||animal1           -> one folder per mouse
% |   |   R1 power 5kHz.mat     -> saved file for R1 session
% |   |   R2 power 5kHz.mat     -> saved file for R2 session
% |   |   R3a power 5kHz.mat    -> saved file for R3a session
% |   
% ||||animal2
%     |   R1 power 5kHz.mat      -> saved file for R1 session
%     |   R2 power 5kHz.mat      -> saved file for R2 session
%     |   R3a power 5kHz.mat     -> saved file for R3a session
% 
% Written and tested in MATLAB 2018b.

%%
clear;
close all;
clc;

%% Loading the sub-folders in a folder 

path = uigetdir;            % choose the folder path using GUI
d = dir(path);              % load the folder contents
dfolders = d([d(:).isdir]); % chose only folders
dfolders = dfolders(3:end); % reject first 2 entries as they are empty

%% Load the contents of all the selected folders

prompt = {'Recording Session:', 'Analysis:', 'Sampling frequency:'};
dlgtitle = 'Metadata';
dims = [1 35];
definput = {'R1', 'power', '5kHz'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

s1 = {dfolders.name}; % extrat names of all subfolders

T = {};

for i = 1:length(s1)-1
     destination = strcat(path, '\', s1{i}, '\');
     Files = dir(fullfile(destination, '*.mat'));
     Filenames = {Files.name};
     file = strcat(answer(1), {' '}, answer(2), {' '}, answer(3), '.mat'); % number of channels recorded from
     load(strcat(destination, file{1}));
     T{1, i} = final_baselineZ;
end

%% end of script

