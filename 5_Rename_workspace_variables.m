%%
% This script renames the variables loaded into the workspace so that the animal ID is concatenated to it. For example, if the loaded variables are 
% tf, frex, tx, etc. and come from the animal G1 C1 0L0R, then this script will rename those as tf_G1_C1_0L0R, frex_G1_C1_0L0R, tx_G1_C1_0L0R, etc. respectively 
% and saves them in the respective folder. This will help in loading the data from all the animals at once.

% Written and tested in MATLAB 2018b.

%%
clear;
close all;
clc;

%% Loading the raw data 

[file, path] = uigetfile; % choose the file using GUI
load(strcat(path, file)); % load the file

%% 
s1 = whos;               % Extract all variables from the workspace
s2 = {s1.name};          % Extract the name of workspace variables

%% Renaming the old varibales such that they have the rat ID

ID = path(end-10:end-1);
ID_new = ID;
ID_new(ID == ' ') = '_'; 
animal_ID = strcat('_', ID_new); % rat ID

% Old varibales to be renamed

oldname = s2;
oldname(strncmpi(oldname,'path', length(oldname))) = [];
oldname(strncmpi(oldname,'file', length(oldname))) = [];
% oldname = {'lfp'; 'tx'};

% Initializing the cell array for storing renamed variables
newname = cell(size(oldname));

% Renaming the old variables
for i = 1:length(oldname)
    
    % Check if the old variables exist in the workspace
    teststr = ['exist(''', oldname{1, i}, ''', ''var'')'];
    result = evalin('base', teststr); 
    
    % If the old variables don't exist, throw an error
    if result ~= 1
        error(['A variable named ''', oldname{1, i}, ''' does not exist in the base workspace'])
    end
    
    % New variable names
    newname{1, i} = strcat(oldname{1, i}, animal_ID);
    
    % Assigning the old values to the renamed variables
    str = [newname{1, i}, ' = ', oldname{1, i}, ';'];
    
    try
        evalin('base',str)
    catch
        error('The rename failed')
    end
end

%% Saving the renamed variables

save(strcat(path, file(1:end-4)), newname{:}); % save the results

%% end of script
