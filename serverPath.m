function [dataLoc, exSessDataLoc, saveFigLoc] = serverPath
% Get server path
%
% [dataLoc, exSessDataLoc, saveFigLoc] = serverPath
%
% This function returns the paths for the main data folder, the example 
% session data folder, and the save figure folder for the current system.
% If the hard-coded paths do not exist, the user will be prompted to
% specify them.
% 
% dataLoc = 'F:\Erinn\EL_NatNeuro_2024_data';
% exSessDataLoc = fullfile(dataLoc,'Example Sessions');
% saveFigLoc = 'F:\Erinn\EL_NatNeuro_2024_Figures';
%
% Copyright (c) by Alan Degenhart and Erinn Grigsby
% Emails: alan.degenhart@gmail.com or erinn.grigsby@gmail.com

dataLoc = 'C:\Users\aland\data\EL_NatNeuro_2024_data';
exSessDataLoc = fullfile(dataLoc,'Example Sessions');
saveFigLoc = 'C:\Users\aland\data\EL_NatNeuro_2024_Figures';

% Check if the dataLoc folder exists, if not give a warning and prompt for 
% folder location
check = 1; % Check is the user exited the folder gui
while exist(dataLoc,'dir') ~= 7 && check ~= 0
    warning(['Location provided for dataLoc does not exist.' ...
        ' Please provide a valid folder location'])
    dataLoc = uigetdir([],'Select directory where the paper data is saved');
    check = dataLoc;
end

% Check if the exSessDataLoc folder exists, if not give a warning and 
% prompt for  folder location
check = 1; % Check is the user exited the folder gui
while exist(exSessDataLoc,'dir') ~= 7 && check ~= 0
    warning(['Location provided for exSessDataLoc does not exist.' ...
        ' Please provide a valid folder location'])
    exSessDataLoc = uigetdir([],['Select directory where the example session' ...
        ' data is saved']);
    check = exSessDataLoc;
end

% Check if the saveFigLoc folder exists, if not give a warning and prompt for 
% folder location
check = 1; % Check is the user exited the folder gui
while exist(saveFigLoc,'dir') ~= 7 & ~isempty(saveFigLoc) && check ~= 0
    warning(['Location provided for saveFigLoc does not exist.' ...
        ' Please provide a valid folder location'])
    saveFigLoc = uigetdir([],['Select directory where you would like to ' ...
        'save the paper figure']);
    check = saveFigLoc;
end
