function [TD,directoryName] = load(TD,directoryName)
% Load method for TrajectoryData class
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

if nargin == 1
    [directoryName] = uigetdir;
end

% Get all files in directory
D = dir(directoryName);
fNames = {D.name};
nFiles = length(D);

% Load trial info data
fileInd = cellfun(@strfind,fNames,repmat({'_trialInfo.mat'},1,nFiles), ...
    'UniformOutput',false);
fileInd = ~cellfun(@isempty,fileInd);
trialInfo = load([directoryName '/' fNames{fileInd}]);
trialInfo = trialInfo.trialInfo;

% Initialize trial structure
TD = repmat(TrajectoryData(),trialInfo.nTrials,1);

% Display message to user
fprintf('Loading %d trial(s) into TrajectoryData object ... ', ...
    trialInfo.nTrials)

% Load trial data.
for i = 1:trialInfo.nTrials
    trialName = sprintf('Trl_%0.4d.mat',trialInfo.trialNo(i));
    tempTrial = load([directoryName '/' trialName]);
    TD(i) = tempTrial.tempTrial;
end

fprintf('done.\n')