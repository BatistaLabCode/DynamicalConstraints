function KD = getKinematicData(TD,varargin)
% [KD] = getKinematicData(TD,kinSource)
%
% Get kinematic data method for TrajectoryData class
%
% Usage:
%   KD = TD.getKinematicData()
%
% Inputs:
%   TD             TrajectoryData class
%
% Optional Inputs:
%   kinSource      Kinematic source to plot ('brain')
%
% Ouputs:
%   KD             KinematicData class
%
% Author:  Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

kinSource = 'brain';         % Kinematic source to plot ('brain')

% Parse optional agruments
assignopts({'kinSource'},varargin);

% Create a template for the KinematicData structure
srcStr = [kinSource 'Kin'];
KDtemp = TD(1).(srcStr);

% Loop over all TrajectoryData elements and get KinematicData
nTrials = length(TD);
KD = repmat(KDtemp,nTrials,1);
for i = 1:nTrials
    % Get kinematic data. This will be the brain-controlled cursor position.
    KD(i) = TD(i).(srcStr);
end