function KD = getKinematicData(TD,varargin)
% [KD] = getKinematicData(TD,kinSource)
%
% Get kinematic data method for TrajectoryData class
%
% Author:  Alan D. Degenhart
% Date Created: 2016/06/21
% Last Updated: 2016/06/21
% Last Update:  Initial version of function
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

kinSource = [];         % Kinematic source to plot ('hand' or 'brain')

% Parse optional agruments
assignopts({'kinSource'},varargin);

% Determine if kinematic source is specified or should be determined
% automatically
useDefaultKinSource = true;
if ~isempty(kinSource)
    useDefaultKinSource = false;
end

% Figure out the class of the desired kinematic data source.  Assume this
% is the same for all trials in the TrajectoryData object
if useDefaultKinSource
    switch TD(1).controlSource
        case {'Neural Decoder','Auto-monkey'}
            kinSource = 'brain';
        case {'Phasespace'}
            kinSource = 'hand';
        case {'Force Cursor'}
            kinSource = 'force';
    end
end
srcStr = [kinSource 'Kin'];
KDtemp = TD(1).(srcStr);

% Loop over all TrajectoryData elements and get KinematicData
nTrials = length(TD);
KD = repmat(KDtemp,nTrials,1);
for i = 1:nTrials
    % Get kinematic data.  This can be either the hand position or the
    % brain-controlled cursor position.  By default, the cursor being
    % used for task progression is plotted.
    KD(i) = TD(i).(srcStr);
end