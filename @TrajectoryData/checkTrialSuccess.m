% checkTrialSuccess         Verify trial success status for TrajectoryData
% object.
%
% TD = TD.checkTrialSuccess
%
% This function verifies that the 'successful' status of the trial(s) in
% the provided TrajectoryData object are consistent with the kinematic
% data.  This is used to handle cases where mismatches occur due to
% data-saving issues, particulary RT errors.
%
% Created:  2018.10.29
% Author:   Alan D. Degenhart
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function TD = checkTrialSuccess(TD,varargin)

% Optional arguments
checkStartPosition = false;
startPosRadius = 10;
checkTargetPosition = true;
dispWarning = true;
debugMode = false;
radiusScale = 1;

assignopts(who,varargin);

nTrials = length(TD);

% Loop over trials
scValid = true(nTrials,1);
dTarg = nan(nTrials,1);
dStart = dTarg;
for i = 1:nTrials
    % Get target position and radius
    tPos = TD(i).targPos;
    sPos = TD(i).startPos;
    tRad = TD(i).targSize;
    K = TD(i).getKinematicData;
    pos = K.pos';
    
    % Check to see if the last sample of the position data lies within the
    % target position
    if checkTargetPosition
        dTemp = pos(:,end) - tPos;
        dTarg(i) = norm(dTemp(1:2));
        if dTarg(i) < tRad*radiusScale; scTarg = true; else; scTarg = false; end
    end
    
    % Check to see if the first sample of the position data lies close to
    % the start position
    scStart = true;
    if checkStartPosition
        dTemp = pos(:,1) - sPos;
        dStart(i) = norm(dTemp(1:2));
        if dStart(i) < startPosRadius; scStart = true; else; scStart = false; end
    end
    
    % If current trial was not successful but was marked as such, update
    % the 'scValid' mask
    if (~scTarg || ~scStart) && logical(TD(i).successful)
        scValid(i) = false;
    end
end

nInvalidTrials = sum(~scValid);
invalidTrialList = [TD(~scValid).trialID];
TD = TD(scValid);

if debugMode
    TDinvalid = TD(~scValid);
    dInvalid = d(~scValid);
    for i = 1:length(TDinvalid)
        F = figure;
        TDinvalid(i).plot
        title(sprintf('d = %0.3f',dInvalid(i)))
        axis equal
        pause
        close(F)
    end
end

if dispWarning && nInvalidTrials > 0
    % Loop over invalid trials and generate string
    s = sprintf(' %d',invalidTrialList);
    warning('%d invalid trial(s) detected:%s. These trials have been removed.\n', ...
        nInvalidTrials,s)
end