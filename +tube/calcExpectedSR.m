% Calculate expected success rate for constrained path experiment
%
% Usage:
%   [sr, sr_all] = tube.calcExpectedSR(TD,TU,r)
% 
% Inputs:
%   TD      TrajectoryData object containing trials to apply tube
%               constraint to
%   TU      Base tube object with path
%   r       Array of tube radii to assess
%
% Optional Inputs:
%   startPos    Start target position for the constrained path task
%   kinSource   Determine which kinSource to calculate the success rate
%   rmFail      Remove failed trials
%   plotSegment Determine we are plotting a segment of the tube around
%                      the path or the entire tube around the path
%
% Output:
%   sr          Success Rate for each tube radius
%   sr_all      Success state for all trials for each radius
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [sr, sr_all] = calcExpectedSR(TD,TU,r,varargin)

% Parse optional arguments
startPos = [];
kinSource = 'brain';
rmFail = 0;
plotSegment = [];

assignopts(who,varargin);

if ~isempty(startPos)
    startPos_all = [TD.startPos]';
    mask = ismember(startPos_all, startPos, 'rows');
    TD = TD(mask);
end

% Remove any failed trials from the trajectory data object
if rmFail
    sc = logical([TD.successful]);
    TD = TD(sc);
end

n_trials = length(TD);
n_tube = length(r);
sr = nan(1,n_tube);

srcStr = [kinSource 'Kin'];
% Loop over tube radii
sr_all = cell(n_tube, 1);
for i = 1:n_tube
    % Create tube object with same path but the desired radius
    TU_temp = TU;
    TU_temp.radius = r(i);
    win_mask = ~isnan(TU_temp.window(:,1));
    
    % Handle case where there was no tube window.  Here, just make the
    % first element the desired tube radius
    if sum(win_mask) == 0
        win_mask(1) = true;
    end
    
    % Update radius in tube object window and re-calcualte boundary
    TU_temp.window(win_mask,1) = r(i);
    TU_temp = TU_temp.calcBoundary;
    
    % Loop over trials and apply tube constraint
    TDtemp = TD;
    for j = 1:n_trials
        pos = TDtemp(j).(srcStr).pos;
        contMask = TU_temp.calcContainment(pos,'plotSegment',plotSegment);

        % If not all points are contained in the tube, update the
        % trajectory data object
        if sum(contMask) ~= length(contMask)
            % Find index of first failed point
            offsetInd = find(~contMask,1,'first') - 1;

            % Truncate kinematics
            TDtemp(j).(srcStr).time = TDtemp(j).(srcStr).time(1:offsetInd);
            TDtemp(j).(srcStr).pos = TDtemp(j).(srcStr).pos(1:offsetInd,:);
            TDtemp(j).(srcStr).vel = TDtemp(j).(srcStr).vel(1:offsetInd,:);

            % Update success code
            TDtemp(j).successful = false;
        end
    end
    
    % Get mask with the success state for all trials
    sr_all{i} = logical([TDtemp.successful]);
    
    % Get expected success rate for radius
    sr(i) = sum(logical([TDtemp.successful]))/n_trials;
end