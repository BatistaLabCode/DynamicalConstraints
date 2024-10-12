function [containsTD] = calcTubeContainment(TU,TD,varargin)
% Calculates whether a tube fully contains a trajectory or not. If the 
% trajectory is not fully contained, state the point at which the
% containment failed. 
%
% Note: This code currently assumes that the tube is fully enclosed.
%
% Usage:
%   TU = TU.calcTubeContainment(TD)
%
% p - tube path
% r - tube radius
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Define default values for the optional arguments. Currently, we don't have
% any optional input.

% load the optional arguments
assignopts('ignorecase',who,varargin);

p = TU.path(:,1:2); % Currently only valid for 2D

% Set radius from window.  Currently assume that the first non-NAN value in
% the 'window' is used for the entire trial.
w = TU.window(:,1);
w = w(~isnan(w));
if length(unique(w)) > 1
    warning('More than one unique tube tolerance window found. Using the first identified window.')
end
r = w(1);

% Initize structure
containsTD = struct('trialID',{},'successful',{},'traj',{},...
    'failedBoundaryPts',{},'distFromTUPath',{});

% Iterate through the trials
ntrials = length(TD);
for j = 1:ntrials
    % Find trial trajectory
    traj = TD(j).pos(:,1:2); % Just 2D
    
    % Set temporary values
    successful = 1;
    failBoundaryPts = [];
    calcDist = cell(0);
    
    % Calculate the distance between the tube path and the trajectory path.
    % If the minimum distance between the two is greater than the tube
    % radius then save the index point as a failure point.
    for k = 1:size(traj,1)
        calcDist{k} = sqrt(sum((traj(k,:)-p).^2,2));
        if min(calcDist{k})>r
            successful = 0;
            failBoundaryPts = [failBoundaryPts; k];
        end
    end
    
    % Save structure values
    containsTD(j).trialID = TD(j).trialID;
    containsTD(j).successful = successful;
    containsTD(j).failedBoundaryPts = failBoundaryPts;
    containsTD(j).distFromTUPath = calcDist';
    containsTD(j).traj = traj;
end
