function [containMask] = calcContainment(TU,pos,varargin)
% Calculates whether a tube fully contains a trajectory or not. If the 
% trajectory is not fully contained, state the point at which the
% containment failed. 
%
% Usage:
%   TU = TU.calcContainment(pos)
%
% Input:
%   TU      TubeOject
%   pos     Position of a data point
%
% Optional Inputs:
%   plotSegment    Determine we are plotting a segment of the tube around
%                      the path or the entire tube around the path
%
% Copyright (C) by Batista Lab
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Define default values for the optional arguments.
plotSegment = [];

% load the optional arguments
assignopts('ignorecase',who,varargin);

% Determine if we are just plotting a segment of the tube.
if ~isempty(plotSegment)
    pathIdx = round(linspace(1,size(TU.path,1),TU.maxSegment+1));
    if plotSegment == 0
        plotSegment = 1;
    end
    TU.path = TU.path(pathIdx(plotSegment):pathIdx(plotSegment+1),:);
    TU = TU.calcBoundary;
end

tubePath = TU.path(:,1:2); % Currently only valid for 2D
pos = pos(:,1:2);

% Loop over positions to evaluate
nPos = size(pos,1);
containMask = false(nPos,1);
for i = 1:nPos
    % Find the shortest distance from the current position to the center of
    % the tube
    d = tubePath - pos(i,:);
    d = sqrt(d(:,1).^2 + d(:,2).^2);
    [dMinPos,I] = min(d);
    
    % Find the shortest distance from the identified tube path point to all
    % boundary points
    pathPos = tubePath(I,:);
    d = TU.boundary - pathPos;
    d = sqrt(d(:,1).^2 + d(:,2).^2);
    dMinBoundary = min(d);
    
    % If the distance to point to evaluate is greater than the distance to 
    % the boundary, flag the point as outside of the tube
    if dMinBoundary > dMinPos
        containMask(i) = true;
    end
end

% % DEBUGGING -- Plotting
% figure; hold on;
% TU.plot
% plot(pos(containMask,1),pos(containMask,2),'g.')
% plot(pos(~containMask,1),pos(~containMask,2),'rx')