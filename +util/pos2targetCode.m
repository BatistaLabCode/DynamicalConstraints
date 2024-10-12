function [tc,M] = pos2targetCode(pos,varargin)
% Convert target positions to unique target codes
%
% Inputs
%   pos         N x nDim matrix of position data
%
% Outputs
%   tc          Target code identifier vector

centerPos = [];
dim = 2;
assignopts(who,varargin);

nTrials = size(pos,1);

% If 2D set the 3rd dimension to zero
if dim == 2
    pos(:,3) = 0;
end
% Find unique targets and re-center
uniTarg = unique(pos,'rows');

% Calculate center position if not specified.
if isempty(centerPos)
    centerPos = mean(uniTarg,1,'omitnan');
end
pos = pos - repmat(centerPos,nTrials,1);

% Convert to radians (to keep target codes ordered appropriately).  Also
% compute distance from center so that targets with the same angle but
% different distances from the center are assigned unique target codes.
ang = atan2d(pos(:,2),pos(:,1));
ang = round(ang);
r = round(sqrt(pos(:,1).^2 + pos(:,2).^2),2);

% Adjust target angles such that theta = 0 is in the positive x-direction.
% This should allow target codes to be compared across datasets with the
% same target orientation relative to the center target.
ang(ang < 0) = ang(ang < 0) + 360;

% Find unique target
[~,~,IC] = unique([ang r],'rows');

% Loop over trajectories and update target code data
nTrials = size(pos,1);
tc = nan(nTrials,1);
for i = 1:nTrials
    % Set unique target code
    tc(i) = IC(i);
end

% Creat position/angle/tc map.  This is useful to resolve target locations
% from the target code information.
pos = pos + repmat(centerPos,nTrials,1);
m = [tc ang pos];
uniM = unique(m,'rows');
M.uniTC = uniM(:,1);
M.uniAng = uniM(:,2);
M.uniPos = uniM(:,3:end);