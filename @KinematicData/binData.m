function [k,binTimes] = binData(K,dataType,binWidth,t0,window,tOffset)
% binData           Bin kinematic data method for KinematicData class
%
% Inputs:
%   dataType        Type of kinematic data to bin (pos, vel, or acc)
%   binWidth        Width of each time bin
%   tOffset         Time offset to apply to data
%   window          Window to bin data over
%   
% Outputs:
%   k               Binned kinematic data
%   binTimes        Bin edges
%
% Author:   Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Verify that 'dataType' is appropriate
dataType = lower(dataType);
switch dataType
    case {'pos','position'}
        dataType = 'pos';
    case {'vel','velocity'}
        dataType = 'vel';
    case {'acc','acceleration'}
        dataType = 'acc';
    otherwise
        error('Invalid dataType specified.')
end

nTrials = length(K);

% Check size of offset.  This should be either a scalar or a vector the
% size of K
if length(tOffset) == 1
    tOffset = ones(nTrials,1) * tOffset;
elseif length(tOffset) ~= nTrials
    error('Length of tOffset inconsistent with the kinematic data provided.')
end

% Generate bin times
[binTimes,binOnset,binOffset] = util.generateBins(binWidth,t0,window);

% Initialize matrices
nBins = length(binTimes);
nKinDim = size(K(1).pos,2);
k = repmat(struct('pos',[],'vel',[],'acc',[]),nTrials,1);

% Loop over data and bin
for i = 1:nTrials
    % Get position data
    t = K(i).time;
    pos = K(i).pos;
    
    % Interpolate position
    posInterp = interp1(t,pos,binTimes);
    
    % Find velocity and acceleration by differentiating position and
    % velocity
    velInterp = diff(posInterp);
    accInterp = diff(velInterp);
    
    % Add data to struture
    k(i).pos = posInterp;
    k(i).vel = nan(size(posInterp));
    k(i).vel(2:end,:) = velInterp;
    k(i).acc = nan(size(posInterp));
    k(i).acc(3:end,:) = accInterp;
end