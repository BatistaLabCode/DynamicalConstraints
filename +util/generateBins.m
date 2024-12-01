%[binTimes,binOnset,binOffset] = generateBins(binWidth,t0,window,binType) 
% generateBins          Create bin arrays for binning data
%
% Usage:
%   util.generateBins(binWidth,t0,window,binType)
% 
% Inputs
%   binWidth            Width of each time bin
%   t0                  Time point to generate bins relative to
%   window              Range of bin times
%   binType             Type of bins ('causal' or 'centered')
%
% Outputs
%   binTimes            Time for each bin
%   binOnset            Time onset for each bin
%   binOffset           Time offset for each bin
%
% binType specifies how data should be binned.  Valid options are 'causal'
% (default), and 'centered'.  'Causal' specifies that the data at bin time
% t come from the interval [t-winWidth t]. 'Centered' specifies that data
% at bin time t come from the interval [t-binWidth/2 t+binWidth/2]. The
% default bin type is 'causal'.
%
% Bins are aligned at t = 0.  If a different type of alignment is desired,
% data should be time-shifted prior to binning.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [binTimes,binOnset,binOffset] = generateBins(binWidth,t0,window,binType)

if nargin < 4
    binType = 'causal';
end

% Create bin times
negWin = fliplr(t0:-binWidth:window(1));
posWin = (t0+binWidth):binWidth:window(2);
binTimes = [negWin posWin];

% Create bin edges
switch lower(binType)
    case 'centered'
        % Bins are centered at each bin time
        binOnset = binTimes - binWidth/2;
        binOffset = binTimes + binWidth/2;
    case 'causal'
        binOnset = binTimes - binWidth;
        binOffset = binTimes;
    otherwise
        error('Invalid bin type specified. Valid bin types are "centered" or "causal".')
end