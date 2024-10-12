function [spikeTimes,nSpikes] = getSpikeTimes(S,tOffset,window)
% getSpikeTimes         Get spike time data method for SpikeData class
%
% Inputs
%   S               SpikeData object
%   tOffset         Time offset for each trial
%   window          Window to collect spike data over
%
% Outputs
%   spikeTimes      Cell array containing spike times
%   nSpikes         Number of spikes per cell
%
% Author:       Alan D. Degenhart
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

nTrials = length(S);

if nargin == 1      % offset and window have not been specified
    tOffset = 0;
    doWindow = 0;
elseif nargin == 2  % window has not been specified
    doWindow = 0;
else
    doWindow = 1;
end

% Determine whether or not to find the analysis window if empty
if (nargin == 3) && isempty(window)
    doWindow = 1;
end

% Create time offset vector if needed.  This allows the offset to be
% specified for each trial individually or globally for all trials.
if length(tOffset) == 1
    tOffset = repmat(double(tOffset),nTrials,1);
end

% Loop over trials and get number of spikes per trial and collect spike
% data in a cell array
nSpikes = nan(nTrials,1);
spikeTimes = cell(nTrials,1);
for i = 1:nTrials
    % Get spikes and correct for offset and window
    tempSpikes = double(S(i).spikeTimes);
    tempSpikes = tempSpikes - tOffset(i);
    
    % Window spike data if desired
    if doWindow
        if length(window)==1 % Assume that the window starts at tOffset
            window = [0 window];
        end
        winMask = tempSpikes >= window(1) & tempSpikes < window(2);
        tempSpikes = tempSpikes(winMask);
    end
    
    % Get number of spikes and put into cell arry
    nSpikes(i) = length(tempSpikes);
    spikeTimes{i} = tempSpikes;
end