function [S] = correctUnits(S)
% [S] = correctUnits(S)
%
% Convert spike time to ms.
%
% This function converts the spike times in the provided SpikeData object
% into ms.  Generally, spike times are provided in units of ms, but
% occasionally higher precision is desired.  This function can then be used
% to keep the spike time data consistent with what is expected by other
% functions.
%
% Author:       Alan D. Degenhart
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Check units.  If ms, do nothing.
tsUnits = S(1).timeStampUnits; % Currently, assume all objects have the same units
switch tsUnits
    case 'ms'
        % Units are already ms, nothing needs to be done.  Return.
        return
    case 'us'
        % Need to convert spike times to ms.  Allow function to continue
    otherwise
        error('Invalid spike timestamp units provided.')
end

% Loop over spike time data and convert to ms
nSpikeData = length(S);
for i = 1:nSpikeData
    % Convert from us to ms
    spikeTimes = double(S(i).spikeTimes);
    spikeTimes = spikeTimes/1000;
    
    % Update SpikeData
    S(i).timeStampUnits = 'ms';
    S(i).spikeTimes = spikeTimes;
end