% calcMoveOnset     Calculate movement onset for kinematic data
%
% TD = TD.calcMoveOnset
%
% Author:   Alan D. Degenhart
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [K] = calcMoveOnset(K,varargin)
% Optional arguments
plotFlag = false;   % Plot speed profile and estimated movement onset
speedThresh = 0.5;  % Speed threshold
dim = [1 2];        % Dimensions to use for speed calculation
fCutoff = 5;        % Cutoff frequency for LP filter
minOnsetTime = 0;   % Minimum onset time

assignopts(who,varargin);

% Check to see if kinematic data is empty
if isempty(K(1).time)
    return
end

% Determine sampling rate and design filter
sR = mean(diff(K(1).time));
Fs = 1000/sR;
Wn = fCutoff/(Fs/2);
[b,a] = butter(4,Wn,'low');

% Make sure minimum onset time is the appropriate size
nTrials = length(K);
if length(minOnsetTime) == 1
    minOnsetTime = ones(1,nTrials)*minOnsetTime;
end
assert(length(minOnsetTime) == nTrials,'Size of provided minimum onset times does not match that of the kinematic data');

% Set up figure if desired
if plotFlag
    figure; hold on;
end

% Loop over trials and calculate movement onset

sMax = nan(nTrials,1);
for i = 1:nTrials
    % Get kinematic data
    v = K(i).vel;
    t = K(i).time;
    
    % If kinematic data is empty, skip
    if isempty(v) | length(v)==1
        continue
    end
    
    % Calculate speed and smooth
    T = length(t);
    s = nan(1,T);
    for j = 1:T
        s(j) = norm(v(j,dim));
    end
    
    % Remove the nan's from the structure, this will deal with marker
    % drops.
    idx = ~isnan(s);
    s = s(idx);
    t = t(idx);
    T = length(t);
    
    % Check to ensure sufficient samples exist to filter. If not, use
    % unsmoothed speed.
    if T > 12 % Currently hard-code sample limit
        s = filtfilt(b,a,s);
    else
        warning('\nInsufficient samples exist to smooth kinematic data.  Using raw speed profile.')
    end
    
    % Set speeds for those time points less than the minimum onset time to
    % 0.  This is used to ensure that speed peaks ocurring before a
    % specific time (e.g., a go cue) are ignored.
    s(t<minOnsetTime(i)) = 0;
    
    % Find point where speed exceeds threshold
    sMax(i) = max(s);
    [sMax,I] = find(s > speedThresh*sMax(i),1,'first');
    K(i).moveOnset = I;
    
    % Plot speed profile and movement onset if desired
    if plotFlag
        t = t - t(I); % Align to estimated movement onset
        plot(t,s,'k')
        plot(t(I),s(I),'ro')
    end
end

% Set axes on figure
if plotFlag
    xlabel('Time (ms)')
    ylabel('Speed (mm/s)')
    set(gca,'TickDir','out','Box','off')
end