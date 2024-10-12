function [pAvg,p,t] = psth(S,varargin)
% Calculate PSTH method for SpikeData class.
%
% [p] = S.psth
%
% Optional arguments:
%
% Author:       Alan D. Degenhart
% Date Created: 2016.08.03
% Last Updated: 2016.08.03
% Last Update:  Initial version of function.

% Optional arguments
t0 = 0;                 % Time offset
preOnsetSamp = 0;       % Pre-onset samples
postOnsetSamp = 1000;   % Post-onset samples
sd = 20;        % Standard deviation of Gaussian kernel used for convolution
validSort = 1;  % Valid sort codes
exclCh = [];    % Channels to exclude
inclCh = [];    % Channels to include in the psth
plotFlag = false;   % Plot PSTH

assignopts(who, varargin);

% Design kernel
tKernel = -100:1:100;
k = normpdf(tKernel,0,sd)*1000;

% Remove unwanted channels and sort codes
srt = [S.sort];
ch = [S.channel];
srtMask = ismember(srt,validSort);
chMask = ~ismember(ch,exclCh) & (isempty(inclCh)+ismember(ch,inclCh));
mask = srtMask & chMask;
S = S(mask);
nTrials = length(S);

% Define time axis.  Currently this is set up such that t = 1 is the first
% time step.  This is consistent with the state timing in MonkeyHost.
t = (-preOnsetSamp + 1):(postOnsetSamp - 1);
%t = (t0 - preOnsetSamp + 1):(t0 + postOnsetSamp - 1);
nTS = length(t);

% Define t0 so that there is a unique value for each trial
if length(t0)==1
    t0 = t0*ones(nTrials,1);
end

% Loop over spike data objects and calculate PSTH
p = nan(nTrials,nTS);
r = zeros(1,nTS);   % "Empty" raster (pre-initialized to save time)
for i = 1:nTrials
    % Get spike times and remove those outside of the time range of
    % interest
    st = double(S(i).spikeTimes) - t0(i);
    stMask = (st >= t(1)) & (st <= t(end));
    st = st(stMask);
    
    % Correct for pre-onset samples.  This converts spike times into
    % indices for the time vector
    %st = st + preOnsetSamp + t0(i);
    st = st + preOnsetSamp;
    rTemp = r;
    rTemp(st) = 1;
    
    % Convolve to calculate PSTH
    r2(i,:) = rTemp;
    p(i,:) = conv(rTemp,k,'same');
end

% Calculate average PSTH
pAvg = nanmean(p,1);

if plotFlag
    % Calculate standard deviation
    s = std(p,[],1);
    uci = pAvg + s;
    lci = pAvg - s;
    
    x = [t fliplr(t)];
    y = [lci fliplr(uci)];
    
    hold on
    patch(x,y,ones(1,3)*.8)
    plot(t,pAvg,'k-','LineWidth',2)
end