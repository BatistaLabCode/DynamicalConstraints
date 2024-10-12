function raster(S,tOffset,window,varargin)
% Plot raster method for SpikeData class
%
% This function generates a raster plot for the spike time data in the
% input SpikeData vector
%
% Author:       Alan D. Degenhart
% Date Created: N/A
% Last Updated: 2016.06.08

% Optional arguments
channel = unique([S.channel]);  % Individual channels to plot
exclCh = [];                    % Channels to exclude from raster
validSrt = 1;                   % Valid sort codes

% Parse optional arguments
varRemain = assignopts('ignorecase', who, varargin);

% Ensure units are in ms
S = S.correctUnits;

% Get channel and sort codes and generate masks
ch = [S.channel];
srt = [S.sort];
exclChMask = ~ismember(ch,exclCh);
chMask = ismember(ch,channel);
srtMask = ismember(srt,validSrt);
dataMask = chMask & exclChMask & srtMask;

% Get rid of undesired data
S = S(dataMask);
ch = ch(dataMask);
srt = srt(dataMask);

% Check channel and sort information.  If there is more than one channel,
% assume a raster is being plotted for the entire array.  In this case,
% zero-pad up to the maximum channel ID.
uniCh = unique(ch);
uniSrt = unique(srt);
if (length(uniCh) > 1) && (length(uniSrt) == 1)
    % Data is from multiple channels with a single sort per channel.  In
    % this case, the y-axis can be the channel number.
    plotMode = 'channel';
    yAxStr = 'Channel';
elseif (length(uniCh) > 1) && (length(uniSrt) > 1)
    % Data is from multiple channels with multiple sorts per channel.  In
    % this case, y-axis is unit number and does not necessarily correspond
    % to the actual channel number
    plotMode = 'trial'; % Plotting raster for a single channel
    yAxStr = 'Unit';
else
    % Data is from a single trial.  In this case, the y-axis is the trial
    % number.
    plotMode = 'trial';
    yAxStr = 'Trial';
end

% Get all spike times for the provided spike data
[spikeTimes,nSpikes] = S.getSpikeTimes(tOffset,window);
nTrials = length(S);

% Create plotting data.  Collapse spike times into a single long array
% along with trial number.
xData = nan(1,sum(nSpikes));
yData = nan(1,length(xData));
onset = 1;
for i = 1:nTrials
    offset = onset + nSpikes(i) - 1;
    xData(onset:offset) = spikeTimes{i};
    
    % Set y-axis value.  If plotting trials, increment y-axis values
    % sequentially.  If plotting channels for a single sort, use the
    % channel number as the y-axis value.  For multiple sorts, use the
    % 'trial' mode and plot channel/sorts sequentially.
    switch plotMode
        case 'trial'
            yVal = i;
        case 'channel'
            yVal = ch(i);
    end
    
    yData(onset:offset) = yVal;
    onset = offset + 1;
end

yLim = [0 (max(yData) + 1)];

plot([xData;xData],[yData + 0.4; yData - 0.4],'k')
ylabel(yAxStr)
xlabel('Time (ms)')
set(gca,'YLim',yLim)