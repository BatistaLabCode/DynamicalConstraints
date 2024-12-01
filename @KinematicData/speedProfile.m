% speedProfile      Method to calculate and plot speed profile for
% KinematicData class
%
% Usage:
%   [s,t,F] = K.speedProfile
%
% Inputs:
%   K             Kinematic Data
%
% Optional Arguments:
%   useMoveOnset    Align to movement onset
%
% Outputs:
%   s               Calculated speed
%   t               time
%   F               Figure that plots the speed profile
%
% Author:   Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [s,t,F] = speedProfile(K,varargin)

% Optional arguments
useMoveOnset = true;

assignopts(who,varargin);

% Get number of samples for each trial
nTrials = length(K);
nSamp = nan(nTrials,1);
tMin = nan(nTrials,1);
tMax = nan(nTrials,1);
tCell = cell(nTrials,1);
sCell = cell(nTrials,1);
for i = 1:nTrials
    % Get time and velocity, calculate speed
    t = K(i).time;
    v = K(i).vel;
    s = sqrt(v(:,1).^2 + v(:,2).^2);
    nSamp(i) = length(t);
    
    % Correct time for movement onset
    if useMoveOnset
        onsetIdx = K(i).moveOnset;
        tOnset = t(onsetIdx);
        t = t-tOnset;
    end
    
    % Interpolate over consistent time scale and add to cell arrays
    tInterp = t(1):t(end); % Use 1ms time grid
    sInterp = interp1(t,s,tInterp);
    tMin(i) = tInterp(1);
    tMax(i) = tInterp(end);
    tCell{i} = tInterp;
    sCell{i} = sInterp;
end

% Determine min and max time.  This is used to align the data for averaging
tMin = min(tMin);
tMax = max(tMax);
tAll = tMin:tMax;
maxSamp = length(tAll);

% Setup figure
nRow = 1;
nCol = 2;
axSp = 50;
axW = 300;
axH = 300;
[fW, fH, Ax] = calcFigureSize(nRow,nCol,axW,axH,axSp);

F = figure('Position',[10 10 fW fH]);
subplotSimple(nRow,nCol,1,'Ax',Ax); hold on;

% Loop over trials and get speed
sMat = nan(nTrials,maxSamp);
for i = 1:nTrials
    % Get data
    t = tCell{i};
    s = sCell{i};
    
    % Plot data
    plot(t,s,'color',ones(1,3)*.75)
    
    % Add data to matrix
    tMask = (tAll >= t(1)) & (tAll <= t(end));
    sMat(i,tMask) = s;
end

% Calculate average speed profile and standard error bars
nPts = sum(~isnan(sMat),1);
sMean = nanmean(sMat,1);
stdDev = nanstd(sMat,0,1);

% Find time limits.  Only show those time points where at least 50% of the
% trials are present
ptMask = nPts > (max(nPts)*.5);
tTemp = tAll(ptMask);
sTemp = sMean(ptMask);
%ci = stdDev(ptMask) ./ sqrt(nPts(ptMask));
ci = stdDev(ptMask);

% Plot mean
plot(tTemp,sTemp,'k-','LineWidth',2)

% Plot movement onset
yLim = get(gca,'YLim');
plot([0 0],yLim,'r--')

% Format plot
xlabel('Time (ms)')
ylabel('Speed (m/s)')
set(gca,'TickDir','out')
set(gca,'XLim',[tTemp(1) tTemp(end)],'YLim',yLim);
title('All trials')

% Plot mean + CI
subplotSimple(nRow,nCol,2,'Ax',Ax); hold on;
x = [tTemp fliplr(tTemp)];
y = [sTemp + ci fliplr(sTemp - ci)];
patch(x,y,ones(1,3)*.75,'EdgeColor','none');
plot(tTemp,sTemp,'k','LineWidth',2)

% Plot movement onset
plot([0 0],get(gca,'YLim'),'r--')

xlabel('Time (ms)')
ylabel('Speed (m/s)')
set(gca,'TickDir','out')
set(gca,'XLim',[tTemp(1) tTemp(end)]);
title('Confidence interval: StdDev')

% Set output variables.  The output will be only the segement of the speed
% profile where data for at least 50% of the trials is available.
t = tTemp;
s = sTemp;