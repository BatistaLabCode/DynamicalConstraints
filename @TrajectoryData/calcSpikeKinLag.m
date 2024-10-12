% calcSpikeKinLag       Calculate lag between spiking and kinematic data
% for TrajectoryData class.
%
% Optional Inputs:
%   taskStr    Task Information
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [lag,F] = calcSpikeKinLag(TD,varargin)

% Optional arguments
taskStr = '';

assignopts(who,varargin);

% Get dataset string (for plot titles and filenames)
dS = datestr(TD(1).date,'yyyymmdd');

% Plot speed profile
K = [TD.handKin];
[sP,t_sP,Fsp] = K.speedProfile;
plotDatasetInfo(TD(1).subject,TD(1).date,[TD.trialID])
titleStr = sprintf('%s%s :: %s :: Speed profile', ...
    TD(1).subject,dS,taskStr);
plotTitle(titleStr)
fileNameStr = sprintf('%s%s_SpeedProfile',TD(1).subject,dS);
set(Fsp,'Name',fileNameStr)

% Re-bin data starting at movement onset
bW = 45; % bin width (ms)
spOffset = 0;
winStart = 1;
nTrials = length(TD);
for i = 1:nTrials
    % Get movement onset time (used to re-bin spikes)
    t = TD(i).handKin.time;
    onsetIdx = TD(i).handKin.moveOnset;
    
    % Get spikes, combine all valid sort codes, and re-bin
    S = TD(i).spikes;
    S = S.combineSorts;
    t0 = t(onsetIdx);
    winEnd = TD(i).trajOffset;
    S = S.binSpikes(bW,t0,spOffset,[winStart winEnd]);
    
    % Update SpikeData object
    TD(i).spikes = S;
end

% Define time bins for combined data
t0 = 0;
window = [-500 1000]; % relative to movement onset
[binTimes,~,~] = generateBins(bW,t0,window);
nBins = length(binTimes);

% Loop over trials and get spike count data
nCh = 96;
validCh = 1:nCh;
validSort = 1;
spikeCountMat = nan(nTrials,nCh,nBins);
for i = 1:nTrials
    % Get spike data for current trial and define mask of valid
    % channels/sorts
    S = TD(i).spikes;
    spikeCh = [S.channel];
    chMask = ismember(spikeCh,validCh);
    sortMask = ismember([S.sort],validSort);
    spikeMask = chMask & sortMask;
    
    % Get movement onset and define time mask.  This should be the same for
    % all channels b/c they are all binned with respect to the same window.
    moveOnsetIdx = TD(i).handKin.moveOnset;
    moveOnset = TD(i).handKin.time(moveOnsetIdx);
    t = S(1).binTimes - moveOnset;
    tMask = ismember(t,binTimes);
    
    % Get spike counts (time x channel/sort) and add to matrix
    spikeCounts = [S.binSpikeCount];
    spikeCounts = spikeCounts(tMask,spikeMask);
    
    % Add spike counts to matrix.  Need to define a new time mask with
    % respect to the bin times for the combined data
    tMask = ismember(binTimes,t);
    spikeCountMat(i,:,tMask) = spikeCounts';
end

% Re-code data by target with peak modulation.  This is done in order to
% make averaging more effective.  For each electrode, averaging should only
% be done across those trials where the target was closest to the preferred
% direction of the unit.
tWin = 250; % Window to use (surrounding movement onset)
tC = [TD.targetCode];
uniTarg = unique(tC);
nTarg = length(uniTarg);
spikeCountTargMat = nan(nTarg,nCh,nBins);
peakSpikeCountMat = nan(nTarg,nCh);
tMask = (binTimes > -tWin/2) & (binTimes < tWin/2); % Focus on 500ms surrounding movement onset
for i = 1:nTarg
    % Get spike counts for current target and average
    targMask = ismember(tC,uniTarg(i));
    spikeCountTarg = spikeCountMat(targMask,:,:);
    spikeCountTarg = squeeze(nanmean(spikeCountTarg,1));
    
    % Find peak spike count over time
    peakSpikeCount = spikeCountTarg(:,tMask);
    peakSpikeCount = nanmean(peakSpikeCount,2);
    
    % Add data to large matrices
    peakSpikeCountMat(i,:) = peakSpikeCount;
    spikeCountTargMat(i,:,:) = spikeCountTarg;
end

% Loop over electrodes and re-align data such that the peak spike count
% across targets occurs for the same target code
peakTargSpikeCountMat = nan(nCh,nBins);
for i = 1:nCh
    % Get average perimovement spike count and find max
    targSpikeCount = peakSpikeCountMat(:,i);
    [~,I] = max(targSpikeCount);
    peakTargSpikeCountMat(i,:) = spikeCountTargMat(I,i,:);
    
    % Re-map spike counts so that the max occurs at the 4th index
    targInd = (I-3):I+4;
    targMask = targInd < 1;
    targInd(targMask) = targInd(targMask) + max(uniTarg);
    targMask = targInd > max(uniTarg);
    targInd(targMask) = targInd(targMask) - max(uniTarg);
    
    spikeCountTargMat(:,i,:) = spikeCountTargMat(targInd,i,:);
end

% Get peak spike counts.  By definition, this will correspond to the 4th
% target
peakSpikeCount = spikeCountTargMat(:,:,tMask);
peakSpikeCount = squeeze(nanmean(peakSpikeCount,3));
[~,Isrt] = sort(peakSpikeCount(4,:));
peakFR = peakSpikeCount/(bW/1000);

% Plot aligned tuning curves -- calculate average response for each neuron
% as a function of target direction, then re-map 

% Setup figure
nRow = 1;
nCol = 4;
axSp = 50;
axW = 300;
axH = 300;
[fW, fH, Ax] = calcFigureSize(nRow,nCol,axW,axH,axSp);

Fneural = figure('Position',[10 10 fW fH]);
fileNameStr = sprintf('%s%s_NeuralPopulationResponse',TD(1).subject,dS);
set(Fneural,'Name',fileNameStr)

subplotSimple(nRow,nCol,1,'Ax',Ax); hold on;

% Setup colors
cHSV_light = [0 1 0.9];
vLight = cHSV_light(3);
vDark = vLight * 0.1;
vAll = linspace(vDark,vLight,nCh);
cHSV_all = repmat(cHSV_light,nCh,1);
cHSV_all(:,3) = vAll;
c = hsv2rgb(cHSV_all);

x = ((1:nTarg) - 4)*(360/nTarg); % Currently assume targets are evenly-distributed
for i = 1:nCh
    % Get data for sorted channel index and tuning curve
    chInd = Isrt(i);
    frData = peakFR(:,chInd);
    plot(x,frData,'color',c(i,:))
end
plot(zeros(1,2),get(gca,'YLim'),'k--')

xlabel('Target direction (degrees, sorted)')
ylabel('Firing rate (sp/s)')
set(gca,'TickDir','out','XTick',x)
title('Aligned tuning curves sorted by peak firing rate')

% Plot population PSTH vs time -- order/color targets by target modulation.
% Find the target-sorted responses for each channel, then average across
% channels for each sorted target.  Thus, the PSTH with the lowest response
% will correspond to the target with the lowest firing rate for each
% neuron.

subplotSimple(nRow,nCol,2,'Ax',Ax); hold on;

% Setup colors
cHSV_light = [204/360 1 0.9];
vLight = cHSV_light(3);
vDark = vLight * 0.1;
vAll = linspace(vDark,vLight,nTarg);
cHSV_all = repmat(cHSV_light,nTarg,1);
cHSV_all(:,3) = vAll;
c = hsv2rgb(cHSV_all);

% Loop over targets
[~,Itarg] = sort(peakSpikeCount); % Sort responses for each target
for i = 1:nTarg
    % Get the appropriate target for each channel
    frTemp = nan(nCh,nBins);
    targIdx = Itarg(i,:);
    for j = 1:nCh
        scTemp(j,:) = squeeze(spikeCountTargMat(targIdx(j),j,:));
    end
    frTemp = nanmean(scTemp,1)/(bW/1000); % Average over neurons and convert to firing rate
    plot(binTimes,frTemp,'color',c(i,:),'LineWidth',2)
end

set(gca,'XLim',[-500 500],'TickDir','Out')
xlabel('Time relative to movement onset (ms)')
ylabel('Firing rate (sp/s)')
title('Population PSTH sorted by target')
yLim = get(gca,'YLim');

% Calculate average firing rate over time for the target condition with the
% highest firing rate.
frMat = peakTargSpikeCountMat/(bW/1000);
avgFR = nanmean(frMat,1);
frStdErr = nanstd(frMat,0,1)/sqrt(nCh);

% Truncate FR data to +/- 500ms
tMask = (binTimes > -500) & (binTimes < 500);
t_fR = binTimes(tMask);
avgFR = avgFR(tMask);
frStdErr = frStdErr(tMask);

% Plot average firing rate
subplotSimple(nRow,nCol,3,'Ax',Ax); hold on;

% Set up data for plotting CI patch
cPatch_HSV = cHSV_light;
cPatch_HSV(2) = .4;
cPatch_RGB = hsv2rgb(cPatch_HSV);

xPatch = [t_fR fliplr(t_fR)];
yPatch = [(avgFR + frStdErr) fliplr(avgFR - frStdErr)];
patch(xPatch,yPatch,cPatch_RGB,'EdgeColor','none')
plot(t_fR,avgFR,'color',c(end,:),'LineWidth',2)

set(gca,'XLim',[-500 500],'TickDir','Out')
xlabel('Time relative to movement onset (ms)')
ylabel('Firing rate (sp/s)')
title('Population PSTH (peak target -- CI:SE)')

% Truncate speed profile to +/- 500ms
tMask = (t_sP > -500) & (t_sP < 500);
t_sP = t_sP(tMask);
sP = sP(tMask);

avgFR = (avgFR - min(avgFR))/(max(avgFR) - min(avgFR));
sP = (sP - min(sP))/(max(sP) - min(sP));

% Calculate lag
[~,I_kin] = max(sP);
tMax_kin = t_sP(I_kin);
[~,I_neural] = max(avgFR);
tMax_neural = t_fR(I_neural);
lag = tMax_neural - tMax_kin;

subplotSimple(nRow,nCol,4,'Ax',Ax); hold on;

% Plot data
ax(1) = plot(t_fR,avgFR,'color',c(end,:),'LineWidth',2);
ax(2) = plot(t_sP,sP,'k','LineWidth',2);
yLim = get(gca,'YLim');
plot(ones(1,2) * tMax_neural,yLim,'LineStyle','--','color',c(end,:))
plot(ones(1,2) * tMax_kin,yLim,'k--')

ylabel('Normalized activity')
xlabel('Time (ms)')
title(sprintf('Lag (neural - kinematics): %0.3f ms',lag))
set(gca,'TickDir','out')
Hl = legend(ax,{'Neural','Kinematics'},'Location','NorthWest');
set(Hl,'box','off')

% Set figure title
plotDatasetInfo(TD(1).subject,TD(1).date,[TD.trialID])
titleStr = sprintf('%s%s :: %s :: Neural population response', ...
    TD(1).subject,dS,taskStr);
plotTitle(titleStr)

% Set output figure handle
F = [Fsp Fneural];