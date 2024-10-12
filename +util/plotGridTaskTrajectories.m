function F = plotGridTaskTrajectories(TD,fName,varargin)
% Analysis script for position Kalman filter grid task
%
% Inputs:
%   TD      TrajectoryData object 
%
% Author:           Alan D. Degenhart
% Date Created:     2015/07/07
% Date Updated:     2016/08/05
% Last Update:      Function now takes TrajectoryData object as an input

centerPos = [];
plotIntTarg = false;
plotScaleFactor = 1.5;
plotScale = [];
plotStates = {'Step 2'};
ColMat = [];
axisOn = 1;
avgMode = 'spatial';
startMarker = 'o';
endMarker = 'o';
trlMarkerSz = 5;
avgMarkerSz = 10;
plotTargets = 1;

assignopts(who,varargin);

% Remove failed trials
sc = logical([TD.successful]);
TD = TD(sc);

% Get all start and end positions
startPos = [TD.startPos]';
endPos = [TD.targPos]';
allPos = [startPos;endPos];

% Calculate center position if not provided.  Calculate center position as
% the average across unique target positions.
if isempty(centerPos)
    uniPos = unique(allPos,'rows');
    centerPos = mean(uniPos); 
end

% Subtract off center position to normalize.  Use of the 'repmat' function
% is no longer necessary in newer (~R2016) versions of MATLAB, but is being
% kept here for compatibility with older versions.
allPos = allPos - repmat(centerPos,size(allPos,1),1);

% Calculate distance from center.  This is used to normalize the plot axes
% appropriately.
if isempty(plotScale)
    r = sqrt(allPos(:,1).^2 + allPos(:,2).^2);
    rMax = max(r) + TD(1).targSize; % Add target radius
    plotScale = rMax*plotScaleFactor;
end

% Get target codes for start and end targets.  Note that the target code
% has to be determined using the combined set of start and end target
% positions, 
tcAll = util.pos2targetCode(allPos,'centerPos',centerPos);
tcStart = tcAll(1:size(startPos,1));
tcEnd = tcAll(size(startPos,1) + 1:end);
[uniTC_all,IA] = unique(tcAll);
uniTargPos_all = allPos(IA,:); % Mapping from target code to position

% Find all target pairs.  Here we assume that the task is *always* an
% AB<->BA task.
targComb = [tcStart tcEnd];
uniTargComb = unique(targComb,'rows');
targPairs = sort(uniTargComb,2);
targPairs = unique(targPairs,'rows');

% Setup plot colors
col = [0 .6 .85;.8 .15 0];
colAvg = [0 .6 .85;.8 .15 0];
colGray = ones(1,3)*.75;

% Determine plot arrangement
nPairs = size(targPairs,1);
nCol = 3;
nRow = nPairs;

% Setup figure
axSp = 75;
axW = 300;
[fW, fH, Ax] = calcFigureSize(nRow,nCol,axW,axW,axSp);

fCtr = 1;
F(fCtr) = figure('Position',[10 10 fW fH]);
set(F(fCtr),'Name',[fName '_TargetPairs'])

plotTitle([fName ': Target Pairs'])
plt.plotDatasetInfo(TD(1).subject,TD(1).date,[TD.trialID])

% Normalize trajectory data
TDnorm = TD.normalize(centerPos);

% Loop over target pairs
for i = 1:nPairs
    % Get target subsets from A->B and B->A.  d1 represents trajectories
    % from A to B, and d2 representes trajectories from B to A
    d1Mask = (tcStart == targPairs(i,1)) & (tcEnd == targPairs(i,2)); % to *B*
    d2Mask = (tcStart == targPairs(i,2)) & (tcEnd == targPairs(i,1)); % to *A*
    startTarg_d1 = uniTargPos_all(targPairs(i,1),:);
    endTarg_d1 = uniTargPos_all(targPairs(i,2),:);
    startTarg_d2 = uniTargPos_all(targPairs(i,2),:);
    endTarg_d2 = uniTargPos_all(targPairs(i,1),:); 
    
    % Get target subset and set target code data.  This target code
    % information is only used to select the colors to use when plotting.
    TD_d1 = TD(d1Mask);
    TD_d1 = TD_d1.setTargetCode(1);
    TD_d2 = TD(d2Mask);
    TD_d2 = TD_d2.setTargetCode(2);
    
    % Merge data and normalize
    TD_merged = [TD_d1 ; TD_d2];
    
    % Get corresponding colors
    colInfo_d1 = util.getColorInfo(ColMat,startTarg_d1,endTarg_d1);
    colInfo_d2 = util.getColorInfo(ColMat,startTarg_d2,endTarg_d2);
    cL = [colInfo_d1{1};colInfo_d2{1}];
    cD = [colInfo_d1{2};colInfo_d2{2}];
    
    % Re-set TC information
    d1Inds = 1:length(TD_d1);
    d2Inds = (length(TD_d1)+1):length(TD_merged);
    all_inds = {d1Inds, d2Inds};

    % Plot trajectories for a single direction
    for j = 1:2
        plot_no = (i-1) * nCol + j;
        subplotSimple(nRow,nCol,plot_no,'Ax',Ax);
        inds_j = all_inds{j};
        TD_j = TD_merged(inds_j);
        
        % Plot targets -- plot both start and end target
        TD_j.plot('color',colGray,'plotTrajectories',0,...
            'plotStates', plotStates, ...
            'plotIntTarg',plotIntTarg, ...
            'plotStartTarg', plotTargets,...
            'plotTargets',plotTargets,...
            'startMarker',startMarker,...
            'endMarker',endMarker)

        % Plot targets and all trajectories.  Trajectories are colored by start
        % target.
        TD_j.plot('color',cL(j,:),'trajColorWeight',1,...
            'plotStates',plotStates,...
            'startMarker',startMarker,...
            'plotTargets',plotTargets,...
            'endMarker',endMarker,...
            'markerSize',trlMarkerSz)

        % Plot average trajectories.  Trajectories are colored by start target.
        TD_merged_avg = TD_j.average('avgMode', avgMode);
        TD_merged_avg.plot('color',cD(j,:),'plotTargets',0,...
            'lineWidth',2,...
            'startMarker',startMarker,...
            'endMarker',endMarker,...
            'markerSize',avgMarkerSz)

        axis([-1 1 -1 1]*plotScale)
        plt.scaleBar(gca,20,'mm')
        axis on
        set(gca,'XTick',[],'YTick',[])

        % Plot number of trials per trajectory.  Text is colored by start
        % target.
        text(-1*plotScale*.95,-1*plotScale*.9,sprintf('N = %d',length(inds_j)), ...
            'color',cD(j,:),'FontSize',14)
    end
    
    % Plot all targets (across target pairs).  This will label the targets
    % by their target codes, which are defined as the unique target
    % positions.
    plot_no = (i-1) * nCol + 3;
    subplotSimple(nRow,nCol,plot_no,'Ax',Ax);
    TDnorm.plot('color',colGray,'plotTrajectories',0,...
        'plotStates',plotStates,...
        'plotTargets',plotTargets,...
        'plotIntTarg',plotIntTarg,...
        'startMarker',startMarker,...
        'endMarker',endMarker)

    % Plot targets and all trajectories.  Trajectories are colored by start
    % target.
    TD_merged.plot('color',cL,'trajColorWeight',1,...
        'plotStates',plotStates,...
        'plotTargets',plotTargets,...
        'markerStyle','none',...
        'startMarker',startMarker,...
        'endMarker',endMarker,...
        'markerSize',trlMarkerSz)
    
    % Plot average trajectories.  Trajectories are colored by start target.
    TD_merged_avg = TD_merged.average('avgMode', avgMode);
    TD_merged_avg.plot('color',cD,'plotTargets',0,'lineWidth',3,...
        'startMarker',startMarker,...
        'endMarker',endMarker,...
        'markerSize',avgMarkerSz)
    
    axis([-1 1 -1 1]*plotScale)
    plt.scaleBar(gca,20,'mm')
    if axisOn
        axis on
    else
        axis off
    end
    set(gca,'XTick',[],'YTick',[])
    
    % Plot number of trials per trajectory.  Text is colored by start
    % target.
    text(-1*plotScale*.95,-1*plotScale*.9,sprintf('N = %d',length(d1Inds)), ...
        'color',cD(1,:),'FontSize',14)
    text(-1*plotScale*.95,-1*plotScale*.8,sprintf('N = %d',length(d2Inds)), ...
        'color',cD(2,:),'FontSize',14)
end
