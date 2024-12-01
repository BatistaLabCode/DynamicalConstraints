% F = plotIntermediateTargetTaskTrajectories(TD,C,fName,varargin)
% 
% Plot cursor trajectories for intermediate target task
%
% Usage:
%   F = util.plotIntermediateTargetTaskTrajectories(TD,C,fName)
%
% Inputs:
%   TD      TrajectoryData object
%   C       Colormap struct
%   fName   Name for title and figure ID
%
% Optional Inputs:
%   createFigure        Create a figure
%   plotIntTarg         Plot the intermediate target
%   plotStates          Determine which task states to plot
%   plotScale           Axis Limits
%   plotTubes           Plot the tube boundary
%   plotSegment         Plot only portions of the tube
%   avgMode             Average method for the trajectories
%   r_factor            Factor used for spatial averaging
%   trialsPerCondition  Number of trials to subselect for plotting
%   setRandSeed         Fix the random number generator for a consistenet
%                           permutation for trial subselection.
%   rmvFail             Remove failed trials
%
% Outputs:
%   F               Figure with 4 panels of trajectories.
%
% Author: Alan Degenhart and Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function F = plotIntermediateTargetTaskTrajectories(TD,C,fName,varargin)

createFigure = true;
plotIntTarg = false;
plotStates = {'Step 2'};
plotScale = [];
plotTubes = false;
plotSegment = [];
avgMode = 'spatial';
r_factor = 0.1;  % Factor used for spatial averaging
trialsPerCondition = [];     % Allow for subselection of trials
setRandSeed = [];
rmvFail = 1;

assignopts(who,varargin);

% Get all start and end positions.  Data should already be normalized
startPos = [TD.startPos]';
uniStartPos = unique(startPos,'rows');
nUniStartPos = size(uniStartPos,1);

% Setup plot colors
colGray = ones(1,3)*.75;

% Set up figure
if createFigure
    nCol = nUniStartPos;
    nRow = 1;
    axSp = 30;
    axW = 300;
    [fW, fH, Ax] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
    F = figure('Position',[100 100 fW fH]);
    set(F,'Name',[fName '_TargetPairs'])
    plt.plotTitle([fName ': Target Pairs'])
    plt.plotDatasetInfo(TD(1).subject,TD(1).date,[TD.trialID])
end

% Loop over target pairs
for i = 1:nUniStartPos
    % Get mask for start position
    startPos_temp = uniStartPos(i,:);
    posMask = ismember(startPos,startPos_temp,'rows');
    TD_temp = TD(posMask);

    % Get color
    [colInfo] = util.getColorInfo(C,startPos_temp,[]);
    cD = colInfo{2};
    cL = colInfo{1};

    % Remove failed trials
    if rmvFail
        scMask = logical([TD_temp.successful]);
    else
        scMask = ones(size(TD_temp));
        scMask = scMask == 1;
    end

    % Plot trajectories (success only)
    if sum(scMask) > 0
        TD_temp = TD_temp(scMask);

        % Create subplot
        if createFigure
            plt.subplotSimple(nRow,nCol,i,'Ax',Ax);
        end

        % Plot tubes (if desired)
        hold on;
        if plotTubes
            TU = [TD_temp.tube];
            TU(1).plot('lineSpec','--','plotSegment',plotSegment)
        end

        % Plot only the number of subselected trials
        if ~isempty(trialsPerCondition)
            if length(TD_temp)<trialsPerCondition
                pos = 1:length(TD_temp);
            else
                if ~isempty(setRandSeed)
                    rng(setRandSeed)
                    pos = randperm(length(TD_temp),trialsPerCondition);
                else
                    pos = randperm(length(TD_temp),trialsPerCondition);
                end
            end
        else
            pos = 1:length(TD_temp);
        end

        % Plot all targets (across target pairs)
        TD_temp(pos).plot('color',colGray,'plotTrajectories',0,...
            'plotStates',plotStates,'plotIntTarg',plotIntTarg, ...
            'plotTargNum',false)

        % Plot targets and all trajectories
        TD_temp(pos).plot('color',cL,'trajColorWeight',1,...
            'plotStates',plotStates,'markerStyle','none','plotStartTarg',true,...
            'plotEndTarg',false)

        % Plot average trajectories
        TD_temp = TD_temp.average('avgMode', avgMode, ...
            'r_factor', r_factor);
        TD_temp.plot('color',cD,'plotTargets',0,'lineWidth',3)
    end

    if ~isempty(plotScale)
        axis([-1 1 -1 1]*plotScale)
    end
    plt.scaleBar(gca,20,'mm')
    axis on
    set(gca,'XTick',[],'YTick',[])

    % Plot number of trials per trajectory.
    text(-1*plotScale*.95,1*plotScale*.1*(9.75),sprintf('N = %d/%d',...
        sum(scMask),length(scMask)),'color',cD,'FontSize',14, ...
        'VerticalAlignment','Top')
end