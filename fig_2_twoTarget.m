function [F] = fig_2_twoTarget(dataLoc,varargin)
% [F] = fig_2_twoTarget(dataLoc) Plots the intuitive mapping (IM) 
% two-target trajectories for a given set of example sessions
%
% Inputs:
%   dataLoc    Paths for the main data folder
%
% Optional Inputs:
%   exampleSess     The example session used in the paper (fig2).
%   saveFig         Determine whether or not to save the data
%   savePathBase    Where to save the figures.
%   taskType        Which task type to plot
%   plotScale       Axis Limits 
%   centerPos       Center of the workspace
%   avgMode         Average method for the trajectories.
%   C               Colormap struct. Default bc_int
% 
% Outputs:
%   F               Figure with 4 panels of trajectories.
%
% Created by Erinn Grigsby
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

%% Load in the D structure and the data
exampleSess = {'20181004'}; % The example session used in the paper (fig2).
saveFig = 0;                % Determine whether or not to save the data,
                            % default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
taskType = {'tt_int_all'};  % Which task type to plot
plotScale = 180;               % Axis Limits plotScal*[-1 1 -1 1]
centerPos = [0 0 0];           % Defined center of the workspace
avgMode = 'samp';              % Average method for the trajectories.
C = [];                        % Colormap struct. Default bc_int
colIdxCor = [1 2 3 4 6 5 8 7]; % If you need to adjust the color scheme, 
                               % default is to match to figure 2 format. 
                               % If empty([]) colored in a clockwise
                               % manner.
pltIdx = [7 3;5 1;8 4; 2 6];   % If you need to adjust the plotting order,
                               % Default is to figure 2, empty plots from
                               % smallest target angle and its pair to
                               % largest.

% Assign the optional inputs
assignopts(who,varargin);

% Determine the sessions with the correct data
load(fullfile(dataLoc,'publicationQualitySessions.mat'));
D = D(ismember({D.dataset},exampleSess));
dir_list = db.get_task_datasets(D,taskType);

if isempty(dir_list)
    F = [];
    warning('Files were not found')
end

if isempty(C)
C = util.defineTaskColormap('bc_int');
end

for k = 1:size(dir_list,1)

    TD_int = util.loadSessionData(dir_list(k));

    % Normalize the data
    if ~ismember(mean(unique(round([TD_int.startPos]'),'rows')),centerPos)
        TD_int = TD_int.normalize(mean(unique([TD_rot.startPos]','rows')));
    else
        TD_int = TD_int.normalize(centerPos);
    end

    % Exclude any trials that were not part of the two target pair.
    uniStartTarg = unique([TD_int.startPos]','rows');
    TD_int = TD_int(ismember([TD_int.targPos]',uniStartTarg,'rows'));

    % Adjust the color target scale to match the targets
    scale = norm(TD_int(1).startPos(1:2))./100;
    C.targPos = C.targPos.*(scale/.9);

    % Adjust the color order to match with presentation in the figure
    if ~isempty(colIdxCor)
    C.col_dark = C.col_dark(colIdxCor,:);
    C.col_light = C.col_light(colIdxCor,:);
    end
    if isempty(pltIdx)
        pltIdx = (1:4)'+[0 4];
    end

    % Create a summary trajectory figures for the session.
    axSp = 75;
    axW = 300;
    [fW, fH, Ax] = plt.calcFigureSize(2,2,axW,axW,axSp);

    % Plot the trajectories
    F(k) = figure('Position',[10 10 fW fH]); % Trajectory plot
    F(k).Name = sprintf('%s%s_TwoTargetSummary',...
        dir_list(k,1).subject, dir_list(k,1).dataset);
    plt.plotTitle(sprintf("Summary Plot: %s %s",dir_list(k,1).subject,...
        dir_list(k,1).dataset));
    for n = 1:4
        figure(F(k))
        plt.subplotSimple(2,2,n,'Ax',Ax); % Set the new axes positions

        % Identify the target pairs
        m1 = ismember([TD_int.targetCode],pltIdx(n,1));
        m2 = ismember([TD_int.targetCode],pltIdx(n,2));

        % Calculate the average
        TDavg1 = TD_int(m1==1).average('avgMode',avgMode);
        TDavg2 = TD_int(m2==1).average('avgMode',avgMode);

        % Plot the individua; trajectories
        TD_int(m1==1).plot('ColMat',C,'trajColorWeight',.3,...
            'plotTargets',0,'markerSize',3);
        TD_int(m2==1).plot('ColMat',C,'trajColorWeight',.3,...
            'plotTargets',0,'markerSize',3);

        % Plot targets
        TDavg1.plot('ColMat',C,'plotTrajectories',0,'plotTargNum',0,...
            'plotStartTarg',1,'plotEndTarg',0,'targColMatIdx',2)
        TDavg2.plot('ColMat',C,'plotTrajectories',0,'plotTargNum',0,...
            'plotStartTarg',1,'plotEndTarg',0,'targColMatIdx',2)

        % Plot average trajectories
        TDavg1.plot('ColMat',C,'endMarker','arrow',...
            'plotTargets',0,'axesLimits',plotScale*[-1 1 -1 1])
        TDavg2.plot('ColMat',C,'endMarker','arrow',...
            'plotTargets',0,'axesLimits',plotScale*[-1 1 -1 1])

        plt.scaleBar(gca,20,'mm')
        axis on,set(gca,'XTick',[],'YTick',[])
    end
end

%% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF(F,savePathBase)
end

end