function [F,h,p] = fig_4_early_vs_late_two_target_trajectories_comparisons(dataLoc,varargin)
% [F,h,p] = fig_4_early_vs_late_two_target_trajectories_comparisons(dataLoc) 
% Plots the histograms of early vs late d' values for all sessions of the
% separation-maximizing two target task and shuffled calculation of d' for
% all session. Also plots example trajectories for the early two-target 
% trials and the late two-target trials from a single session.
%
% Inputs:
%   dataLoc    Paths for the main data folder
%
% Optional Inputs:
%   exampleSess         The example session used in the paper (fig4).
%   saveFig             Determine whether or not to save the data
%   savePathBase        Where to save the figures.
%   plotScale           Axis Limits
%   colShuff            Color map for the shuffle data
%   colSM               Color map for the early vs late data
%   markerSize          Marker size
%   condIdx             Condition that we are running the comparison on:
%                         There are two possible pairings:
%                         1) IM trial - IM proj
%                         2) SM trial - SM proj [Default]
%
% Outputs:
%   F                   Figures for the histograms and the example trials
%   h                   Hypothesis rejection value for comparing the 
%                           shuffle and early vs late sessions.
%   p                   P-value for comparing the shuffle and early vs
%                           late sessions.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Create a comparison of the early vs late trajectories for the rotated
% mapping, two target task.
exampleSess = {'20190719'};% Figure 4 example session
saveFig = 0;                % Determine whether or not to save the data,
                            % default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
plotScale = 180;            % Axis scale for the plot panels
colShuff = [180 178 180]./255; % Color map for the shuffle data
colSM = [219 240 245]./255;    % Color map for the early vs late data
markerSize = 2;                % Marker size
condIdx = 2; % The condition that we are running the comparison for, there 
             % are two possible pairings, 1) IM trial - IM proj and
             % 2) SM trial - SM proj. The default for this function (2).

% Assign the optional inputs
assignopts(who,varargin);

load(fullfile(dataLoc,'dPrime_EarlyLate+ShuffledData.mat'))

% Determine the sessions with the correct data
load(fullfile(dataLoc,'exampleDatasetCatalog.mat'))
dir_list = db.get_task_datasets(D, {'tt_int','tt_rot'});

% Plot the d' values for each condition (early, late, shuff1, and shuff2)
maskExSess = find(ismember({dir_list(:,1).dataset},exampleSess));

% Collect the d' information 
dp_E = [dat_dPrime_EL(:,condIdx).dP_E];
dp_L = [dat_dPrime_EL(:,condIdx).dP_L];
dp_sh1 = [dat_dPrime_EL(:,condIdx).dP_shuff1];
dp_sh2 = [dat_dPrime_EL(:,condIdx).dP_shuff2];

% Create a histogram of all animals
F = figure('Position',[112 50 1200 750]);hold on
addLeg = 0;

% Plot the difference (Early - Late)
histogram(dp_E-dp_L,'BinWidth',0.25,'FaceColor',colSM)
plot(mean(dp_E-dp_L,'omitnan'),0,'kv','MarkerFaceColor',colSM)
[hEL,pEL]= ttest(dp_E-dp_L);

% Plot the difference (Shuffle 1 - Shuffle 2)
histogram(dp_sh1-dp_sh2,'BinWidth',0.25,'FaceColor',colShuff)
plot(mean(dp_sh1-dp_sh2,'omitnan'),0,'kv','MarkerFaceColor',colShuff)
[hS12,pS12]= ttest(dp_sh1-dp_sh2);

% Plot the example sessions
if ~isempty(maskExSess)
    plot(dp_E(maskExSess)-dp_L(maskExSess),0,'ko','MarkerFaceColor',colSM)
    plot(dp_sh1(maskExSess)-dp_sh2(maskExSess),0,'ko','MarkerFaceColor',colShuff)
    addLeg = 1; % Add to the legend
end

% Set the plot details
title("d' difference for (Early - Late)/(Shuffle 1-2)")
xlabel("d'_{Early} - d'_{Late}"),ylabel("# Experiments")
if addLeg
    legend(sprintf('Early - Late, h=%0.0f,p=%d',hEL,pEL),...
    sprintf('Early - Late, Mean= %0.3f, std = %0.3f',...
        mean(dp_E-dp_L,'omitnan'),std(dp_E-dp_L,'omitnan')),...
    sprintf('Shuffle 1 - 2, h=%0.0f,p=%d',hS12,pS12),...
    sprintf('Shuffle 1 - 2, Mean= %0.3f, std = %0.3f',...
        mean(dp_sh1-dp_sh2,'omitnan'),std(dp_sh1-dp_sh2,'omitnan')),...
    'ExSess: Early - Late','ExSess: Shuffle 1 - 2')
else
    legend(sprintf('Early - Late, h=%0.0f,p=%d',hEL,pEL),...
    sprintf('Early - Late, Mean= %0.3f, std = %0.3f',...
        mean(dp_E-dp_L,'omitnan'),std(dp_E-dp_L,'omitnan')),...
    sprintf('Shuffle 1 - 2, h=%0.0f,p=%d',hS12,pS12),...
    sprintf('Shuffle 1 - 2, Mean= %0.3f, std = %0.3f',...
        mean(dp_sh1-dp_sh2,'omitnan'),std(dp_sh1-dp_sh2,'omitnan')))
end
ylim([0 20])
plt.plotTitle('Early vs Late comparison, single run')
F.Name = 'dPrime_early_v_late_singleRun_histogram';

% Compare the shuffle and early vs late sessions.
[h,p] = ttest(dp_E-dp_L,dp_sh1-dp_sh2);

% Create a plot of the example session
% Load in the trajectory data of the example session
[TD, ~, ~] = util.loadSessionData(dir_list(maskExSess,2));
TD = TD.normalize([0 0 0]);
TD = TD(ismember([TD.targPos]',unique([TD.startPos]','rows'),'rows'));

colMat = [0.06565 0.3067 0.41; 0.44	0.0572	0.06996];
targPair = [[TD.startPos]' [TD.targPos]'];
[uniTarg,~,idxTarg] =unique(targPair,'rows');

F(2) = figure('Position',[65 80 1200 675]);
F(2).Name = [exampleSess{1} '_early_vs_Late_example'];
avgTD = TD.average();

for n = 1:size(uniTarg,1)
    tmpTD = TD(idxTarg==n);
    szTD = length(tmpTD);
    midPt = round(szTD/2);

    subplot(1,2,1)
    title('Early Trajectories: A $\Longleftrightarrow$ B','Interpreter','latex')
    if n == 1
        avgTD.plot('plotTrajectories',0,'plotTargets',1,'plotTargNum',0,...
            'plotStartTarg',1,'plotEndTarg',1,'color',.9*[1 1 1])
    end
    tmpTD(1:midPt).plot('plotTargets',0,'color',colMat(n,:),...
        'trajColorWeight',.65,'startMarker',[],'markerSize',markerSize)

    subplot(1,2,2)
    title('Late Trajectories: A $\Longleftrightarrow$ B','Interpreter','latex')
    if n == 1
        avgTD.plot('plotTrajectories',0,'plotTargets',1,'plotTargNum',0,...
            'plotStartTarg',1,'plotEndTarg',1,'color',.9*[1 1 1])
    end
    tmpTD(midPt+1:end).plot('plotTargets',0,'color',colMat(n,:),...
        'startMarker',[],'markerSize',markerSize)
end

% Set the plot details
for n = 1:2
    subplot(1,2,n)
    axis([-1 1 -1 1]*plotScale)
    plt.scaleBar(gca,20,'mm')
    axis on
    set(gca,'XTick',[],'YTick',[])
end
plt.matchAxis(F(2));
plt.plotTitle(['Example Session: ' exampleSess])

% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF([F(:)],savePathBase)
end
end