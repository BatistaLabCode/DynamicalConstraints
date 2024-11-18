function [F2d,F,h2di,p2di] = fig_3(dataLoc,varargin)
% [F2d,F,h2di,p2di] = fig_3(dataLoc) Plots the intuitive mapping (IM) 
% two-target trajectories in the separating maximum project (SM). The code 
% also assess the d' in the workspace dimension for all session.
%
% Inputs:
%   dataLoc    Paths for the main data folder
%
% Optional Inputs:
%   exampleSess     The example session used in the paper (fig2).
%   saveFig         Determine whether or not to save the data
%   savePathBase    Where to save the figures.
%   colorMode       Trajectory color mode to use
%   plotPreOnset    Plot pre-onset data
%   xSpec           Plot the orthonormalize GPFA/latent space
%   markerScale     Scale of the start and end trajectory markers
%   markSize        Size of the midpoint markers (ind trls and mean)
%   plotDim         Dimensions to plot
% 
% Outputs:
%   F2d               Figure SM latent of the two target trials.
%   F                 Histogram of the d' values across all sessions,
%                       comparing d' for SM vs IM
%   h2di              t-test results for the 2D space
%   p2di              p-value for t-test for the 2D space
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


% Load in the D structure and the data
exampleSess = {'20190719'}; % The example session that we plot the neural
                            % space trajectories.
saveFig = 0;                % Determine whether or not to save the data,
                            % default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
colorMode = 'specified';    % Trajectory color mode to use
plotPreOnset = false;       % Plot pre-onset data
xSpec = 'xorth';            % Plot the orthonormalize GPFA space
markerScale = [0 1.5];      % Scale of the start and end trajectory markers
markSize = [3 10];          % Size of the midpoint markers (ind trls and mean)
plotDim = [1 2];            % Which dimensions are being plotting

% Assign the optional inputs
assignopts(who,varargin);

% Determine the sessions with the correct data
load(fullfile(dataLoc,'publicationQualitySessions.mat'));
dir_list = db.get_task_datasets(D, {'tt_int','tt_rot'});
load(fullfile(dataLoc,'dPrime_workspace+neuralSpace_fig3.mat'))

% Identify the example session
maskSessEx = find(ismember({dir_list(:,1).dataset},exampleSess));

% Create the 2D GPFA plots of the example session

% Load the decoder
[TD, ~, resultInt] = util.loadSessionData(dir_list(maskSessEx,1));
TD = TD.normalize([0 0 0]);
TD = TD(ismember([TD.targPos]',unique([TD.startPos]','rows'),'rows'));
[uniTarg,~,condition] = unique([TD.targPos]','rows');
[~, Prot,resultRot] = util.loadSessionData(dir_list(maskSessEx,2));
[~,~,TTrot] = orthogonalize(zeros(Prot.xDim,1),resultRot.estParams.C);
%[~,~,TT] = orthogonalize(zeros(P.xDim,1),result.estParams.C);

% Collect the optimization algorthim data
aRotTD = util.predictDecodeState_GPFA(TD,Prot,'spikeCountSrc','decodeSpikeCounts',...
    'useT',1,'TT',TTrot,'predictPos',1,'trunGPFA',1);
gp2D = [aRotTD.GPFA];
optMet = opt.get_data({gp2D.xorth},condition);
aRotTD = opt.applyProjection(aRotTD,Prot.ExpInfo.OptProj.M);

% Create the averages
gpa = [aRotTD.GPFA];
gpfaAA = gpa(condition==1).average();
gpfaBA = gpa(condition==2).average();

% Create the gpfa variables
gpfaA = gpa(condition==1);
gpfaB = gpa(condition==2);

% Create the optional inputs
optArg = {'plotAxLabels',false,...
    'plotPreOnset', plotPreOnset, ...
    'xSpec', xSpec, ...
    'plotDim', [1 2], ...
    'colorMode',colorMode,...
    'lineStr', '-'};

% Define the color maps and set up the figure
C = util.defineTaskColormap('bc_int');

F2d = figure('Position',[129 354.5000 1046 420]);

% Determine the indices for different color schemes
aIdx = find(C.endTC == find(ismember(round(C.targPos),round(uniTarg(1,:)),'rows')));
bIdx = find(C.endTC == find(ismember(round(C.targPos),round(uniTarg(2,:)),'rows')));

% Apply the projection to the midpoint data
optMet.x_midProj = Prot.ExpInfo.OptProj.M'*optMet.x_mid;
optMet.mu_ABproj = Prot.ExpInfo.OptProj.M'*optMet.mu_AB;
optMet.mu_BAproj = Prot.ExpInfo.OptProj.M'*optMet.mu_BA;

% Calculate the projection sigma
xMat_ABProj = Prot.ExpInfo.OptProj.M'*optMet.xMat_AB;
xMat_BAProj = Prot.ExpInfo.OptProj.M'*optMet.xMat_BA;
sig_AB = cov(xMat_ABProj');
sig_BA = cov(xMat_BAProj');

% Iterate through each subplot to plot the desired information
for n = 1:3
    subplot(1,3,n)
    hold on

    % Plot the individual trajectories for the 1st and 2nd panel
    if ismember(n,[1 2])
        gpfaA.trajectory('col',C.col_light(aIdx,:),'lineWidthavg',0.25,...
            'lineWidth',.25,'markerSize',3,'markerScale', markerScale(1),optArg{:})
        gpfaB.trajectory('col',C.col_light(bIdx,:),'lineWidthavg',0.25,...
            'lineWidth',.25,'markerSize',3,'markerScale', markerScale(1),optArg{:})
    end

    % Plot the average trajectories with the desire color scheme
    gpfaAA.trajectory('col',C.col_dark(aIdx,:),'lineWidthavg',2,...
        'markerSize',3,'markerScale', markerScale(2),optArg{:})
    gpfaBA.trajectory('col',C.col_dark(bIdx,:),'lineWidthavg',2,...
        'markerSize',3,'markerScale', markerScale(2),optArg{:})

    % Plot the individual trial midpoints
    if ismember(n,[2 3])
        plot(optMet.x_midProj(plotDim(1),condition==1),...
            optMet.x_midProj(plotDim(2),condition==1),...
            'd','MarkerSize',markSize(1),'MarkerFaceColor',C.col_dark(aIdx,:),...
            'MarkerEdgeColor',C.col_dark(aIdx,:))
        plot(optMet.x_midProj(plotDim(1),condition==2),...
            optMet.x_midProj(plotDim(2),condition==2),...
            'd','MarkerSize',markSize(1),'MarkerFaceColor',C.col_dark(bIdx,:),...
            'MarkerEdgeColor',C.col_dark(bIdx,:))
    end

    % Plot the midpoint mean, covariance, and d' line between the points
    if n == 3
        plt.plotCovEllipse(optMet.mu_ABproj,sig_AB,C.col_dark(aIdx,:),...
            'Marker','d','MarkerSize',markSize(2),'LineWidth',1.5,...
            'dim',plotDim,'LineStyle','-');
        plt.plotCovEllipse(optMet.mu_BAproj,sig_BA,C.col_dark(bIdx,:),...
            'Marker','d','MarkerSize',markSize(2),'LineWidth',1.5,...
            'dim',plotDim,'LineStyle','-');
        plot([optMet.mu_ABproj(plotDim(1)) optMet.mu_BAproj(plotDim(1))],...
            [optMet.mu_ABproj(plotDim(2)) optMet.mu_BAproj(plotDim(2))],...
            ':','Color','k','LineWidth',3)
    end

    % Set the axis view and label information
    axis square
    grid on
    xlabel('Neural dimension 2')
    ylabel('Neural dimension 1')
    zlabel('Neural dimension 3')
    view(90,90) % Rotates the space to align with the view in the paper.
end
plt.matchAxis(F2d);
plt.plotTitle(sprintf('%s%s GPFA Trajectories 2D',dir_list(1,1).subject,dir_list(1,1).dataset));
F2d.Name = sprintf('%s%s_gpfaTrajectories_2D',dir_list(1,1).subject,dir_list(1,1).dataset);

% Create a histogram of the different d' conditions
colIM = [180 178 180]./255;
colSM = [219 240 245]./255;

F = figure;
hold on
histogram(dP_WS(:,1),'BinWidth',0.25,'FaceColor',colIM)
histogram(dP_WS(:,2),'BinWidth',0.25,'FaceColor',colSM)
plot(mean(dP_WS(:,1)),0,'kv','MarkerFaceColor',colIM)
plot(mean(dP_WS(:,2)),0,'kv','MarkerFaceColor',colSM)
xlabel("d'"),ylabel("# Experiments")
title("Intuitive Mapping control trials - Workspace",'Interpreter','none')
ylim([0 30])
axis square
legend('IM trials-decoder_Int','IM trials-decoder_SM',...
    'Mean IM trials-decoder-Int','Mean IM trials-decoder-SM', ...
    'interpreter', 'none')
F.Name = 'workspace_dPrime_histograms';

% Paired t-test for the 2D space
[h2di,p2di] = ttest(dP_WS(:,1),dP_WS(:,2));

% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF([F F2d],savePathBase)
end
end
