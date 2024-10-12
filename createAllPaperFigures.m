% Run this code to create all main figures in the paper.The default is to
% run all figures without saving the data, but there is an option to
% save the figures and even subselect the desired plots.
%
% Note: that there are additional varables you can adjust for each figure
% function. Review the code of each figure to determine the option options
% that you can modify.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

clc, clear, close all

% VARABLES to modify
plotFig = []; % Figures you would like plot
saveFig = 0;  % Determine if you want to save the figures or not, default is no.

if isempty(plotFig)
    plotFig = 2:6;
end

[dataLoc, exSessDataLoc, saveFigLoc] = serverPath;
figTime = []
tic
strTime = toc;
% Create the plots for figure 2
if ismember(2,plotFig)
    [F2] = fig_2_twoTarget(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime = [figTime; toc];
end

% Create the plots for figure 3
if ismember(3,plotFig)
    tic,
    [F3_2d,F3,h3_2di,p3_2di] = fig_3(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime = [figTime; toc];
end

% Create the plots for figure 4
if ismember(4,plotFig)
    tic,
    [F4_proj] = fig_4_plot_View_vs_Projections(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime(3) = toc;
    tic,
    [F4_el,h4,p4] = fig_4_early_vs_late_TT_trajectories_comparisons(dataLoc,...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime(4) = toc;
end

% Create the plots for figure 5
if ismember(5,plotFig)
    tic,
    [F5, F5_ff, F5_ffcb, F5_stat] = fig_5_plot_FlowField_Analysis(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime = [figTime; toc];
end

% Create the plots for figure 6 and 7
if sum(ismember([6 7],plotFig))
    tic,
    [F6_hist,F6_intAng,F7_tube,F7_intAngTube,F6_EL] = fig_6_fig7(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime = [figTime; toc];
end

% Total run time
totTime = sum(figTime) - strTime;
sprintf('Total Run Time: %0.3f s',totTime)