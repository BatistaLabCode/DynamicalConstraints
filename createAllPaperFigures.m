% Run this script to create all main figures in the paper. The default is
% to create all figures without saving the data, but there is an option to
% save the figures and even subselect the desired plots.
%
% Note: that there are additional varables you can adjust for each figure
% function. Review the code of each figure to determine the option options
% that you can modify.
%
% NOTE for fig_6_fig_7: That for plots where they are plotting a subset of
% trial trajectories. You can adjust the number trials for plotted using the
% variable <trialsPerCondition>. Also when selecting the subset of trials
% there are options to chose. 1) Use the hardcoding values that will
% exactly plot the trials in the paper. Vaiabile <useExample>. 2) Fix the 
% random number generator for a consistenet permutation for trial 
% subselection. Variable <setRandSeed>. 3) Let the set generation be 
% unfixed, meaning that each time you run, you will likely recieve a 
% different set of trials.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

clc, clear, close all

% Add the paper code to the path
pathName = pwd;
addpath(genpath(pathName))

% VARABLES to modify
plotFig = []; % Figures you would like plot
saveFig = 0;  % Determine if you want to save the figures or not, default is no.

if isempty(plotFig)
    plotFig = 2:6;
end

[dataLoc, exSessDataLoc, saveFigLoc] = serverPath;
figTime = [];
tic
strTime = toc;

% Create the plots for figure 2
if ismember(2,plotFig)
    [F2] = fig_2_two_target(dataLoc, ...
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
    [F4_proj] = fig_4_plot_view_vs_projections(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime(3) = toc;
    tic,
    [F4_el,h4,p4] = fig_4_early_vs_late_two_target_trajectories_comparisons(dataLoc,...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime(4) = toc;
end

% Create the plots for figure 5
if ismember(5,plotFig)
    tic,
    [F5_ff, F5_stat] = fig_5_plot_flow_field_analysis(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime = [figTime; toc];
end

% Create the plots for figure 6 and 7
if sum(ismember([6 7],plotFig))
    tic,
    [F6_hist,F6_intAng,F7_tube,F7_intAngTube,F6_EL] = fig_6_fig_7(dataLoc, ...
        'saveFig',saveFig,'savePathBase',saveFigLoc);
    figTime = [figTime; toc];
end

% Total run time
totTime = sum(figTime) - strTime;
sprintf('Total Run Time: %0.3f s',totTime)