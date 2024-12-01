% flow.plot_session_stats  Plot error distributions for flow analysis.
%
% Usage:
%   [h] = flow.plot_session_stats(A_int_rot, A_pred_rot, stats)
%
% This function plots error distributions of a single session for two flow
%   field comparisons.
%
% Plots the following:
% - mse distributions
% - Note: It is possible to calculate other distribution by adding the
%       calculation to flow.compare_flow_fields.
%
% Inputs:
%   A_int_rot   Flow field comparision between the intuitive and rotated
%   A_pred_rot  Flow field comparision between the predicted and rotated
%   stats       Structure containing results of statistical tests
%
% Optional Inputs:
%   subject                 Subject ID
%   dataset                 Experiment ID
%   edges                   Bin edges for error distribution
%   metric_str              Metric name in the structure (default: mse)
%
% Outputs:
%   F   Figure handle
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [h] = plot_session_stats(A_int_rot, A_pred_rot, stats,varargin)
subject = '';           % Subject ID
dataset = '';           % Experiment ID
edges = 0:50:1500;      % Bin edges for error distribution
metric_str = 'mse';     % Metric name in the structure 

assignopts (who, varargin);

% Filled the undefined varables if information is available in the data
if isfield(A_int_rot,'F_base') 
    if isempty(subject)
        subject = A_int_rot.F_base.subject;
    end
    if isempty(dataset)
        dataset = A_int_rot.F_base.dataset;
    end
end

% Setup figure
n_row = 1;
n_col = 1;
axSz = 300;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(n_row,n_col,axSz,axSz,axSp);

h = figure('Position',[10 10 fW fH]);
fName = sprintf('%s%s_FlowAnalysis_ErrorDist', subject, dataset);
set(h,'Name',fName)

% Plot distributions - mean squared error
plt.subplotSimple(n_row, n_col, 1, 'Ax', Ax); hold on;
plot_metric_dist( ...
    A_int_rot.(metric_str).valid, ...
    A_pred_rot.(metric_str).valid, ...
    stats.pvals.(metric_str).name, ...
    edges, ...
    stats.pvals.(metric_str).p)

% Plot figure title
title_str = sprintf('%s %s :: Flow analysis :: Error distributions', ...
    subject, dataset);
plt.plotTitle(title_str)
end


% Function to plot distributions of each metric
function plot_metric_dist(x_int_rot, x_pred_rot, name_str, edges, pval)

% Define font sizes
font_size.label = 14;
font_size.tick = 14;
font_size.title = 14;

% Plot histogram
col = {ones(1, 3) * 0.5, [66, 135, 245]/255};
ax_h_int_rot = histogram(x_int_rot, edges, ...
    'DisplayStyle', 'stairs', ...
    'EdgeColor', col{1}, ...
    'LineWidth', 2);
ax_h_pred_rot = histogram(x_pred_rot, edges, ...
    'DisplayStyle', 'stairs', ...
    'EdgeColor', col{2}, ...
    'LineWidth', 2);

% Plot median
plot(median(x_int_rot) * ones(1, 2), get(gca, 'YLim'),  ...
    '--', 'color', 'k', 'LineWidth', 2)
plot(median(x_pred_rot) * ones(1, 2), get(gca, 'YLim'), ...
    '--', 'color', col{2}, 'LineWidth', 2)

% Set plot title
legend([ax_h_int_rot, ax_h_pred_rot], ...
    {'Intuitive vs. rotated', 'Predicted vs rotated'}, 'Box', 'off')
set(gca, 'XLim', [edges(1), edges(end)], 'TickDir', 'out')
ylabel('Counts (voxels)', 'FontSize', font_size.label)
xlabel(name_str, 'FontSize', font_size.label)
title(sprintf('%s (p = %0.3e)', name_str, pval))

set(get(gca,'XAxis'), 'FontSize', font_size.tick)
set(get(gca,'YAxis'), 'FontSize', font_size.tick)

end