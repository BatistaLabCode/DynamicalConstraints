% flow.plot_session_stats  Plot error distributions for flow analysis.
%
% Usage:
%   [h] = flow.plot_session_stats(A_int_rot, A_pred_rot, stats)
%
% Inputs:
%   A_int_rot   Flow field comparision between the intuitive and rotated
%   A_pred_rot  Flow field comparision between the predicted and rotated
%   stats       Structure containing results of statistical tests
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [h] = plot_session_stats(A_int_rot, A_pred_rot, stats)

% Setup figure
n_row = 1;
n_col = 3;
axSz = 300;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(n_row,n_col,axSz,axSz,axSp);

h = figure('Position',[10 10 fW fH]);
subject = A_int_rot.F_base.subject;
dataset = A_int_rot.F_base.dataset;
fig_str = sprintf('%s%s_FlowAnalysis_ErrorDist', subject, dataset);
fName = sprintf(fig_str);
set(h,'Name',fName)

% Plot distributions - angular error
plt.subplotSimple(n_row, n_col, 1, 'Ax', Ax); hold on;
edges = 0:15:180;
plot_metric_dist( ...
    A_int_rot.ang.valid, ...
    A_pred_rot.ang.valid, ...
    stats.pvals.ang.name, ...
    edges, ...
    stats.pvals.ang.p)

% Plot distributions - magnitude error
plt.subplotSimple(n_row, n_col, 2, 'Ax', Ax); hold on;
edges = 0:5:100;
plot_metric_dist( ...
    A_int_rot.mag.valid, ...
    A_pred_rot.mag.valid, ...
    stats.pvals.mag.name, ...
    edges, ...
    stats.pvals.mag.p)

% Plot distributions - mean squared error
plt.subplotSimple(n_row, n_col, 3, 'Ax', Ax); hold on;
edges = 0:50:1000;
plot_metric_dist( ...
    A_int_rot.mse.valid, ...
    A_pred_rot.mse.valid, ...
    stats.pvals.mse.name, ...
    edges, ...
    stats.pvals.mse.p)

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

% Get bin counts
hist_cts_int_rot = histcounts(x_int_rot, edges);
hist_cts_pred_rot = histcounts(x_pred_rot, edges);

% Plot histogram
col = {ones(1, 3) * 0.5, [66, 135, 245]/255};
ax_h_int_rot = histogram('BinEdges', edges, 'BinCounts', hist_cts_int_rot, ...
    'DisplayStyle', 'stairs', ...
    'EdgeColor', col{1}, ...
    'LineWidth', 2);
ax_h_pred_rot = histogram('BinEdges', edges, 'BinCounts', hist_cts_pred_rot, ...
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