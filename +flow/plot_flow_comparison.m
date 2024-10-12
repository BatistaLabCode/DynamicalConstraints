% plot_flow_comparison   Plot flow field comparison for two conditions
%
% Usage:
%   [F] = flow.plot_flow_comparison(A)
%
% Inputs:
%   A   Output of 'flow.compare_flow_fields' function
%
% Outputs:
%   F   Figure handle
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [F] = plot_flow_comparison(A)

subject = A.F_base.subject;
dataset = A.F_base.dataset;
cond_str = [A.F_base.condition '_' A.F_comp.condition];

% Set up figure
nRow = 3;
nCol = 2;
axSz = 300;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(nRow,nCol,axSz,axSz,axSp);

F = figure('Position',[100 100 fW fH]);
fig_str = sprintf('%s%s_FlowAnalysis_FlowComp_%s', subject, dataset, cond_str);
fName = sprintf(fig_str);
set(F,'Name',fName)

edges = 0:15:180;
plot_error_grid(A.ang.grid, A.F_base.grid.grid, 1, nRow, nCol, Ax, ...
    'Angular error', edges)

edges = 0:10:100;
plot_error_grid(A.mag.grid, A.F_base.grid.grid, 2, nRow, nCol, Ax, ...
    'Magnitude error', edges)

edges = 0:150:2000;
plot_error_grid(A.mse.grid, A.F_base.grid.grid, 3, nRow, nCol, Ax, ...
    'Mean squared error', edges)

% Plot title and set figure name
title_str = sprintf('%s %s :: Flow analysis :: Flow field comparison', ...
    subject, dataset);
plt.plotTitle(title_str)
end


% Function to plot heatmap and error
function plot_error_grid(grid_data, x, row_num, n_row, n_col, Ax, plot_title, edges)

plot_no = (row_num - 1) * 2 + 1;

% Subplot 1 -- heatmap of angular difference
plt.subplotSimple(n_row, n_col, plot_no, 'Ax', Ax); hold on;

imagesc(x, x, grid_data')
plt.colormapNew('red');
ax_lim = [-200 200];
set(gca, 'CLim', [0 edges(end)], 'XLim', ax_lim, 'YLim', ax_lim, ...
    'TickDir', 'out')
xlabel('X')
ylabel('Y')
title(sprintf('Voxel-wise difference: %s', plot_title))

% Subplot 2 -- heatmap of angular difference
plt.subplotSimple(n_row, n_col, plot_no + 1, 'Ax', Ax); hold on;

% Calculate average angle between vectors
grid_data = reshape(grid_data, 1, []);
grid_data = grid_data(~isnan(grid_data));
histogram(grid_data, edges, 'FaceColor', ones(1, 3) * 0.5)
set(gca, 'XLim', [0 edges(end)], 'TickDir', 'out')
plot(ones(1, 2) * median(grid_data), get(gca, 'YLim'), 'k--', 'LineWidth', 2)

xlabel('Difference')
ylabel('Counts')
title(sprintf('Distribution: %s', plot_title))

end