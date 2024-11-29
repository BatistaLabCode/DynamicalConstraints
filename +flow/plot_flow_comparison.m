% plot_flow_comparison   Plot flow field comparison for two conditions
%
% Usage:
%   [F] = flow.plot_flow_comparison(A)
%
% This function creates two plots of the error distribution for a single 
% flow field comparison.
%      1) A heatmap of the error distribution values per voxel
%      2) A histogram of the error distribution values
%
% Inputs:
%   A   Output of 'flow.compare_flow_fields' function
%
% Optional Inputs:
%   subject                 Subject ID
%   dataset                 Experiment ID
%   grid                    Grid that defines the flow field
%   ax_lim                  Plotting limits of the flow field
%   edges                   Bin edges for error distribution
%   metric_str              Metric name in the structure (default: mse)
%   metric_plot_title       Metric name for printing text 
%
% Outputs:
%   F   Figure handle
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [F] = plot_flow_comparison(A,varargin)

% Optionally inputs
subject = '';           % Subject ID
dataset = '';           % Experiment ID
grid = [];              % Grid that defines the flow field
ax_lim = [-200 200];    % Plotting limits of the flow field
edges = 0:50:1500;      % Bin edges for error distribution
metric_str = 'mse';     % Metric name in the structure
metric_plot_title = 'Mean squared error';   % Metric name for printing text 

assignopts (who, varargin);

cond_str = [A.condition];

% Filled the undefined varables if information is available in the data
if isfield(A,'F_base') 
    if isempty(grid)
        grid = A.F_base.grid.grid;
    end
    if isempty(subject)
        subject = A.F_base.subject;
    end
    if isempty(dataset)
        dataset = A.F_base.dataset;
    end
else
    if isempty(grid)
        grid = linspace(-190,190,size(A.(metric_str).grid,1)+1);
        warning('No grid information, assuming the grid covers -190 to 190.')
    end
end

% Set up figure
nRow = 1;
nCol = 2;
axSz = 300;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(nRow,nCol,1.25*axSz,axSz,axSp);

F = figure('Position',[100 100 fW fH]);
fig_str = sprintf('%s%s_FlowAnalysis_FlowComp_%s', subject, dataset, cond_str);
fName = sprintf(fig_str);
set(F,'Name',fName)

plot_error_grid(A.(metric_str).grid, grid, 1, nRow, nCol, Ax, ...
    metric_plot_title, edges,ax_lim)

% Plot title and set figure name
title_str = sprintf('%s %s :: Flow analysis :: Flow field comparison :: %s', ...
    subject, dataset,cond_str);
plt.plotTitle(title_str)
end

% Function to plot heatmap and error
function plot_error_grid(grid_data, x, row_num, n_row, n_col, Ax, plot_title, edges,ax_lim)

plot_no = (row_num - 1) * 2 + 1;

% Subplot 1 -- heatmap of angular difference
plt.subplotSimple(n_row, n_col, plot_no, 'Ax', Ax); hold on;

imagesc(x, x, grid_data')
plt.colormapNew('red');
axis square
h = colorbar;
h.Label.String = ['Difference: ' plot_title];
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
axis square

xlabel('Difference')
ylabel('Counts')
title(sprintf('Distribution: %s', plot_title))

end