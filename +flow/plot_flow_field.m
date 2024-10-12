% plot_flow_field   Plot estimated flow field
%
% [H] = flow.plot_flow_field(FF)
%
% Inputs:
%   FF      Flow field results structure (from el.flow.calc_flow_field)
%
% Outputs:
%   H       Vector of figure handles
%


function [h] = plot_flow_field(FF,varargin)

% Get color info
C = util.defineTaskColormap('bc');
arrow_size = 7; % If empty then no arrowhead is plotted.
pltGrid = 0;
assignopts (who, varargin);

% Set up figure
nRow = 3;
nCol = 3;
axSz = 300;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(nRow,nCol,axSz,axSz,axSp);

h = figure('Position',[100 100 fW fH]);
subject = FF.subject;
dataset = FF.dataset;
cond_str = FF.condition;
fig_str = sprintf('%s%s_FlowAnalysis_FlowField_%s', subject, dataset, cond_str);
fName = sprintf(fig_str);
set(h,'Name',fName)

% Subplot 1 -- scatter of all points
plt.subplotSimple(nRow,nCol,1,'Ax',Ax); hold on;

% For each starting position, calculate the distribution of point-to-point
% distances
ax_lim = FF.kin.ax_lim;
%ax_lim = [-1, 1] * 150;
start_pos_uni = FF.kin.start_pos;
n_start_pos = size(start_pos_uni, 1);
col_info_all = cell(n_start_pos, 1);
for i = 1:n_start_pos
    % Get color info for condition
    [col_info] = util.getColorInfo(C,start_pos_uni(i,:),[]);
    col_info_all{i} = col_info;

    % Scatter plot of positions
    pos = FF.kin.all_pos{i};
    scatter(pos(:,1), pos(:,2), 15, ...
        'MarkerFaceColor', col_info{1}, ...
        'MarkerEdgeColor', col_info{2})
    if pltGrid
        set(gca, 'XTick', FF.grid.grid, 'YTick', FF.grid.grid, 'TickDir', 'out')
        grid on
    else
        set(gca, 'XLim', ax_lim, 'YLim', ax_lim, 'TickDir', 'out')
    end
end
title('Cursor positions')
xlabel('X')
ylabel('Y')

% Subplot 2 -- histograms of distance between points
plt.subplotSimple(nRow,nCol,2,'Ax',Ax); hold on;

d_peak = FF.kin.d_peak;
hist_cts = FF.kin.hist_counts;
delta_hist = FF.kin.delta_hist;
edges = FF.kin.dist_edges;
for i = 1:n_start_pos
    % Plot data
    cts = hist_cts{i};
    x = edges(1:end-1) + delta_hist/2;
    plot(x, cts, 'color', col_info_all{i}{2}, 'LineWidth', 2)

    % Find peak value and plot
    y = max(cts);
    plot(ones(1, 2) * d_peak(i), [0, y], '--', ...
        'color', col_info_all{i}{2}, 'LineWidth', 2)
end

% Format plot
set(gca, 'TickDir', 'out')
title('Distance between points')
xlabel('Distance (mm)')
ylabel('Probability')

% Loop over grid and plot
grid_cts = FF.grid.grid_cts;
x_grid = FF.grid.grid;
d_grid = FF.grid.d_grid;
for i = 1:n_start_pos
    % Define colormap
    c = col_info_all{i}{2};
    c_hsv = rgb2hsv(c);
    cmap_hsv = repmat(c_hsv, 64, 1);
    s = linspace(c_hsv(2), 0, 64);
    cmap_hsv(:,2) = s;
    cmap = hsv2rgb(flipud(cmap_hsv));
    cmap(1,:) = ones(1, 3);

    % Setup subplot and plot map
    plt.subplotSimple(nRow,nCol,3+i,'Ax',Ax); hold on;
    C = squeeze(grid_cts(i,:,:))';  % Need to flip b/c the 'image' function flips x and y
    x = x_grid(1:end-1) + d_grid/2;
    y = x;
    imagesc(x,y,C)
    colormap(gca, cmap)
    set(gca,'CLim', [0 max(max(C))], 'TickDir', 'out')

    xlabel('X')
    ylabel('Y')
    title_str = sprintf('Occupancy histogram: Start position %d', i);
    title(title_str)
end

% Plot occupancy grid for combined data
plt.subplotSimple(nRow,nCol,6,'Ax',Ax); hold on;
C = squeeze(sum(grid_cts, 1))';  % Need to flip b/c the 'image' function flips x and y
x = x_grid(1:end-1) + d_grid/2;
y = x;
imagesc(x,y,C)
cmap = linspace(1,0,64);
cmap = [cmap', cmap', cmap'];
colormap(gca, cmap)
set(gca,'CLim', [0 max(max(C))], 'TickDir', 'out')

xlabel('X')
ylabel('Y')
title_str = sprintf('Occupancy histogram: Start position %d', i);
title(title_str)


% Plot map of counts per voxel (combined)
plt.subplotSimple(nRow,nCol,3,'Ax',Ax); hold on;

% Get hist data
grid_hist = FF.grid.grid_hist;
grid_hist_edges = FF.grid.grid_hist_edges;

plot(grid_hist_edges(1:end-1), grid_hist(1,:), ...
    'color', col_info_all{1}{2});
plot(grid_hist_edges(1:end-1), grid_hist(2,:), ...
    'color', col_info_all{2}{2});
plot(grid_hist_edges(1:end-1), sum(grid_hist, 1), 'k', ...
    'LineWidth', 2)
set(gca, 'YLim', [0 75], 'XLim', [1 grid_hist_edges(end-1)], ...
    'TickDir', 'out')
xlabel('Number of points per voxel')
ylabel('Counts (# voxels)')
title('Number of data points per voxel')

% Plot overall flow field
plt.subplotSimple(nRow,nCol,9,'Ax',Ax); hold on;

% Get data from structure
num_pts = FF.grid.num_pts;
X = FF.grid.X;
Y = FF.grid.Y;
VX_combined = FF.grid.VX;
VY_combined = FF.grid.VY;

n = reshape(num_pts, 1, []);
c = get_scaled_color(log(n), [0, 0, 0]);
x = reshape(X, 1, []);
y = reshape(Y, 1, []);
dx = reshape(VX_combined, 1, []);
dy = reshape(VY_combined, 1, []);
max_len = max(vecnorm([dx; dy]));
scale_factor = 2*d_grid/max_len;

plt.quiver_new(x, y, dx, dy, ...
    'scale_factor', scale_factor, ...
    'color', c, ...
    'line_width', 1.1,...
    'arrow_size',arrow_size)
colorbar
colormap(gca, flipud(unique(c,'rows')))
if pltGrid
    set(gca, 'XTick', FF.grid.grid, 'YTick', FF.grid.grid, 'TickDir', 'out', 'Box', 'on')
    grid on
else
    set(gca, 'XLim', ax_lim, 'YLim', ax_lim, 'TickDir', 'out', 'Box', 'on')
end
xlabel('X')
ylabel('Y')
title('Flow field - Combined')

% Plot flow field for each start position
num_pts_cond = FF.grid.num_pts_start_pos;
VX_start_pos = FF.grid.VX_start_pos;
VY_start_pos = FF.grid.VY_start_pos;
for i = 1:n_start_pos
    plt.subplotSimple(nRow,nCol,6+i,'Ax',Ax); hold on;

    n = reshape(num_pts_cond(:, :, i), 1, []);
    c = get_scaled_color(log(n), col_info_all{i}{2});
    dx = reshape(VX_start_pos(:,:,i), 1, []);
    dy = reshape(VY_start_pos(:,:,i), 1, []);
    max_len = max(vecnorm([dx; dy]));
    scale_factor = 2*d_grid/max_len;

    plt.quiver_new(x, y, dx, dy, ...
        'scale_factor', scale_factor, ...
        'color', c, ...
        'line_width', 1.1,...
        'arrow_size',arrow_size)
    if pltGrid
        set(gca, 'XTick', FF.grid.grid, 'YTick', FF.grid.grid, 'TickDir', 'out', 'Box', 'on')
        grid on
    else
        set(gca, 'XLim', ax_lim, 'YLim', ax_lim, 'TickDir', 'out', 'Box', 'on')
    end
    plt.scaleBar(gca,20,'mm')
    colorbar
    colormap(gca, flipud(unique(c,'rows')))
    xlabel('X')
    ylabel('Y')
    title(sprintf('Flow field - Start position %d', i))
end

title_str = sprintf('Energy landscape :: Flow field analysis :: %s %s, %s mapping', ...
    subject, dataset, cond_str);
plt.plotTitle(title_str)

end


% Function to get scaled colors for observation counts
function c = get_scaled_color(n, base_color)

% Find max number of observations -- this will define the range of colors
% that needs to be produced
max_n = log(90);

% Define colormap.  This should transition from gray to the base color by
% varying the saturation. No need to include white option as the data
% points won't be plotted anyway.
if sum(base_color == [0, 0, 0]) ~= 3
    % A color was specified
    base_color_hsv = rgb2hsv(base_color);
    col_map_hsv = repmat(base_color_hsv, 64, 1);
    s = linspace(0, base_color_hsv(2), 64);
    col_map_hsv(:, 2) = s;
    col_map = hsv2rgb(col_map_hsv);
else
    % Black was specified
    v = linspace(1, 0, 64);
    col_map = [v', v', v'];
end

% Get color for each point
n_pts = length(n);
c = nan(n_pts, 3);
for i = 1:n_pts
    c(i, :) = plt.getColorMapData(n(i), [0, max_n], col_map);
end

end