function [FF] = calc_flow_field(TD, d_grid, cond_str, varargin)
% Calcuates the flow field of each unique start and end target pair. IE one
% flow field from A to B and another flow field for target B to A.
%
% Usage:
%   flow.calc_flow_field(TD, d_grid, cond_str)
%
% Inputs:
%   dataLoc    Paths for the main data folder
%
% Optional Inputs:
%   min_pts_per_voxel     Minimum number of voxels to include for
%                           calculations
%   ax_lim                Grid size
%   delta_hist            Distance between edges.
% 
% Outputs:
%   FF               Flow field structure for a given projection
%
% FF structure is:
%   subject
%   dataset 
%   condition:     Calculated flow field condition
%   kin:           FF for the SepMax trials
%     start_pos:      Start target position
%     all_pos:        Concatenated cursor positions, separated by start_pos
%     all_vel:        Concatenated cursor velocities, separated by start_pos
%     hist_counts:    Distribution of distance between cursor position at 
%                       time i and i+1. (Normalized)               
%     delta_hist:     Delta step for distance histogram
%     dist_edges:     Edges for distance histogram
%     d_peak:         Max distance for each target position (This is used
%                       to automatically determine grid size if not set)
%     ax_lim:         Axis limits
%   grid:          Flow field grid data
%     d_grid:            Number of voxel grids.
%     grid:              Grid edges
%     min_pts_per_vowel: Minimum number of voxels to include for calculations
%     num_pts:           Number of points per voxel for all conditions
%     num_pts_start_pos: Number of points per voxel per condition
%     X:                 Grid of x-direction position for voxels
%     Y:                 Grid of y-direction position for voxels
%     VX_start_pos:      Velocity flow field x-direction per condition
%     VY_start_pos:      Velocity flow field y-direction per condition
%     VX:                Flow field x-direction averaged across conditions
%     VY:                Flow field x-direction averaged across conditions
%     grid_cts:          Number of points in the per voxel per condition
%     grid_hist_edges:   Bins of the number of points per voxel
%     grid_hist:         Distribution of number of points per voxel
%
% Created by Erinn Grigsby and Alan Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Parameters
min_pts_per_voxel = 3;  % Define the minimum number of datapoints per voxel to use
ax_lim = [-200, 200];
delta_hist = 2.5;

assignopts(who, varargin);

% Preprocess -- two-target intuitive
center_pos = [0, 0, 0];
[TD,~] = util.preprocessGridTaskTrajData(TD,'centerPos',center_pos);

% Remove any non two-target trials
start_pos = [TD.startPos]';
end_pos = [TD.targPos]';
start_pos_uni = unique(start_pos, 'rows');
end_pos_mask = ismember(end_pos, start_pos_uni, 'rows');
TD = TD(end_pos_mask);
start_pos = start_pos(end_pos_mask, :);

% Define analysis parameters
edges = 0:delta_hist:ax_lim(2);  % Edges for distance histogram

% For each starting position, calculate the distribution of point-to-point
% distances
n_start_pos = size(start_pos_uni, 1);
hist_cts = cell(n_start_pos, 1);
all_pos = cell(n_start_pos, 1);
all_vel = cell(n_start_pos, 1);
for i = 1:n_start_pos
    % Get all cursor positions.  Loop over all trials and calculate
    % velocity from positions.
    start_pos_mask = ismember(start_pos, start_pos_uni(i,:), 'rows');
    tmpVRCurDat = [TD(start_pos_mask).brainKin];
    pos = {tmpVRCurDat.pos}';
    n_trials = length(pos);
    vel = cell(n_trials, 1);
    for j = 1:length(pos)
        pos_temp = pos{j}(:,1:2);
        vel_temp = diff(pos_temp,[],1);
        vel{j} = vel_temp;
        pos{j} = pos_temp(1:end-1,:);
    end
    pos = cell2mat(pos);
    vel = cell2mat(vel);
    all_pos{i} = pos;
    all_vel{i} = vel;
    
    % Calculate distance between points
    num_samp = size(pos, 1);
    d_all = cell(num_samp, 1);
    for j = 1:(size(pos, 1) - 1)
        pos_temp = pos(j+1:end,:);  % Get all points after current point
        d = pos_temp - pos(j,:);  % Vector from current point to other points
        d = vecnorm(d');
        d_all{j} = d;
    end 
    
    % Calcuate histogram
    d = cell2mat(d_all');
    cts = histcounts(d, edges);
    [hist_cts{i}] = cts/sum(cts);  % Normalize
end

d_peak = nan(n_start_pos, 1);
for i = 1:n_start_pos
    % Plot data
    cts = hist_cts{i};
    x = edges(1:end-1) + delta_hist/2;
    
    % Find peak value and plot
    [~, I] = max(cts);
    d_peak(i) = x(I);
end

% Calculate grid size.  Define this as half the peak distance
if isempty(d_grid)
    d_grid = mean(d_peak)/2;
end
x_grid = d_grid/2:d_grid:ax_lim(2);  % Center the grid at 0
x_grid = [-fliplr(x_grid) x_grid];

% Loop over grid and calculate the number of points per voxel
n_grid = length(x_grid) - 1;
grid_cts = nan(n_start_pos, n_grid, n_grid);
for i = 1:n_start_pos
    % Get all positions for current target
    pos = all_pos{i};
    
    % Loop over grid and find the number of points in each voxel
    for j = 1:n_grid
        x_lim = [x_grid(j), x_grid(j+1)];
        for k = 1:n_grid
            % Find points in the current voxel
            y_lim = [x_grid(k), x_grid(k+1)];
            x_mask = (pos(:, 1) >= x_lim(1)) & (pos(:, 1) < x_lim(2));
            y_mask = (pos(:, 2) >= y_lim(1)) & (pos(:, 2) < y_lim(2));
            mask = x_mask & y_mask;
            
            % Calcualte the number of points in the current voxel
            grid_cts(i,j,k) = sum(mask);
        end
    end
end

% Determine the distribution of voxels by point count
[~,grid_hist_edges] = histcounts(reshape(grid_cts,1,[]));
grid_hist = nan(n_start_pos, length(grid_hist_edges) - 1);
for i = 1:n_start_pos
    % Determine counts per bin
    temp_cts = reshape(grid_cts(i,:,:), 1, []);
    grid_hist(i,:) = histcounts(temp_cts, grid_hist_edges);
end

% Loop over grid and calculate flow field
VX_start_pos = nan(n_grid, n_grid, n_start_pos);
VY_start_pos = nan(n_grid, n_grid, n_start_pos);
VX_combined = nan(n_grid, n_grid);
VY_combined = nan(n_grid, n_grid);
X = nan(n_grid, n_grid);
Y = nan(n_grid, n_grid);
num_pts = nan(n_grid, n_grid);
num_pts_cond = nan(n_grid, n_grid, n_start_pos);
for i = 1:n_grid
    x_lim = [x_grid(i), x_grid(i+1)];
    for j = 1:n_grid
        y_lim = [x_grid(j), x_grid(j+1)];
        vel_combined = cell(n_start_pos, 1);
        for k = 1:n_start_pos
            % Find the valid points for the current voxel
            pos = all_pos{k};
            vel = all_vel{k};
            x_mask = (pos(:, 1) >= x_lim(1)) & (pos(:, 1) < x_lim(2));
            y_mask = (pos(:, 2) >= y_lim(1)) & (pos(:, 2) < y_lim(2));
            mask = x_mask & y_mask;
            
            % Get all valid velocity vectors
            if sum(mask) > 0
                vel_temp = vel(mask,:);
                vel_combined{k} = vel_temp;
                
                % Only include the point in the per-start condition average
                % if sufficient samples exist
                if sum(mask) >= min_pts_per_voxel
                    if sum(mask) > 1  % Handle case where there is only 1 element
                        vel_temp = mean(vel_temp, 1);
                    end
                else
                    vel_temp = [nan, nan];
                end
            else
                vel_temp = [nan, nan];
            end
            
            % Add velocity vectors to matrix
            VX_start_pos(i, j, k) = vel_temp(1);
            VY_start_pos(i, j, k) = vel_temp(2);
            num_pts_cond(i, j, k) = sum(mask);
        end
        
        % Combine velocity vectors across start positions
        vel_temp = cell2mat(vel_combined);
        num_pts(i, j) = size(vel_temp, 1);
        
        % If sufficient samples exist, average and add to combined flow
        % field.
        if size(vel_temp, 1) >= min_pts_per_voxel
            if size(vel_temp, 1) > 1
                vel_temp = mean(vel_temp, 1);
            end
            VX_combined(i, j) = vel_temp(1);
            VY_combined(i, j) = vel_temp(2);
        end
        X(i,j) = x_grid(i);
        Y(i,j) = x_grid(j);
    end
end

% Pack info up into results structure
FF.subject = TD(1).subject;
FF.dataset = datestr(TD(1).date, 'yyyymmdd');
FF.condition = cond_str;

% Kinematic data
FF.kin.start_pos = start_pos_uni;
FF.kin.all_pos = all_pos;
FF.kin.all_vel = all_vel;
FF.kin.hist_counts = hist_cts;
FF.kin.delta_hist = delta_hist;
FF.kin.dist_edges = edges;
FF.kin.d_peak = d_peak;
FF.kin.ax_lim = ax_lim;

% Flow field grid data
FF.grid.d_grid = d_grid;
FF.grid.grid = x_grid;
FF.grid.min_pts_per_voxel = min_pts_per_voxel;
FF.grid.num_pts = num_pts;
FF.grid.num_pts_start_pos = num_pts_cond;
FF.grid.X = X;
FF.grid.Y = Y;
FF.grid.VX = VX_combined;
FF.grid.VY = VY_combined;
FF.grid.VX_start_pos = VX_start_pos;
FF.grid.VY_start_pos = VY_start_pos;
FF.grid.grid_cts = grid_cts;
FF.grid.grid_hist_edges = grid_hist_edges;
FF.grid.grid_hist = grid_hist;

end

