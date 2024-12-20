% defineTaskColormap       Get colormap for brain control analysis
%
% This function returns a colormap providing target position, target label,
% target code for both start and end target (TC), and two color maps. 
% col_dark is the highest saturation and col_light is its lighter counter 
% part.
%
% Usage:
%   util.defineTaskColormap(task)
% 
% Inputs:
%   task           Task Condition with distinct color schemes. Options
%                       include: bc, bc_int, bc_rot, bc_grid, bc_fine,
%                       and centerOut.
%
% Optional Inputs:
%   target_radius  The distance of the targets from the center point
% 
% Outputs:
%   C              Colormap structure for the given task
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function C = defineTaskColormap(task, varargin)

% Parse optional arguments
target_radius = 90;

assignopts(who, varargin);

if isempty(task)
    task = 'bc';
end

switch task
    case 'bc'
        % Color map for the standard two-target conditional grid task
        
        % Define all valid target positions and their corresponding labels
        targPos = [ ...
            1 0 0; ...
            [1 1 0] * sqrt(2)/2; ...
            0 1 0; ...
            [-1 1 0] * sqrt(2)/2; ...
            -1 0 0; ...
            [-1 -1 0] * sqrt(2)/2; ...
            0 -1 0; ...
            [1 -1 0] * sqrt(2)/2];
        
        targLabel = { ...
            'T1', ...
            'T2', ...
            'T3', ...
            'T4', ...
            'T5', ...
            'T6', ...
            'T7', ...
            'T8'};
        
        % Define start and end target codes, which correspond to the
        % 'targPos' matrix defined above
        startTC = [1, 5, 2, 6, 3, 7, 4, 8];
        endTC = [5, 1, 6, 2, 7, 3, 8, 4];
        
        % Define colors for start/end target pairs
        col_dark = [ ...
            11 100 80; ... % red
            197 100 70; ... % blue
            118 100 65; ... % green
            283 100 85; ... % purple
            40 100 80; ...   % yellow-orange
            220 100 80; ...  % dark blue
            85 100 80; ...   % yellow-green
            330 100 80];   % pink
        
        col_light = [ ...
            11 30 80; ... % red
            197 30 80; ... % blue
            118 30 80; ... % green
            283 30 85; ... % purple
            40 30 80; ...   % yellow-orange
            220 30 80; ...  % dark blue
            85 30 80; ...   % yellow-green
            330 30 80];   % pink
        
        normalize_flag = true;
        
    case 'bc_int'
        % Color map for the standard two-target conditional grid task, this
        % sets all the colors the samne for the two target task for the
        % intuitive mapping.
        
        % Define all valid target positions and their corresponding labels
        targPos = [ ...
            1 0 0; ...
            [1 1 0] * sqrt(2)/2; ...
            0 1 0; ...
            [-1 1 0] * sqrt(2)/2; ...
            -1 0 0; ...
            [-1 -1 0] * sqrt(2)/2; ...
            0 -1 0; ...
            [1 -1 0] * sqrt(2)/2];
        
        targLabel = { ...
            'T1', ...
            'T2', ...
            'T3', ...
            'T4', ...
            'T5', ...
            'T6', ...
            'T7', ...
            'T8'};
        
        % Define start and end target codes, which correspond to the
        % 'targPos' matrix defined above
        startTC = [1, 5, 2, 6, 3, 7, 4, 8];
        endTC = [5, 1, 6, 2, 7, 3, 8, 4];
        
        % Define colors for start/end target pairs
        col_dark = [ ...
            11 100 80; ...  % red
            197 100 70; ... % blue
            11 100 80; ...  % red
            197 100 70; ... % blue
            11 100 80; ...  % red
            197 100 70; ... % blue
            11 100 80; ...  % red
            197 100 70];    % blue
        
        col_light = [ ...
            11 30 80; ...  % red
            197 30 80; ... % blue
            11 30 80; ...  % red
            197 30 80; ... % blue
            11 30 80; ...  % red
            197 30 80; ... % blue
            11 30 80; ...  % red
            197 30 80];    % blue
        
        normalize_flag = true;
    case 'bc_rot'
        % Color map for the standard two-target conditional grid task, this
        % sets all the colors the samne for the two target task for the
        % rotated mapping.
        
        % Define all valid target positions and their corresponding labels
        targPos = [ ...
            1 0 0; ...
            [1 1 0] * sqrt(2)/2; ...
            0 1 0; ...
            [-1 1 0] * sqrt(2)/2; ...
            -1 0 0; ...
            [-1 -1 0] * sqrt(2)/2; ...
            0 -1 0; ...
            [1 -1 0] * sqrt(2)/2];
        
        targLabel = { ...
            'T1', ...
            'T2', ...
            'T3', ...
            'T4', ...
            'T5', ...
            'T6', ...
            'T7', ...
            'T8'};
        
        % Define start and end target codes, which correspond to the
        % 'targPos' matrix defined above
        startTC = [1, 5, 2, 6, 3, 7, 4, 8];
        endTC = [5, 1, 6, 2, 7, 3, 8, 4];
        
        % Define colors for start/end target pairs
        col_dark = [ ...
            358 87 44; ...   % dark red
            198 84 41; ...   % navy
            358 87 44; ...   % dark red
            198 84 41; ...   % navy
            358 87 44; ...   % dark red
            198 84 41; ...   % navy
            358 87 44; ...   % dark red
            198 84 41]; ...  % navy
        
        col_light = [ ...
            6 62 62; ...   % rust
            206 51.3 61.2; ...  % dark french blue
            6 62 62; ...   % rust
            206 51.3 61.2; ...  % dark french blue
            6 62 62; ...   % rust
            206 51.3 61.2; ...  % dark french blue
            6 62 62; ...   % rust
            206 51.3 61.2]; ... % dark french blue
        
        normalize_flag = true;
    case 'bc_grid'
        % Color map for the standard two-target conditional grid task, this
        % sets all the colors the same for the two target task for the
        % rotated mapping.
        
        % Define all valid target positions and their corresponding labels
        targPos = [ ...
            1 0 0; ...
            [1 1 0] * sqrt(2)/2; ...
            0 1 0; ...
            [-1 1 0] * sqrt(2)/2; ...
            -1 0 0; ...
            [-1 -1 0] * sqrt(2)/2; ...
            0 -1 0; ...
            [1 -1 0] * sqrt(2)/2];
        
        targLabel = { ...
            'T1', ...
            'T2', ...
            'T3', ...
            'T4', ...
            'T5', ...
            'T6', ...
            'T7', ...
            'T8'};
        
        % Define start and end target codes, which correspond to the
        % 'targPos' matrix defined above
        startTC = [1, 5, 1, 1, 5, 5];
        endTC = [5, 1, 3, 7, 3, 7];
        
        % Define colors for start/end target pairs
        col_dark = [ ...
            11 100 80; ...  % red
            197 100 70; ... % blue
            40 100 80; ...   % yellow-orange
            40 100 80; ...   % yellow-orange
            220 100 80; ...  % dark blue
            220 100 80];     % dark blue
        
        col_light = [ ...
            11 30 80; ...  % red
            197 30 80; ... % blue
            40 30 80; ...   % yellow-orange
            40 30 80; ...   % yellow-orange
            220 30 80; ...  % dark blue
            220 30 80];     % dark blue
        
        normalize_flag = true;
    case 'centerOut'
        % 8-target center-out
        
        % Define all valid target positions and their corresponding labels
        targPos = [ ...
            0 0 0; ...
            1 0 0; ...
            [1 1 0] * sqrt(2)/2; ...
            0 1 0; ...
            [-1 1 0] * sqrt(2)/2; ...
            -1 0 0; ...
            [-1 -1 0] * sqrt(2)/2; ...
            0 -1 0; ...
            [1 -1 0] * sqrt(2)/2];
        
        targLabel = { ...
            'C', ...
            'T1', ...
            'T2', ...
            'T3', ...
            'T4', ...
            'T5', ...
            'T6', ...
            'T7', ...
            'T8'};
        
        % Define start and end target codes, which correspond to the
        % 'targPos' matrix defined above
        startTC = [1, 1, 1, 1, 1, 1, 1, 1];
        endTC = [2, 3, 4, 5, 6, 7, 8, 9];
        
        % Define colors for start/end target pairs
        col_dark = interpColor(8, 1, 'hsv', 1);
        col_dark = col_dark([1, 5, 2, 6, 3, 7, 4, 8], :);  % alternate color
        col_dark = rgb2hsv(col_dark);
        col_light = col_dark;
        col_light(:,2) = col_light(:,2) * 0.5;
        
        normalize_flag = false;
        
    case 'bc_fine'
        % Color map for the standard two-target conditional grid task
        
        % Define all valid target positions and their corresponding labels
        targPos = [ ...
            1 0 0; ...
            [1 1 0] * sqrt(2)/2; ...
            0 1 0; ...
            [-1 1 0] * sqrt(2)/2; ...
            -1 0 0; ...
            [-1 -1 0] * sqrt(2)/2; ...
            0 -1 0; ...
            [1 -1 0] * sqrt(2)/2;...
            [cosd(22.5:45:360)' sind(22.5:45:360)' [0 0 0 0 0 0 0 0]']];
        
        targLabel = { ...
            'T1', ...
            'T2', ...
            'T3', ...
            'T4', ...
            'T5', ...
            'T6', ...
            'T7', ...
            'T8', ...
            'T9', ...
            'T10',...
            'T11', ...
            'T12', ...
            'T13', ...
            'T14', ...
            'T15', ...
            'T16'};
        
        % Define start and end target codes, which correspond to the
        % 'targPos' matrix defined above
        startTC = [1, 5, 2, 6, 3, 7, 4, 8, 9, 13, 10, 14, 11, 15, 12, 16];
        endTC = [5, 1, 6, 2, 7, 3, 8, 4, 13, 9, 14, 10, 15, 11, 16, 12];
        
        % Define colors for start/end target pairs
        col_dark = [ ...
            11 100 80; ... % red
            197 100 70; ... % blue
            118 100 65; ... % green
            283 100 85; ... % purple
            40 100 80; ...   % yellow-orange
            220 100 80; ...  % dark blue
            85 100 80; ...   % yellow-green
            330 100 80; ...   % pink
            22 100 80; ... % ??
            210 100 70; ... %? blue
            80 100 65; ... %? green
            240 100 85; ... % ?purple
            60 100 80; ...   % ?yellow-orange
            180 100 80; ...  %? dark blue
            100 100 80; ...   %? yellow-green
            350 100 80];   % ?pink
        
        col_light = [ ...
            11 30 80; ... % red
            197 30 80; ... % blue
            118 30 80; ... % green
            283 30 85; ... % purple
            40 30 80; ...   % yellow-orange
            220 30 80; ...  % dark blue
            85 30 80; ...   % yellow-green
            330 30 80; ...
            22 30 80; ... % ??
            210 30 70; ... %? blue
            80 30 65; ... %? green
            240 30 85; ... % ?purple
            60 30 80; ...   % ?yellow-orange
            180 30 80; ...  %? dark blue
            100 30 80; ...   %? yellow-green
            350 30 80];   % ?pink
        
        normalize_flag = true;
end

% Scale target position
targPos = targPos * target_radius;

% Specify undefined colors
undef_col_dark = 0.5 * ones(1,3);
undef_col_light = 0.75 * ones(1,3);

% Normalize color arrays
if normalize_flag
    nF = [360 100 100]; % Normalization factor
    col_light = col_light./nF;
    col_dark = col_dark./nF;
end

% Convert hsv colors to rgb
col_light = hsv2rgb(col_light);
col_dark = hsv2rgb(col_dark);

% Pack up data in structure
C.task = task;
C.targPos = targPos;
C.targLabel = targLabel;
C.startTC = startTC;
C.endTC = endTC;
C.col_light = col_light;
C.col_dark = col_dark;
C.undef_col_light = undef_col_light;
C.undef_col_dark = undef_col_dark;