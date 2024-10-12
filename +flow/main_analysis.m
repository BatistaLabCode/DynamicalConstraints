% Main analysis function for flow analysis
%
% This function serves as a testbed for the neural flow analysis used to
% characterize similarity in dynamics between conditions.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [FlowResults] = main_analysis(base_dir, intuitive_dir, rotated_dir, varargin)

% Parse optional arguments.
d_grid = 20;  % Allow the grid to be pre-specified
min_pts_per_voxel = 2;
% Define locations to save results structure
data_save_loc = 'F:\Erinn\testSave\FlowComparison\mat';
saveData = 1;

assignopts (who, varargin);


    
% Define dataset locations.  Eventually we will want to replace this with
% some sort of loop over valid datasets
D.base = base_dir;
D.two_target_intuitive = intuitive_dir;
D.two_target_rotated = rotated_dir;

% Load data
TD_int = TrajectoryData().load(D.two_target_intuitive);
TD_rot = TrajectoryData().load(D.two_target_rotated);

if ~ismember(mean(unique([TD_rot.startPos]')),[0 0 0])
TD_int = TD_int.normalize(mean(unique([TD_rot.startPos]','rows')));
TD_rot = TD_rot.normalize(mean(unique([TD_rot.startPos]','rows')));
end

% Get subject, dataset, and decoder info and load decoder
subject = TD_int(1).subject;
dataset = datestr(TD_int(1).date, 'yyyymmdd');
D.rotated_mapping = [TD_rot(1).decoderName '.mat'];

% Run analysis for intuitive mapping and estimate flow field
cond_str = 'Intuitive';
[FF_int] = flow.calc_flow_field(TD_int, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Run analysis for the reverse direction of the intuitive mapping and
% estimate flow field
TD_int_rev = reverseDirection(TD_int);
cond_str = 'IntuitiveReverse';
[FF_int_rev] = flow.calc_flow_field(TD_int_rev, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Load rotated mapping
rot_map_path = fullfile(D.base, D.rotated_mapping);
M = load(rot_map_path);

% Apply rotated mapping to intuitive trajectories
TD_rot_pred = util.predictDecodeState_GPFA(TD_int, M.bci_params, ...
    'spikeCountSrc', 'decodeSpikeCounts');
cond_str = 'Rotated_Predicted';
[FF_pred] = flow.calc_flow_field(TD_rot_pred, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Run analysis for rotated two-target trials
cond_str = 'Rotated';
[FF_rot] = flow.calc_flow_field(TD_rot, d_grid, cond_str);

% Compare flow fields
[A_int_rot] = flow.compare_flow_fields(FF_int, FF_rot, ...
    'subject', subject, 'dataset', dataset, 'cond_str', 'IntVsRot');
[A_pred_rot] = flow.compare_flow_fields(FF_pred, FF_rot, ...
    'subject', subject, 'dataset', dataset, 'cond_str', 'PredVsRot');

% Summarize session (statistical tests)
[stats] = flow.calc_session_stats(A_int_rot, A_pred_rot);

% Add data to results structure
FlowResults.subject = subject;
FlowResults.dataset = dataset;
FlowResults.params.d_grid = d_grid;
FlowResults.params.min_pts_per_voxel = min_pts_per_voxel;
FlowResults.FF_int = FF_int;
FlowResults.FF_int_rev = FF_int_rev;
FlowResults.FF_pred = FF_pred;
FlowResults.FF_rot = FF_rot;
FlowResults.A_int_rot = A_int_rot;
FlowResults.A_pred_rot = A_pred_rot;
FlowResults.stats = stats;

% Save results
intPos = regexp(intuitive_dir,'_');
intPos = [intPos(1)+1:intPos(2)-1];
rotPos = regexp(rotated_dir,'_');
rotPos = [rotPos(1)+1:rotPos(2)-1];
if saveData
    mkdir(data_save_loc)
    save_filename = sprintf('%s%s_int%s_rot%s_FlowResults.mat', subject,...
        dataset,intuitive_dir(rotPos),rotated_dir(rotPos));
    save_filename = fullfile(data_save_loc, save_filename);
    save(save_filename, 'FlowResults')
end
end
%% Functions
function [revTD] = reverseDirection(TD)
% Reverse the direction of the position data, need for calculating a
% portion of the flow field.
revTD = TD;

for n = 1: size(TD,1)
    revTD(n).pos = flipud(TD(n).pos);
end
end