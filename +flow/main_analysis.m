% flow.main_analysis Main analysis function for flow analysis
%
% Usage:
%   flow.main_analysis(base_dir, intuitive_dir, rotated_dir)
%
% This function serves as a testbed for the neural flow analysis used to
% characterize similarity in dynamics between conditions.
%
% Inputs:
%   base_dir        Path for the main data folder
%   intuitive_dir   Path for the movement intuitive (MovInt) data folder
%   rotated_dir     Path for the separation maximizing (SepMax) data folder
%
% Optional Inputs:
%   d_grid             The size of the each workspace voxel grid.
%   min_pts_per_voxel  Minimum number of voxels to include for calculations
%   data_save_loc      Location to save flow field session data.
%   saveData           Determines whether save a FlowResult structure in
%                           data_save_loc folder.
% 
% Outputs:
%   FlowResult          Structure that includes all the summary flow fields
%                           used in the main paper figures and analysis for
%                           a given session.
% FlowResult structure is:
%   FF_int:      Flow Field (FF) for the MovInt trials
%   FF_int_rev:  FF for the time reversed MovInt trials  
%   FF_pred:     FF for MovInt trials projected into the SepMax workspace
%   FF_rot:      FF for the SepMax trials
%   FF_pred_Int: FF for SepMax trials projected into the MovInt workspace
%   A_int_rot:   Comparison of FF_int and FF_rot
%   A_pred_rot:  Comparison of FF_pred and FF_rot
%   stat:        Stats comparing error distribution for A_int_rot and A_pred_rot
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [FlowResults] = main_analysis(base_dir, intuitive_dir, rotated_dir, varargin)

% Parse optional arguments.
d_grid = 20;           % Allow the grid to be pre-specified
min_pts_per_voxel = 2; % Minmum number of points contained in a valid voxel
centerPos = [0 0 0];   % Defined center of the workspace
saveRaw = 0;           % Determine whether to save the raw grid data to the 
                       % comparison.
% Define locations to save results structure
data_save_loc = [];
saveData = 1;

assignopts (who, varargin);

dataLoc = serverPath;
if isempty(data_save_loc)
    data_save_loc = fullfile(dataLoc,'flowAnalysis','mat');
end

% Define dataset locations.
D.base = base_dir;
D.two_target_intuitive = intuitive_dir;
D.two_target_rotated = rotated_dir;

% Load data
TD_int = TrajectoryData().load(D.two_target_intuitive);
TD_rot = TrajectoryData().load(D.two_target_rotated);

if ~ismember(mean(unique([TD_rot.startPos]')),centerPos)
    TD_int = TD_int.normalize(mean(unique([TD_rot.startPos]','rows')));
    TD_rot = TD_rot.normalize(mean(unique([TD_rot.startPos]','rows')));
else
    TD_int = TD_int.normalize(centerPos);
    TD_rot = TD_rot.normalize(centerPos);
end

% Get subject, dataset, and decoder info and load decoder
if ismember(TD_int(1).subject,{'Monkey E'})
    subject = 'monkeyE';
    TD_int(1).subject = 'Monkey E';
    TD_rot(1).subject = 'Monkey E';
elseif ismember(TD_int(1).subject,{'Monkey D'})
    subject = 'monkeyD';
    TD_int(1).subject = 'Monkey D';
    TD_rot(1).subject = 'Monkey D';
elseif ismember(TD_int(1).subject,{'Monkey Q'})
    subject = 'monkeyQ';
    TD_int(1).subject = 'Monkey Q';
    TD_rot(1).subject = 'Monkey Q';
else
    subject = TD_int(1).subject;
end
dataset = datestr(TD_int(1).date, 'yyyymmdd');
D.rotated_mapping = [TD_rot(1).decoderName '.mat'];
D.intMov_mapping = [TD_int(1).decoderName '.mat'];

% Run analysis for intuitive mapping and estimate flow field
cond_str = 'Intuitive Movement';
[FF_int] = flow.calc_flow_field(TD_int, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Run analysis for the reverse direction of the intuitive mapping and
% estimate flow field
TD_int_rev = reverseDirection(TD_int);
cond_str = 'Intuitive Movement Reverse';
[FF_int_rev] = flow.calc_flow_field(TD_int_rev, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Load the rotated mapping
rot_map_path = fullfile(D.base, D.rotated_mapping);
M_rot = load(rot_map_path);

% Apply rotated mapping to intuitive trajectories
TD_rot_pred = util.predictDecodeState_GPFA(TD_int, M_rot.bci_params, ...
    'spikeCountSrc', 'decodeSpikeCounts');
cond_str = 'SepMax_Predicted';
[FF_pred] = flow.calc_flow_field(TD_rot_pred, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Run analysis for rotated two-target trials
cond_str = 'Separation Max';
[FF_rot] = flow.calc_flow_field(TD_rot, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Load the intuitive mapping
int_map_path = fullfile(D.base, D.intMov_mapping);
M = load(int_map_path);

% Apply intuitive mapping to rot trajectories
TD_int_pred = util.predictDecodeState_GPFA(TD_rot, M.bci_params, ...
    'spikeCountSrc', 'decodeSpikeCounts');
cond_str = 'IntMov_Predicted';
[FF_pred_Int] = flow.calc_flow_field(TD_int_pred, d_grid, cond_str, ...
    'min_pts_per_voxel', min_pts_per_voxel);

% Compare flow fields
[A_int_rot] = flow.compare_flow_fields(FF_int, FF_rot,'saveRaw',saveRaw,...
    'subject', subject, 'dataset', dataset, 'cond_str', 'IntMovVsSepMax');
[A_pred_rot] = flow.compare_flow_fields(FF_pred, FF_rot,'saveRaw',saveRaw,...
    'subject', subject, 'dataset', dataset, 'cond_str', 'PredVsSepMax');

% Summarize session (statistical tests)
[stats] = flow.calc_session_stats(A_int_rot, A_pred_rot);

% Add data to results structure
FlowResults.subject = TD_int(1).subject;
FlowResults.dataset = dataset;
FlowResults.FF_int = FF_int;
FlowResults.FF_int_rev = FF_int_rev;
FlowResults.FF_pred = FF_pred;
FlowResults.FF_rot = FF_rot;
FlowResults.FF_pred_Int = FF_pred_Int;
FlowResults.A_int_rot = A_int_rot;
FlowResults.A_pred_rot = A_pred_rot;
FlowResults.stats = stats;

% Save results
idxFile = length(base_dir);
intPos = regexp(intuitive_dir(idxFile+1:end),'_');
intPos = idxFile+(intPos(1)+1:intPos(2)-1);
rotPos = regexp(rotated_dir(idxFile+1:end),'_');
rotPos = idxFile+(rotPos(1)+1:rotPos(2)-1);
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
    revTD(n).brainKin.pos = flipud(TD(n).brainKin.pos);
end
end