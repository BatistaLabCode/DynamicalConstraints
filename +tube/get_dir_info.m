% el.tube.get_dir_info  Get data and figure saving paths for constrained
% path (tube) analysis
%

function [dir_info] = get_dir_info(subject, dataset,varargin)

% Define base save path (might update this to use the 'serverPath' function
% to allow data to be saved on the CMU file system)
save_path = '/Volumes/Samsung_T5/Analysis/EnergyLandscape/ConstrainedPathAnalysis';

optArg = assignopts(who,varargin);

% Define data path (subject-independent)
save_data_path = fullfile(save_path, 'mat');
int_targ_data_path = fullfile(save_data_path, 'int_targ_data');
suc_data_path = fullfile(save_data_path, 'success_rate_data');
[~, ~] = mkdir(int_targ_data_path);
[~, ~] = mkdir(suc_data_path);

% Define save paths and create directories (if necessary)
if nargin == 2
    fig_path = fullfile(save_path, 'fig', [subject, dataset]);
    [~, ~] = mkdir(fig_path);
else
    fig_path = [];
end

% Define paths to data
dir_info.save_path = save_path;
dir_info.fig_path = fig_path;
dir_info.int_targ_data_path = int_targ_data_path;
dir_info.suc_data_path = suc_data_path;