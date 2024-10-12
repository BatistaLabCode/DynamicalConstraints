% flow.batch_analysis  Run flow analysis on all datasets
%
% Usage:
%   flow.batch_analysis
%
% This function iterates over all two-target intuitive and rotated datasets
% and runs the 'flow.main_analysis' function.
%
% Optional Inputs:
%   data_save_loc  Location to save flow field batch data.
%   D              Structure with all the valid experimental data. Save as
%                   filename <<publicationQualitySessions.mat>>
%   task_list      Tasks to calculate the flow field for.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function batch_analysis(varargin)

% Define locations to save results structure
data_save_loc = 'F:\Erinn\testSave\FlowComparison\mat';
D = [];
task_list = { ...
    'tt_int',... % Shorthand flag for Two target intuitive
    'tt_rot'};   % Shorthand flag for Two target rotated
assignopts (who, varargin);

% Load dataset info
if isempty(D)
    [File_D, loc_D] = uigetfile('*.mat','Load the file containing the summary of sessions');
    load([loc_D File_D])
end

% Collect the directory information for valid files.
dir_list = db.get_task_datasets(D, task_list);

% Determine the number of sessions with a volid location for the data
mask = cellfun(@(x) length(x)>0,[dir_list(:,1).trajectory]);
dir_list = dir_list(mask,:);

% Iterate over task datasets
n_ds = size(dir_list, 1);
invalid_list = false(1, n_ds);
invalid_ds_list = cell(1, n_ds);
for i = 1:n_ds
    fprintf('Analyzing dataset %d of %d:\n', i, n_ds)
    try
        % Run analysis
        flow.main_analysis( ...
            dir_list(i, 1).base, ...
            dir_list(i, 1).trajectory{1}, ...
            dir_list(i, 2).trajectory{1}, ...
            'data_save_loc',data_save_loc);
        fprintf('done.\n')
    catch
        fprintf('Error encountered ... skipping dataset.\n\n')
        invalid_list(i) = true;
        invalid_ds_list{i} = dir_list(i, 1).base;
    end
end

% Display list of invalid datasets
n_invalid = sum(invalid_list);
if n_invalid > 0
    fprintf('\n\nThe following data directories were skipped:\n')
    invalid_ds_list = invalid_ds_list(invalid_list);
    for i = 1:n_invalid
        fprintf('%s\n', invalid_ds_list{i})
    end
end