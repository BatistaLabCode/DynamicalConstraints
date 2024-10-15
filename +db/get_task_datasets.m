% db.get_task_datasets  Get set of datasets that have all the listed 
% desired tasks
%
% Usage:
%   db.get_task_datasets(D,task_list)
%
% Inputs:
%   D             Structure with all the valid experimental data. Save as
%                   filename <<publicationQualitySessions.mat>>
%   task_list     Cell struct listing the task_names to group by session.
% 
% Outputs:
%   dir_struct    Structure with the location information of all files.
%   task_datasets List of all valid experiment sessions.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function [dir_struct,task_datasets] = get_task_datasets(D, task_list)

% Get all subjects, datasets, tasks, and conditions
n = length(D);
all_datasets = {D.dataset};
all_tasks = {D.task};
all_conditions = {D.condition};

% Get unique datasets
uni_datasets = unique(all_datasets);

% Create masks for each task
n_tasks = length(task_list);
valid_task_mask = false(n, n_tasks);
for i = 1:n_tasks
    valid_task_mask(:, i) = ismember(all_tasks, task_list(i));
end

% To get all datasets, define a matrix of indices that can be applied to
% the entire dataset list.
task_datasets = nan(n, n_tasks);
idx = 0;
for ds = uni_datasets
    % Get mask for current dataset
    ds_mask = ismember(all_datasets, ds);
    
    % Find all possible conditions for the current dataset
    valid_mask = ds_mask & ismember(all_tasks, task_list);
    valid_cond = unique(all_conditions(valid_mask));
    n_cond = length(valid_cond);
    
    % Iterate over conditions and get tasks for each
    for i = 1:n_cond
        % Increment counter
        idx = idx + 1;
        for j = 1:n_tasks
            mask = ds_mask & ismember(all_tasks, task_list(j)) & ...
                ismember(all_conditions, valid_cond(i));
            task_idx = find(mask);
            if ~isempty(task_idx)
                task_datasets(idx, j) = task_idx;
            end
        end
    end
end

% Remove any datasets where all tasks were not present
nan_mask = sum(~isnan(task_datasets), 2) == n_tasks;
task_datasets = task_datasets(nan_mask, :);
n_ds = size(task_datasets, 1);

% Get structure of task datasets
D_task = repmat(D(1), n_ds, n_tasks);
for i = 1:n_ds
    for j = 1:n_tasks
        D_task(i, j) = D(task_datasets(i, j));
    end
end

% Get list of directories
dir_struct = db.get_dataset_dirs(D_task,'showWarning',0);