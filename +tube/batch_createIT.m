% tube.batch_createIT  Batch analysis function creating the ID data
% structure
%
% Usage:
%   tube.batch_createIT
%
% This function creates a summary structure
%
% See: flow.main_analysis for description of the FlowResult structure.
%
% Optional Inputs:
%   data_save_loc  Location to save flow field batch data.
%   D              Structure with all the valid experimental data. Save as
%                   filename <<publicationQualitySessions.mat>>
%   task_list      Tasks to calculate the flow field for.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function invalid_list = batch_createIT(varargin)
save_path = [];
D = [];
task = 'rot_constr';

optArg = assignopts(who,varargin);

%if isempty(save_path)
%    save_path = uigetdir;
%end

if isempty(D)
    load('C:\Users\emg27\Dropbox\github\DynamicalConstraint\DynamicalConstraints_NatNeuro_2024_data\publicationQualitySessions.mat')
end

% Get all the tube sessions
maskTask = contains({D.task},task);
tD = D(maskTask);
[dataset,idx] = unique({tD.dataset});
subject = {tD(idx).subject};
ds_info = [subject' dataset'];
n_ds = size(ds_info, 1);

% Iterate over experiments and run analysis
invalid_mask = false(1, n_ds);
invalid_list = cell(1, n_ds);
for i = 1:n_ds
    try
        % Run analysis
        fprintf('Processing dataset %s %s ... ', ...
            ds_info{i, 1}, ds_info{i, 2})

        if ismember(ds_info{i,1},{'Quincy','monkeyQ'})
            if str2num(ds_info{i,2})>20210101
                centerPos = [-65 -330 0];
            else
                centerPos = [60,-270,0];
            end
        else
            centerPos = [0 0 0];
        end
        % Get IT object for dataset
        IT = tube.get_dataset_info(ds_info{i, 1}, ds_info{i, 2},'D',D);

        % Load data
        dir_list = db.get_dataset_dirs(D(ismember({D.dataset},dataset(i))));
        IT = tube.get_IT_data(IT,'centerPos',centerPos,...
            'dir_TD',dir_list(1));

        % Define data path (subject-independent)
        save_data_path = fullfile(save_path, 'mat');
        int_targ_data_path = fullfile(save_data_path, 'int_targ_data');
        [~, ~] = mkdir(int_targ_data_path);

        % Save results
        f_name_int_targ = [ds_info{i, 1}, ds_info{i, 2}, '_int_targ.mat'];
        save(fullfile(int_targ_data_path, f_name_int_targ), 'IT')

        fprintf('done.\n')
    catch ERR
        fprintf('error encountered.  Skipping ... \n')
        invalid_mask(i) = true;
        invalid_list{i} = sprintf('%s %s', ds_info{i, 1}, ds_info{i, 2});
    end
end

% List any datasets that were skipped due to errors
n_invalid = sum(invalid_mask);
if n_invalid > 0
    fprintf('\nThe following datasets were skipped due to errors:\n\n')
    invalid_list = invalid_list(invalid_mask);
    for i = 1:n_invalid
        fprintf('%s\n', invalid_list{i})
    end
end