% flow.batch_create_session_figs  Create figs for all session
% experiments.
%
% Usage:
%   flow.batch_create_session_figs()
%
% Optional Inputs:
%   data_save_loc  Location of the saved flow field data.
%   fig_save_loc   Location to save the flow field figures.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function batch_create_session_figs(varargin)

% Define locations to save results structure
data_save_loc = [];
fig_save_loc = [];

assignopts (who, varargin);

% Load all data
if ~isempty(data_save_loc)
    FR = flow.load_session_results('data_save_loc',data_save_loc);
else
    FR = flow.load_session_results();
end

% Iterate over sessions and create plots
n_ds = length(FR);
invalid_mask = false(1, n_ds);
for i = 1:n_ds
    try
        fprintf('Creating figures for session %d of %d.\n', i, n_ds)
        if ~isempty(fig_save_loc)
            h = flow.create_session_figs(FR(i),'fig_save_loc',fig_save_loc);
        else
            h = flow.create_session_figs(FR(i));
        end
        close all
    catch
        invalid_mask(i) = true;
    end
end

% List sessions that generated errors
n_invalid = sum(invalid_mask);
if n_invalid > 0
    % Get list of invalid datasets
    all_subject = {FR.subject};
    all_dataset = {FR.dataset};
    invalid_subject = all_subject(invalid_mask);
    invalid_dataset = all_dataset(invalid_mask);
    
    % Iterate over invalid datasets and list
    fprintf('Could not create figures for the following session(s):\n')
    for i = 1:n_invalid
        fprintf('%s %s\n', invalid_subject{i}, invalid_dataset{i})
    end
end