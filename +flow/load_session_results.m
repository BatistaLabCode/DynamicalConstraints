% flow.load_session_results  Load all session results data from disk.
%
% Usage:
%   [FR] = flow.load_session_results()
%
% This function returns a structure array of 'FlowResults' structures
% containing the results of the flow analysis run on each dataset.  Each
% element of FR is a structure of the type returned by the
% flow.main_analysis function.
% 
% Optional Inputs:
%   data_save_loc Location of the saved flow field data.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function [FR] = load_session_results(varargin)
data_save_loc = [];
assignopts (who, varargin);

if isempty(data_save_loc)
    dataLoc = serverPath;
    data_save_loc = fullfile(dataLoc,'flowAnalysis','mat');
end

% Define data location and get all *.mat files
valid_files = util.findDirContents(data_save_loc, '_FlowResults.mat');

% Load the first dataset to initialize the structure
FlowResults = load(fullfile(data_save_loc, valid_files{1}));
FlowResults = FlowResults.FlowResults;

% Iterate over datasets, load, and add to structure array
n_ds = length(valid_files);
FR = repmat(FlowResults, 1, n_ds);
for i = 1:n_ds
    FlowResults = load(fullfile(data_save_loc, valid_files{i}));
    FR(i) = FlowResults.FlowResults;
end