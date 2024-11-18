% flow.calc_session_stats  Calculate statistics for flow field
% comparison.
%
% Usage:
%   [S] = flow.calc_session_stats(A_int_rot, A_pred_rot)
%
% This function compares the error distributions of the provided datasets,
% which themselves are a summary of a flow field comparison.
%
% Compare the following:
% - mse distributions
% - It is possible to calculate other distribution by adding the
%       calculation to flow.compare_flow_fields.
%
% For each of these, compute the p-value
%
%
% Inputs:
%   A_int_rot        Flow field comparison for IntMov vs SepMax
%   A_pred_rot       Flow field comparison for predicted SepMax vs SepMax
%
% Outputs:
%   stats            Structure that includes the rank sum p-value between
%                       condition.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [stats] = calc_session_stats(A_int_rot, A_pred_rot)

% Keep track of the test used
stats.test = 'Wicoxon rank sum';

metric_fields = {'mse'};
metric_names = {'Mean squared error'};
n_metrics = length(metric_fields);

% Iterate over each metric and compute p-value
for i = 1:n_metrics
    field_str = metric_fields{i};
    p = ranksum(A_int_rot.(field_str).valid, A_pred_rot.(field_str).valid);
    stats.pvals.(field_str).name = metric_names{i};
    stats.pvals.(field_str).p = p;
end