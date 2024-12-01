% flow.create_session_figs  Create all session figures for flow
% comparison analysis.
%
% Usage:
%   [h] = flow.create_session_figs(FlowResults)
%
% This function takes the single-session summary structure FlowResults and
% creates all of the figures summarizing performance for that session.
%
% Optional Inputs:
%   arrow_size    Size of the arrow heads for the flow field
%   fig_save_loc  Location to save the flow field figures. 
%
% Output:
%   h             figure handles
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function [h] = create_session_figs(FlowResults, varargin)

% Specify location where figures are saved
arrow_size = 7;
fig_save_loc = [];
assignopts (who, varargin);

C_int = util.defineTaskColormap('bc_int');
C_rot = util.defineTaskColormap('bc_rot');
    
% Plot individual flow fields
[h_int] = flow.plot_flow_field(FlowResults.FF_int,'C',C_int,...
    'arrow_size',arrow_size);
[h_pred] = flow.plot_flow_field(FlowResults.FF_pred,'C',C_int,...
    'arrow_size',arrow_size);
[h_rot] = flow.plot_flow_field(FlowResults.FF_rot,'C',C_rot,...
    'arrow_size',arrow_size);

% Plot flow field comparision
[h_int_rot] = flow.plot_flow_comparison(FlowResults.A_int_rot,...
    'subject',FlowResults.subject,'dataset',FlowResults.dataset,...
    'grid',FlowResults.FF_int.grid.grid);
[h_pred_rot] = flow.plot_flow_comparison(FlowResults.A_pred_rot,...
    'subject',FlowResults.subject,'dataset',FlowResults.dataset,...
    'grid',FlowResults.FF_int.grid.grid);

% Plot comparison error
[h_session] = flow.plot_session_stats(FlowResults.A_int_rot, ...
    FlowResults.A_pred_rot,FlowResults.stats,...
    'subject',FlowResults.subject,'dataset',FlowResults.dataset);

% Save figures -- all figures for each session will be saved in an
% individual directory
if ~isempty(fig_save_loc)
    save_dir = fullfile(fig_save_loc, ...
        [FlowResults.subject, FlowResults.dataset]);
    mkdir(save_dir)
    h = [h_int, h_pred, h_rot, h_int_rot, h_pred_rot, h_session];
    saveFigurePDF(h, save_dir)
end