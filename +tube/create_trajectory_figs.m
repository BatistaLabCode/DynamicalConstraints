% tube.create_trajectory_figs  Create trajectory figures for constrained
% path (tube) experiment.
%
% Usage:
%   tube.create_trajectory_figs(IT)
%
% This function iterates through the two stages of the constrained path
% task (unconstrained and constrained) and creates summary plots of 
% the behavioral trajectories in the tubes.
%
% Input:
%   IT                  IntTargExp object
%
% Optional Inputs:
%   plot_scale          Sets the plotting scale
%   center_pos          Center of the workspace
%   r_factor            Factor used for spatial averaging
%   plotStates          Determine which task states to include  
%   avgMode             Averaging method
%   plotSegment         Determine whether to plot only portions of the tube
%   C                   Define color map
%   trialsPerCondition  Allow for subselection of trials
%   setRandSeed         Fix the random number generator for a consistenet
%                           permutation for trial subselection. 
%
% Output:
%   fh          Figure handle
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [fh] = create_trajectory_figs(IT, varargin)

plot_scale = 200;
center_pos = [0, 0, 0];
r_factor = 0.15;
plotStates = {'Step 2'};
avgMode = 'spatial';
plotSegment = [];
C = [];
trialsPerCondition = [];     % Allow for subselection of trials
setRandSeed = [];

assignopts(who, varargin);

% Get color info for plotting
if isempty(C)
    C = util.defineTaskColormap('bc');
    C.targPos = C.targPos.*(norm(IT.startTargPos)/100)/.9;
end

% Create figures for intermediate target - no tube
saveNameBase = IT.unconstrainedBlockDir(1:end-18);
f_unconst = util.plotIntermediateTargetTaskTrajectories( ...
    IT.TDunconstrained, C, saveNameBase, ...
        'centerPos', center_pos, ...
        'plotTubes', false, ...
        'plotScale', plot_scale, ...
        'avgMode', avgMode, ...
        'r_factor', r_factor,...
        'plotStates',plotStates,...
        'plotSegment',plotSegment,...
        'trialsPerCondition',trialsPerCondition,...
        'setRandSeed',setRandSeed);

% Create figures for intermediate target - tube
nBlocks = length(IT.constrainedBlockDir);
f_const = nan(1,nBlocks);
for j = 1:nBlocks
    saveNameBase = IT.constrainedBlockDir{j}(1:end-18);
    f_const(j) = util.plotIntermediateTargetTaskTrajectories( ...
        IT.TDconstrained{j},C,saveNameBase, ...
        'centerPos', center_pos, ...
        'plotTubes', true, ...
        'plotScale', plot_scale, ...
        'avgMode', avgMode, ...
        'r_factor', r_factor,...
        'plotStates',plotStates,...
        'plotSegment',plotSegment,...
        'trialsPerCondition',trialsPerCondition,...
        'setRandSeed',setRandSeed);
end

% Collect figure handles
fh = [f_unconst, f_const];