function [tPos,sz] = getTargetPosition(trial,stateName,targetName,varargin)
% getTargetPosition         Get position for specified target
%
% This function returns the target position of the specified target for the
% indicated state.
%
% Inputs:
%   trial       Trial structure
%   stateName   Name of state with the desired target
%   targetName  Name of target to return position for
%
% Optional Input:
%   cursorIndex Position source data. 
%
% Outputs:
%   tPos        Position of specified target
%   sz          Size of specified target (radius?)
%
% Author:  Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional arguments
cursorIndex = 1;

assignopts(who,varargin);

SI = trial.Parameters.states;

% Get state index for specified state
stateIdx = SI.findStateIndex(stateName);

% Get target position
targetNames = SI(stateIdx).Cursor(cursorIndex).name;
targIdx = find(strcmp(targetName,targetNames),1,'first');
tPos = SI(stateIdx).Cursor(cursorIndex).window(targIdx,1:3);
sz = SI(stateIdx).Cursor(cursorIndex).window(targIdx,4);