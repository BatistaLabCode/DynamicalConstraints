% findDirContents       Find items within directory
%
% [dList] = findDirContents(baseDir,itemName)
%
% This function returns the valid items in the directory specified by
% BASEDIR with the name specified by ITEMNAME.
%
% Usage:
%   util.findDirContents(baseDir,itemName)
% 
% Inputs:
%   baseDir         Base directory to identified all dependent content
%   itemName        File identifier to subselect for specific file content.              
% 
% Outputs:
%   dList           List of files with valid elements.
%
% Author:   Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com
function dList = findDirContents(baseDir,itemName)

% Get valid directory elements
D = dir(baseDir);
dName = {D.name};
fMask = contains(dName,itemName);

% Get name(s) for valid elements
dName = {D.name};
dList = dName(fMask);