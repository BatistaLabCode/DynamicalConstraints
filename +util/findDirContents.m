function dList = findDirContents(baseDir,itemName)
% findDirContents       Find items within directory
%
% [dList] = findDirContents(baseDir,itemName)
%
% This function returns the valid items in the directory specified by
% BASEDIR with the name specified by ITEMNAME.
%
% Author:   Alan D. Degenhart
% Created:  2018.08.21
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Get valid directory elements
D = dir(baseDir);
dName = {D.name};
fMask = contains(dName,itemName);

% Get name(s) for valid elements
dName = {D.name};
dList = dName(fMask);