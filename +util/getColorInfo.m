% getColorInfo      Get color information for associated target pair
%
% [colInfo] = util.getColorInfo(C, startPos, endPos)
%
% Usage:
%   util.getColorInfo(C, startPos, endPos)
% 
% Inputs:
%   C           'ColMat' structure of position-to-color definitions
%   startPos    Start target position
%   endPos      End target position
%
% Output:
%   colInfo{1} - light color
%   colInfo{2} - dark color
%   colInfo{3} - target pair string
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [colInfo] = getColorInfo(C,startPos,endPos)

% Get number of start/end positions to find colors for
nPos = max(size(startPos, 1), size(endPos, 1));

startPos = round(startPos);
endPos = round(endPos);

colInfo = cell(nPos,3);
for i = 1:nPos
    % Find associated target codes for start and end position.  To do this,
    % first round the target position to find the closest match
    targPos = round(C.targPos);
    tcStart = [];
    tcEnd = [];
    if ~isempty(startPos)
        tcStart = find(ismember(targPos,startPos(i,:),'rows'));
    end
    if ~isempty(endPos)
        tcEnd = find(ismember(targPos,endPos(i,:),'rows'));
    end
    
    % Handle case where start or end target isn't provided.  If at least
    % one target is available, this function should be able to return a
    % valid color.
    if isempty(tcEnd) && ~isempty(tcStart)  % Only start target specified
        condStr = sprintf('Start: %s',C.targLabel{tcStart});
        colMask = ismember(C.startTC',tcStart);
    elseif isempty(tcStart) && ~isempty(tcEnd)  % Only end target specified
        condStr = sprintf('End: %s',C.targLabel{tcEnd});
        colMask = ismember(C.endTC',tcEnd,'rows');
    else  % Both start and end target specified
        condStr = sprintf('%s -> %s',C.targLabel{tcStart},C.targLabel{tcEnd});
        colMask = ismember([C.startTC' C.endTC'],[tcStart tcEnd],'rows');
    end
    
    % Throw error if target was not found
    if isempty(tcStart) && isempty(tcEnd)
        error('Could not find matching color information for provided target.')
    end

    % Check to see if target pair is found in the list of defined start/end
    % target pairs for the task.  If so, get the associated color pair.  If
    % not, get the 'undefined' target pair color.
    if sum(colMask) > 0
        cL = C.col_light(colMask,:);
        cD = C.col_dark(colMask,:);
    else
        cL = C.undef_col_light;
        cD = C.undef_col_dark;
    end
    
    % Add color and condition string info to cell
    colInfo(i,:) = {cL,cD,condStr};
end