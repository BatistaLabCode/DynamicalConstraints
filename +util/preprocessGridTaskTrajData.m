function [TD,targInfo] = preprocessGridTaskTrajData(TD,varargin)
% [TD] = preprocessGridTaskTrajData(TD,startPos) 
%
% Preprocess trajectory data objects from grid task data.
%
% This function preprocesses the TrajectoryData objects in TD to ensure
% that the resultant data is normalized appropriately and has the correct
% target code info.
%
% Input:
%   TD             TrajectoryData object
%
% Optional Input:
%   centerPos      Location of the center ofthe workspace.
%
% Output:
%   TD             Modified TrajectoryData object
%   targInfo       Structure with target information and target codes
%
% Author:       Alan D. Degenhart and Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional arguments
centerPos = [0 0 0];

assignopts(who,varargin);

% Normalize TrajectoryData object
TD = TD.normalize(centerPos);

% Get target codes for start and end targets
startPos = [TD.startPos]';
endPos = [TD.targPos]';
[tcStart,M] = util.pos2targetCode([startPos; endPos]);
tcEnd = tcStart(size(startPos,1)+1:end);
tcStart = tcStart(1:size(startPos,1));
uniTCStart = unique(tcStart);
uniTCEnd = unique(tcEnd);
numTCStart = length(uniTCStart);
numTCEnd = length(uniTCEnd);

% Get unique target combinations
tcComb = [tcStart tcEnd];
targComb = unique(tcComb,'rows');
nComb = nchoosek(length(M.uniTC),2);

% Remove invalid starting points
if length(M.uniTC)-numTCStart>=2
    nComb = nComb - nchoosek(length(M.uniTC)-numTCStart,2);
end
% Remove invalid ending points
if length(M.uniTC)-numTCEnd>=2
    nComb = nComb - nchoosek(length(M.uniTC)-numTCEnd,2);
end
if numTCStart == numTCEnd
    nComb = ceil(size(targComb,1)/2);
end

% Update target code to reflect target pairs.  The total number of target
% codes will be nTarg * 2, and will be ordered such that odd target codes
% represent A->B trajectories and even target codes represent B->A
% trajectories (e.g., targ 1>2 will have a target code of 1, and targ 2->1
% will have a target code of 2.
tcCombList = nan(nComb*2,2);
tcMask = true(size(targComb,1),1);
tempTC = 0;
for i = 1:size(targComb,1)
    if tcMask(i)
        % Get trials from A->B
        targMask = ismember(tcComb,targComb(i,:),'rows');
        tempTC = tempTC+1;
        tcCombList(tempTC,:) = targComb(i,:);
        TD(targMask) = TD(targMask).setTargetCode(tempTC);
        
        % Get trials form B->A
        targMask = ismember(tcComb,fliplr(targComb(i,:)),'rows');
        tempTC = tempTC+1;
        tcCombList(tempTC,:) = fliplr(targComb(i,:));
        TD(targMask) = TD(targMask).setTargetCode(tempTC);
        
        % Remove the flipped target combination from the targComb matrix.
        flipMask = ismember(targComb,fliplr(targComb(i,:)),'rows');
        tcMask(flipMask,:) = false;
    end
end

% Pack up target code mapping info
targInfo.startPos = startPos;
targInfo.endPos = endPos;
targInfo.tcStart = tcStart;
targInfo.tcEnd = tcEnd;
targInfo.tcCombList = tcCombList;
targInfo.M = M;