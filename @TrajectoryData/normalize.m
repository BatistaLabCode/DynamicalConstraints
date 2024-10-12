function TD = normalize(TD,explicitStartPos)
% normalize         Normalize trajectory data
%
% This function sets the target code information and normalizes trajectory
% data relative to the center target.  Unique target code IDs are
% determined for each target position.
%
% Author:  Alan D. Degenhart
% Date Created: 2015/01/14
% Last Updated: 2016/06/21
% Last Update:  Updated to normalize both brain and hand position data.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

explicitStartPosDefined = 0;
if nargin > 1
    explicitStartPosDefined = 1;
end

% Get number of trials
nTrials = length(TD);

% Get target size information
nDim = length(TD(1).targPos);

% Get target and start positions
targPos = nan(nTrials,nDim);
startPos = nan(nTrials,nDim);
calcStartPos = nan(nTrials,1);
for i = 1:nTrials
    calcStartPos(i) = ~TD(i).startPosDefined;
    
    % Get target and start positions
    targPos(i,:) = TD(i).targPos;
    if explicitStartPosDefined
        startPos(i,:) = explicitStartPos;
        calcStartPos(i) = 0;     % Do not need to calculate start pos
    elseif ~calcStartPos(i)
        startPos(i,:) = TD(i).startPos;
    end
    
end

% If needed, calculate average target position to define a start location
calcStartPos = logical(calcStartPos);
nCalcStartPos = sum(calcStartPos);
if nCalcStartPos > 0 % At least one trial needs a start position
    uniTarg = unique(targPos,'rows');
    definedStartPos = nanmean(uniTarg,1);
    startPos(calcStartPos,:) = repmat(definedStartPos,nCalcStartPos,1);
end
nTargDim = size(targPos,2);

% Initialize array for target positions, loop over all trials and get
% target info.
for i = 1:nTrials
    % Normalize target and cursor with respect to start target
    tempTargPos = targPos(i,:) - startPos(i,:);
    TD(i).targPos = tempTargPos';
    
    % Center start position if it was not defined.  If the start position
    % is not defined (either explicitly or automatically), it will *always*
    % be zero.
    if TD(i).startPosDefined
        TD(i).startPos = TD(i).startPos - startPos(i,:)';
    else
        TD(i).startPos = zeros(fliplr(size(startPos(i,:))));
    end
    
    % Get kinematic data.  Do this for both the hand data and the brain
    % data.
    kinSource = {'brain','hand','force','perturb'};
    for j = 1:length(kinSource)
        KD = TD(i).getKinematicData('kinSource',kinSource{j});
        
        % Get position data and re-center
        if ~isempty(KD.time)
            pos = KD.pos;
            nSamples = size(pos,1);
            
            % Check to make sure position is the same size as the target
            % information.  If not, zero-pad
            nKinDim = size(pos,2);
            if nKinDim < nTargDim
                pos(:,(nKinDim+1):nTargDim) = zeros(nSamples,(nTargDim - nKinDim));
            end
            pos = pos - repmat(startPos(i,:),nSamples,1);
            KD.pos = pos;
            TD(i) = TD(i).setKinematicData(KD);
        end
    end
    
    % Normalize the position data for the primary cursor
    primePos = TD(i).pos;
    TD(i).pos = primePos - repmat(startPos(i,1:size(primePos,2)),size(primePos,1),1);
    
    targPos(i,:) = tempTargPos;
    
    % Normalize intermediate target positions if included
    if ~isempty(TD(i).intTargPos)
        intTargPos = TD(i).intTargPos;
        intTargPos = intTargPos - repmat(startPos(i,:)',1,size(intTargPos,2));
        TD(i).intTargPos = intTargPos;
    end
    
    % Normalize the tube path and boundary if included
    if ~isempty(TD(i).tube.path)
        % Normalize the tube path
        pathTube = TD(i).tube.path;
        nSamples = size(pathTube,1);
        pathTube = pathTube - repmat(startPos(i,:),nSamples,1);
        TD(i).tube.path = pathTube;
        % Normalize the tube boundary. Flexible for 2D and 3D boundaries 
        boundTube = TD(i).tube.boundary;
        [boundL, tubeDimen] = size(boundTube);
        boundTube = boundTube - repmat(startPos(i,1:tubeDimen),boundL,1);
        TD(i).tube.boundary = boundTube;
    end
end

% Calculate target code.  A center position of [0 0 0] is used because the
% data is assumed to be centered at the origin.
tC = util.pos2targetCode(targPos,'centerPos',[0 0 0]);

% Loop over trajectories and update target code data
for i = 1:nTrials
    % Set unique target code
    TD(i).targetCode = tC(i);
    
    % Set normalized status
    TD(i).normalized = true;
end