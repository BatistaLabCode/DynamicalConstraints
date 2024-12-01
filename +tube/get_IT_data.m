% tube.get_IT_data  Get intermediate target data for constrained path
% (tube) analysis.
%
% Usage:
%   tube.get_IT_data(IT)
%
% 
% This will automatically collect and fill in the intermediate target 
% experiment data and success rates based on the filepath locations. 
%
% Input:
%   IT          IntTargExp object
%
% Optional Input:
%   dir_TD      Directory where the TrajectoryData is saved
%   centerPos   Center of the workspace
%   intTargNum  Intermediate target number
%   blockSize   Largest number of trials between sucess rate evaluation
%
% Optional Input:
%   IT          IntTargExp object
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function [IT] = get_IT_data(IT,varargin)
dir_TD = [];
centerPos = [0 0 0];
intTargNum = 1;
blockSize = 50;

assignopts(who,varargin);

% Automatically identify filenames for datasets to process/analyze
A.trajDataDir = fullfile(dir_TD(1).base,'translated','trajectoryData');

% Get data: intermediate target - no tube
TD = TrajectoryData().load(fullfile(A.trajDataDir, ...
    IT.unconstrainedBlockDir));
[TDnorm,~] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

% Remove the GPFA data
for i = 1:length(TDnorm)
    TDnorm(i).GPFA = [];
end

% Add TrajectoryData object to IT structure
IT.TDunconstrained = TDnorm;

% Get data: two target rotated
nBlocks = length(IT.ttRotBlockDir);
if nBlocks~=0
    for i = 1:nBlocks
        % Load data and preprocess
        TD = TrajectoryData().load(fullfile(A.trajDataDir,...
            IT.ttRotBlockDir{nBlocks}));
        [TDnorm,~] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

        % Remove the GPFA data
        for j = 1:length(TDnorm)
            TDnorm(j).GPFA = [];
        end

        % Save to the structure
        IT.TD_tt_rot{i} = TDnorm;
    end
end

% Get data: intermediate target - tube
nBlocks = length(IT.constrainedBlockDir);
IT.TDconstrained = cell(nBlocks,1);
for i = 1:nBlocks
    % Load data and preprocess
    TD = TrajectoryData().load(fullfile(A.trajDataDir, ...
        IT.constrainedBlockDir{i}));
    [TDnorm,~] = util.preprocessGridTaskTrajData(TD, ...
        'centerPos',centerPos);

    % Remove GPFA data
    for j = 1:length(TDnorm)
        TDnorm(j).GPFA = [];
    end

    % Determin the number of tubes in the trajectory data and separate them
    % into different cells
    TU = [TDnorm.tube];
    TU = TU.setTubeCode;
    uniTU = unique([TU.tubeCode],'stable');
    for k = 1:length(uniTU)
        % Create a tube mask
        mask = uniTU(k) == [TU.tubeCode];

        % Add TrajectoryData object to IT structure
        IT.TDconstrained{i}{k,1} = TDnorm(mask);
    end
end

% Iterate through the TDconstrained structure so there is one tube size per
% cell row.
tempTDconstrained = [];
for i = 1:size(IT.TDconstrained,1)
    tempTDconstrained = [tempTDconstrained; IT.TDconstrained{i}];
end
IT.TDconstrained = tempTDconstrained;

% Find the start position of the constrained trials. Here we use the first 
% constrained tube block. 
TD = IT.TDconstrained{1};
startPos = [TD.startPos]';
uniStartPos = unique(startPos,'rows');
nStartPos = size(uniStartPos,1);

% Get unique start positions and associated tubes
TUall = repmat(TubeObject(),nBlocks,nStartPos);
tubeTrialsAll = nan(nBlocks,nStartPos);
tubeTrialsSuc = nan(nBlocks,nStartPos);
tubeRadius = nan(nBlocks,nStartPos);
for i = 1:nBlocks
    % Get trajectory and tube data, set tube code.
    TDtemp = IT.TDconstrained{i};
    TU = [TDtemp.tube];
    TU = TU.setTubeCode;
    
    % Find start positions and get tube objects for each position
    tubeStartPos = [TDtemp.startPos]';
    for j = 1:nStartPos
        % Get tube object
        startPosMask = ismember(tubeStartPos,uniStartPos(j,:),'rows');
        TUtemp = TU(startPosMask);
        TUall(i,j) = TUtemp(1);
        
        % Get tube radius
        r = TUtemp(1).window(:,1);
        r = r(~isnan(r));
        tubeRadius(i,j) = r(1);
        TUall(i,j).radius = r(1);
        
        % Calculate acutal success rate
        TDtube = TDtemp(startPosMask);
        tubeTrialsAll(i,j) = length(TDtube);
        tubeTrialsSuc(i,j) = sum([TDtube.successful]);
    end
end
tubeSR = tubeTrialsSuc./tubeTrialsAll;

% Calculate the predicted success rate of the unconstrained trials through
% the tube.
mask = ismember([IT.TDunconstrained.startPos]',startPos,'rows');
TDuncon = IT.TDunconstrained(mask==1);
condSR = tube.calcExpectedSR(TDuncon,TUall(1),tubeRadius)

% Add data to IT structure
IT.successRate = tubeSR;
IT.predictedSuccessRate = condSR;
IT.constrainedTubeRadius = tubeRadius;
IT.startTargPos = unique(startPos,'rows');
IT.intermediate_target_num = intTargNum;
IT.block_size = blockSize;
IT.tubeObject = TUall;
