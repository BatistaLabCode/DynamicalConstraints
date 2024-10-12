% el.tube.get_IT_data  Get intermediate target data for constrained path
% (tube) analysis.
%
% 

function [IT] = get_IT_data(IT,varargin)
dir_TD = [];
centerPos = [0 0 0];

assignopts(who,varargin);

% Automatically identify filenames for datasets to process/analyze
A.trajDataDir = fullfile(dir_TD(1).base,'translated','trajectoryData');

% Get data: intermediate target - no tube
TD = TrajectoryData().load(fullfile(A.trajDataDir, ...
    IT.unconstrainedBlockDir));
[TDnorm,~] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

% Remove spike data and GPFA data
for i = 1:length(TDnorm)
    TDnorm(i).spikes = [];
    TDnorm(i).GPFA = [];
end

% Add TrajectoryData object to IT structure
IT.TDunconstrained = TDnorm;

% Get data: two target intuitive
nBlocks = length(IT.ttIntBlockDir);
if nBlocks~=0
    for i = 1:nBlocks
        % Load data and preprocess
        TD = TrajectoryData().load(fullfile(A.trajDataDir,...
            IT.ttIntBlockDir{nBlocks}));
        [TDnorm,~] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);
        
        % Remove the spike data and GPFA data
        for j = 1:length(TDnorm)
            TDnorm(j).spikes = [];
            TDnorm(j).GPFA = [];
        end
        
        % Determine the number of target pairs in the trajectory data and
        % separate them into different cell -TODO
        
        IT.TD_tt_int{i} = TDnorm;
    end
end

% Get data: two target rotated
nBlocks = length(IT.ttRotBlockDir);
if nBlocks~=0
    for i = 1:nBlocks
        % Load data and preprocess
        TD = TrajectoryData().load(fullfile(A.trajDataDir,...
            IT.ttRotBlockDir{nBlocks}));
        [TDnorm,~] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);
        
        % Remove the spike data and GPFA data
        for j = 1:length(TDnorm)
            TDnorm(j).spikes = [];
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

    % Remove spike data and GPFA data
    for j = 1:length(TDnorm)
        TDnorm(j).spikes = [];
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

% Get data: intermediate target - no tube (washout)
if ~isempty(IT.unconstrainedWashoutBlockDir)
    % Load data and preprocess
    TD = TrajectoryData().load(fullfile(A.trajDataDir, ...
        IT.unconstrainedWashoutBlockDir));
    [TDnorm,~] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);
    
    % Remove spike data and GPFA data
    for i = 1:length(TDnorm)
        TDnorm(i).spikes = [];
        TDnorm(i).GPFA = [];
    end
    
    % Add TrajectoryData object to IT structure
    IT.TDunconstrainedWashout = TDnorm;
end
