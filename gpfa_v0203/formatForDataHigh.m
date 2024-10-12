function [Dall,Davg] = formatForDataHigh(G,varargin)

conditionPrefix = 'Target';
alignOption = 'time0';

assignopts(who,varargin);

TD = G.TD;
TDavg = G.TDavg;
GP = [TD.GPFA];
GPavg = [TDavg.GPFA];

nTrials = length(TD);
nTrialsAvg = length(TDavg);

tC = [TD.targetCode];
tCavg = [TDavg.targetCode];
uniTC = unique(tC);
nTarg = length(uniTC);
col = interpColor(nTarg,1,'hsv',1);

% Initialize structure array for DataHigh data
D.type = 'traj';
D.epochStarts = [1];
D.epochColors = [0 0 0];
D.condition = '';
D.data = [];
Davg = repmat(D,nTrialsAvg,1); % 2 directions per target combination
Dall = repmat(D,nTrials,1);

% Get DataHigh structure for average GPFA data
for i = 1:nTrialsAvg
    % Get average trajectory
    x = GPavg(i).xorth;
    onset = GPavg(i).onsetIdx;
    Davg(i).data = x(:,onset:end);
    
    % Set color and condition string
    Davg(i).epochColors = col(i,:);
    condStr = sprintf('Cond %d',tCavg(i));
    Davg(i).condition = condStr;
end

% Get DataHigh structure for all trials
for i = 1:nTrials
    % Get trajectory
    x = GP(i).xorth;
    
    % Get onset
    switch alignOption
        case 'time0'
            onset = 1;
        case 'moveOnset'
            onset = GP(i).onsetIdx;
    end

    % Get data
    Dall(i).data = x(:,onset:end);
    % Set color and condition string
    Dall(i).epochColors = col(tC(i),:);
    condStr = sprintf('%s %d',conditionPrefix,tC(i));
    Dall(i).condition = condStr;
end