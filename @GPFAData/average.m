function GPavg = average(GP,varargin)
% GPavg = average(GP,varargin) or GPavg = GP.average(varargin)
% Average method for GPFAData class
%
% Inputs:
%   GP      GPFAData class
%
% Optional Inputs:
%   alignMode        Align method: start (default), end, or index 
%   sampQuantCutoff  Quantile to determine sample cutoff ('samp' 
%                    average mode only). Default is 0.5
%   nTS              Specify the number of data points to average over. 
%                    Currently only valid for the start and end align modes.
%
% Authors:  Alan Degenhart and Erinn Grigsby
% Created:  2017.02.16
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


% Optional arguments
alignMode = 'start';    % Align method (start, end, or index)
sampQuantCutoff = 0.5;  % Quantile to determine sample cutoff ('samp' 
                        % average mode only)
nTS = [];               % Specify the number of data points to average over. 
                        % Only valid for the start and end align modes.

assignopts(who, varargin);

% Initialize structure array for average trajectories
nTrials = length(GP);
nDimxSm = size(GP(1).xsm,1);
nDimxOrth = size(GP(1).xorth,1);
nUnits = size(GP(1).y,1);

% Get trajectory (movement) onset indices
alignInd = [GP.onsetIdx];

% If the alignment mode is set to 'index', ensure the appropriate number of
% indices are provided.
if strcmp(alignMode,'index')
    assert(length(alignInd) == nTrials, ...
        'Length of alignment indices does not equal the number of available trials.')
end

% Determine maximum number of pre-onset and post-onset samples and the
% total number of samples for the entire dataset.
TS = [GP.T];
switch alignMode
    case 'start'
        if isempty(nTS)
            nTS = max(TS);
        end
        indAlign = 1;
    case 'end'
        if isempty(nTS)
            nTS = max(TS);
        end
        indAlign = nTS;
    case 'index'
        preOnsetSamp = alignInd - 1;
        postOnsetSamp = TS - preOnsetSamp;
        nTS = max(preOnsetSamp) + max(postOnsetSamp);
        indAlign = max(preOnsetSamp) + 1;
end

% Initialize data matrices
xSmMat = nan(nDimxSm,nTS,nTrials);
xOrthMat = nan(nDimxOrth,nTS,nTrials);
yMat = nan(nUnits,nTS,nTrials);

% Loop over trials and get data
for i = 1:nTrials
    % Determine onset and offset index depending on method
    switch alignMode
        case 'start'
            % Align to start of trial/trajectory
            onset = 1;
            offset = TS(i);
        case 'end'
            % Align to end of trial/trajectory
            offset = nTS;
            onset = nTS - TS(i) + 1;
        case 'index'
            % Align to index provided
            onset = indAlign - preOnsetSamp(i);
            offset = indAlign + postOnsetSamp(i) - 1;
    end
    
    % Get data for current trial
    xSmMat(:,onset:offset,i) = GP(i).xsm;
    xOrthMat(:,onset:offset,i) = GP(i).xorth;
    yMat(:,onset:offset,i) = GP(i).y;
end

% Average across trials
xSm = nanmean(xSmMat,3);
xOrth = nanmean(xOrthMat,3);
y = squeeze(nanmean(yMat,3));
nanMat = xSmMat(1,:,:);
nanMat = squeeze(nanMat);
nTrials = sum(~isnan(nanMat),2);

% % In order to keep the average from being too noisy, remove those time
% % points with fewer than 50% of the total trials contributing to the
% % average
% sampleMask = nTrials >= floor(length([GP.trialId]));

% Find number of points
sampleMask = round(quantile(~isnan(nanMat'),sampQuantCutoff));

% Find index corresponding to first instance of sufficient data.  Remove
% invalid average timepoints and adjust onset accordingly.
trajOnset = find(sampleMask == 1,1,'first');
trajOffset = find(sampleMask == 1,1,'last');
xSm = xSm(:,trajOnset:trajOffset);
xOrth = xOrth(:,trajOnset:trajOffset);
y = y(:,trajOnset:trajOffset);
indAlign = indAlign - trajOnset + 1;

% Get bin width.  This should be the same for all trials
binWidth = [GP.binWidth];
uniBinWidth = unique(binWidth);
if length(uniBinWidth) > 1
    warning('GPFA bin width is not consistent across trials.')
end

% Pack up into GPFAData object
GPavg = GPFAData();
GPavg.trialId = [GP.trialId];           % All trial IDs in average
GPavg.T = size(xSm,2);
GPavg.y = y;
GPavg.Vsm = [];
GPavg.VsmGP = [];
GPavg.xsm = xSm;                        % Averaged neural trajectory
GPavg.xorth = xOrth;
GPavg.nAvgTrials = nTrials(sampleMask == 1);  % Number of trials in average per timestep
GPavg.avgAlignInd = indAlign;           % Index corresponding to alignment point
GPavg.binWidth = uniBinWidth(1);        % Bin width used for GPFA
GPavg.onsetIdx = indAlign;
GPavg.onset = (indAlign-1) * uniBinWidth;
GPavg.averaged = true;