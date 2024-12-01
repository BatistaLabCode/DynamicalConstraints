% TD = predictDecodeState_GPFA(TD, P, varargin)
%
% predictDecodeState_GPFA       Run GPFA decoder on TrajectoryData object
%
% Usage:
%   [TD] = predictDecodeState_GPFA(TD,P)
%
% Inputs:
%   TD          TrajectoryData
%   P           Decoder parameter structure
%
% Optional Inputs:
%   spikeCountSrc Structure to collect the spiking data from
%   trunGPFA      Determine whether or not to truncate the trajectory data
%   TT            Condition to apply orthonormalization
%   predictPos    Update the cursor position
%   useT          Use recorded time vectors
%   timeStep      Offset for the temporal alignment
%   noZeroPad     Truncate the decoded spikes to the number of ON channels 
% 
% Outputs:
%   TD          Updated TrajectoryData
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function TD = predictDecodeState_GPFA(TD, P, varargin)

% Parse optional arguments
spikeCountSrc = 'GPFA';  % Either 'GPFA' or 'decodeSpikeCounts'
trunGPFA = 0; % Determine whether or not to truncate the trajectory data
TT = [];  % Orthonormalization matrix
predictPos = true;
useT = false; 
timeStep = []; % Offset for the temporal alignment
noZeroPad = 0; % Truncate the decoded spikes to the number of channels that are on

assignopts(who, varargin);

nTrials = length(TD);

% Get GPFA parameters
nBins = P.nBins;
xDim = P.xDim;
M = P.M;
CRinv = P.CRinv;
d = P.d;
W = P.W;    % Decoding weights
c = P.c;    % Decode kinematic offset

% Loop over trials
for i = 1:nTrials
    % Get spike counts
    switch spikeCountSrc
        case 'GPFA'  % Spike counts binned using GPFA (offline)
            yAll = TD(i).GPFA.y;
        case 'decodeSpikeCounts'  % Spike counts binned online
            % Note: in this case spike counts will need to start from the
            % beginning of the trial, which isn't always the case with the
            % TrajectoryData object. The values may be slightly 
            yAll = TD(i).decodeSpikeCounts';
            if noZeroPad
                yAll(2:2:end,:) = [];
                yAll(P.exclCh,:) = [];
            end
        case 'binnedSpikes'
            yAll = TD(i).binnedSpikes;
            useT = true;
            if noZeroPad
                yAll(2:2:end,:) = [];
                yAll(P.exclCh,:) = [];
            end
    end
    
    % Initialize
    t1 = zeros(xDim,nBins);
    
    % Loop over timepoints
    nSamp = size(yAll,2);
    x = nan(2,nSamp);
    u = nan(xDim,nSamp);
    for j = 1:nSamp
        % Get spike counts for current timestep
        y = yAll(:,j);
        
        % Shift elements of t1
        t1(:,1:(nBins-1)) = t1(:,2:end);
        
        % Subtract off mean from current timestep and reshape
        dif = y - d;
        t1(:,end) = CRinv * dif;
        term1 = reshape(t1,xDim*nBins,1);
        
        % Calculate final smoothed latent state and predicted kinematics
        u(:,j) = M * term1;
%         if strcmp(xSpec,'xorth') & ~isempty(TT)
%             x(:,j) = W*TT*u(:,j) + c;
%         else
            x(:,j) = W*u(:,j) + c;
%         end
    end
    
    % Calculate time vector and truncate trajectory
    if useT
        t = TD(i).decodeTime;          
    else
        t = (1:nSamp)*P.binwidth;
    end
    if isempty(timeStep)
        tMask = (t >= (TD(i).trajOnset)) & (t <= TD(i).trajOffset);
    else
        tMask = (t > (TD(i).trajOnset - timeStep*TD(i).decoderBinWidth))...
            & (t <= TD(i).trajOffset);
    end
    t = t(tMask);
    x = x(:,tMask);
    
    % Invert the y-axis if using online binned spike counts.
    if ismember(spikeCountSrc, {'decodeSpikeCounts','binnedSpikes'}) & ~noZeroPad
        x(2,:) = -x(2,:);
    end
    
    if isempty(TD(i).GPFA)
        TD(i).GPFA = GPFAData;
    end
    
    % Orthonormalize the data
    if ~isempty(TT)
        TD(i).GPFA.xorth = TT*u;
        % Save the GPFA data
        if trunGPFA
            TD(i).GPFA.xorth = TD(i).GPFA.xorth(:,tMask);
        end
    else % Note this will save a value, but the data will not be orthonormalized
        % Save the GPFA data
        if trunGPFA
            TD(i).GPFA.xorth = u(:,tMask);
        else
            TD(i).GPFA.xorth = u;
        end
    end
    
    % Put data into TD structure
    if predictPos
        TD(i).brainKin.time = t';
        TD(i).brainKin.pos = x';
    end
    
    % Save the GPFA data
    if trunGPFA
        TD(i).GPFA.xsm = u(:,tMask);
        TD(i).GPFA.onsetIdx = 1;
    else
        TD(i).GPFA.xsm = u;
        TD(i).GPFA.onsetIdx = find(tMask,1,'first');
    end
    TD(i).GPFA.trialId = TD(i).trialID;
    TD(i).GPFA.T = size(TD(i).GPFA.xsm,2);
    TD(i).GPFA.onset = TD(i).trajOnset;
    TD(i).GPFA.binWidth = TD(i).decoderBinWidth;  
end