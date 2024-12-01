function TDavg = average(TD,varargin)
% Average method for TrajectoryData class
%
% Usage:
%   TDavg = average(TD,varargin) or TDavg = TD.average(varargin)
%
% This function averages the trajectories contained in the provided
% trajectory data object array.
%
% Inputs:
%   TD      TrajectoryData class
%
% Optional Inputs:
%   N               Number of data points to interpolate over.
%   nDim            Number of dimension for averaging. Default is 2
%   avgMode         Method to use to perform averaging. Three options:
%                       interp, samp, spatial. Interp: interpolates each 
%                       trajectory to the same length. Samp: Averages over 
%                       samples. Spatial: Averages according to a velocity 
%                       flow field.
%   sampQuantCutoff Quantile to determine sample cutoff ('samp' average 
%                       mode only). Default is 0.5
%   r_factor        Scale factor used for spatial averaging (fraction of 
%                       target-to-target distance). Default is 0.15
%   kinSource       The kinematic data to calculate averages for.
%   avgStates       Which task states to include in the average. If empty
%                       than the average is over all provided data.
%   alignStates     Align method when using task states information: start, 
%                       end, or [] (default). If empty the data is aligned 
%                       to start and end of the given task states.
%                       Otherwise the start or end of the trial.
%
% Outputs:
%   TDavg      Average trajectoryData object
%
% Authors:  Alan Degenhart and Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional arguments
N = 100; % Number of data points to interpolate over.
nDim = 2; % Currently only use 2D
avgMode = 'interp';     % Method to use to perform averaging.
sampQuantCutoff = 0.5;  % Quantile to determine sample cutoff ('samp' average mode only)
r_factor = 0.15;  % Scale factor used for spatial averaging (fraction of target-to-target distance)
kinSource = {'brain'};
avgStates = [];
alignStates = [];

assignopts(who,varargin);

% Find unique target codes and targets
tC = [TD.targetCode];
uniTarg = unique(tC);
nTarg = length(uniTarg);

% Create average trajectory data object array
TDavg = repmat(TrajectoryData,nTarg,1);

% Loop over targets and average data
for i = 1:nTarg
    % Get all trajectories to the target
    tcMask = (tC == uniTarg(i));
    TDtemp = TD(tcMask);
    
    nTrials = length(TDtemp);
    
    % Loop over kinematic sources
    for j = 1:length(kinSource)
        trajData = nan(nTrials,nDim,N);

        % Loop over all trials and interpolate trajectory
        posMat = cell(nTrials,1);
        nSamp = nan(nTrials,1);
        for k = 1:nTrials
            KD = TDtemp(k).getKinematicData('kinSource',kinSource{j});
            
            % Make sure kinematic data exists.  It may be possible for the
            % kinematic data object to be empty (e.g., the phasespace
            % marker drops out).
            if isempty(KD.time)
                continue
            end
            
            % Subselect the states to average over.
            if ~isempty(avgStates)
                % Identify the relevant states and the earliest onset and
                % latest offset
                indState = ismember(TDtemp(k).states,avgStates);
                if sum(indState)==0
                    continue
                end
                
                minT = min(TDtemp(k).stateOnset(indState == 1));
                maxT = max(TDtemp(k).stateOffset(indState == 1));
                
                % Create the time index for the data
                if isempty(alignStates)
                    index = find(KD.time>=minT & KD.time<maxT);
                elseif strcmp(alignStates,'endTrl')
                    index = find(KD.time>=minT);
                elseif strcmp(alignStates,'startTrl')
                    index = find(KD.time<maxT); 
                end
                tempPos = KD.pos(index,1:nDim);
            else 
                tempPos = KD.pos(:,1:nDim);     % Only use X and Y positions
            end
            if size(tempPos,1)>0
                nSamp(k) = size(tempPos,1);     % Number of time points to interpolate over
            end
            posMat{k} = tempPos;
        end
        nSampMax = max(nSamp);
        
        if isnan(nSampMax)
            continue
        end
        
        switch avgMode
            case 'interp'
                trajData = nan(nTrials,nDim,N);
                for k = 1:nTrials
                    tempPos = posMat{k};
                    x = linspace(1,N,nSamp(k));

                    % Loop over dimensions and perform interpolation
                    if size(tempPos,1)>1
                        for l = 1:nDim
                            trajData(k,l,:) = interp1(x,tempPos(:,l),1:N,'spline');
                        end
                    end
                end
                
                % Average trajectory data
                trajData = squeeze(mean(trajData,1,'omitnan'));
                
            case 'samp' % Do not time warp
                trajData = nan(nTrials,nDim,nSampMax);
                for k = 1:nTrials
                    tempPos = posMat{k};
                    % Loop over dimensions
                    if size(tempPos,1)>1
                        for l = 1:nDim
                            trajData(k,l,1:nSamp(k)) = tempPos(:,l);
                        end
                    end
                end
                
                % Find number of points
                nSampCutoff = round(quantile(nSamp,sampQuantCutoff));
                trajData = trajData(:,:,1:nSampCutoff);
                
                % Average
                trajData = squeeze(mean(trajData,1,'omitnan'));
            case 'spatial'  % Spatial average (flow field)
                % Collect trajectory data
                trajData = nan(nTrials,nDim,nSampMax);
                for k = 1:nTrials
                    tempPos = posMat{k};
                    % Loop over dimensions
                    if size(tempPos,1)>1
                        for l = 1:nDim
                            trajData(k,l,1:nSamp(k)) = tempPos(:,l);
                        end
                    end
                end
                
                trajData = util.spatialAverage(trajData, 'r_factor', r_factor);
        end
        
        % Put data into new KinematicData object
        KD = KinematicData();
        KD.time = (1:size(trajData,2))';
        KD.pos = trajData';
        KD.source = kinSource{j};
        TDavg(i) = TDavg(i).setKinematicData(KD);
    end
    
    % Place average data into new trajectory object array
    TDavg(i).successful = true;
    TDavg(i).startPos = TDtemp(1).startPos;
    TDavg(i).targPos = TDtemp(1).targPos;
    TDavg(i).targSize = TDtemp(1).targSize;
    TDavg(i).targetCode = TDtemp(1).targetCode;
    TDavg(i).averaged = 1;
end