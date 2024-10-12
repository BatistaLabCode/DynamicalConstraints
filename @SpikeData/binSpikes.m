function [S,spikeCounts,binTimes] = binSpikes(S,binWidth,t0,tOffset,window)
% binSpikes         Bin spike data method for SpikeData class
%
% Inputs:
%   S           SpikeData class object
%   binWidth    Bin width (ms) to bin spike count data over
%   t0          Time point to bin relative to
%   tOffset     Time shift to apply to spike times
%   window      Time window to bin data over.

% Probably should add parse_varargin to allow variable arguments to be
% specified.  Both 'binType' and 'window' should be optional inputs, and
% might be specified independently.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

if ~S(1).multiBit
    % Get spike time cell array
    spikeTimes = S.getSpikeTimes(tOffset);
    
    % Generate bin times
    [binTimes,binOnset,binOffset] = util.generateBins(binWidth,t0,window);
    
    % Initialize spike count matrix
    nUnits = length(S);
    nBins = length(binOnset);
    spikeCounts = nan(nUnits,nBins);
    
    % Loop over trials and bin data
    for i = 1:nUnits
        st = spikeTimes{i};
        % Loop over bins and get spike counts -- exclude data from t ==
        % binOnset and include data from t == binOffset.  This is done
        % because we typically want to bin in a causal manner, and the bin
        % times are associated with the end of the bin.  Thus, a bin time
        % of 0 will include t=0.
        for j = 1:nBins
            spikeCounts(i,j) = sum((st > binOnset(j)) & (st <= binOffset(j)));
        end
        
        S(i).binSpikeCount = spikeCounts(i,:)'; % Transpose here to allow for array concatenating
        S(i).binTimes = binTimes';
    end
else %S(1).multiBit == true, spikes were binned online using Synapse
    binIdx = (S(1).binTimes-tOffset)>=window(1) & ...
        (S(1).binTimes-tOffset)<=window(2); %offset binTimes within window?
    binTimes = S(1).binTimes(binIdx)-tOffset;
    spikeCounts = [S.binSpikeCount];
    spikeCounts = spikeCounts(binIdx,:)';
end