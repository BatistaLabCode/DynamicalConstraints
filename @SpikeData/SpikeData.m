classdef SpikeData < HSTData
    properties
        channel = [];           % Channel number
        sort = [];              % Sort ID
        source = [];            % Spike source
        timeStampUnits = [];    % Time stamp units
        
        spikeTimes = [];        % Array of spike times
        waveform = [];          % Matrix of spike waveforms
        threshold = [];         % Spike detection threshold
        
        binSpikeCount = [];     % Binned spike counts
        binTimes = [];          % times for binned spikes
        
        multiBit = false;       % multiple bits per sort
        
        includesSort = [];      % Data includes sort
        includesWaveforms = []; % Data includes waveforms
    end % End of properties
    
    methods
        function S = SpikeData(); end       % Class constructor 
        
        % Plotting methods
        raster(S,tOffset,window,varargin)   % Plot raster
        [pAvg,p,t] = psth(S,varargin)       % Plot PSTH
        
        % Processing methods
        [S] = correctUnits(S)                 % Convert spike times to ms
        [spikeTimes,nSpikes] = getSpikeTimes(S,tOffset,window)  % Get spike time data
        [S,spikeCounts,binTimes] = binSpikes(S,binWidth,t0,tOffset,window)
        [S] = combineSorts(S,varargin);
    end % End of methods
end % End of classdef