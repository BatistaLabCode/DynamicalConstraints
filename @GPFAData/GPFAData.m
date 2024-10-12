classdef GPFAData < HSTData
    properties
        % Data
        trialId = [];       % Unique trial number
        T = [];             % Number of time steps
        y = [];             % Spike counts
        
        xsm = [];           % Latent state (smoothed, non-orthonormalized)
        Vsm = [];           % Variance (smoothed)
        VsmGP = [];         % Variance (smoothed-GP)       
        xorth = [];         % Orthonormalized latent state (smoothed)
        
        onset = [];         % Onset time
        onsetIdx = [];      % Index for onset event (movement onset, etc.)
        
        % Internal
        causal = [];        % Non-causal or causal GPFA
        averaged = false;   % Specify if data has been averaged
        nAvgTrials = [];    % Number of trials going into each timestep
        avgAlignInd = [];   % Index used for alignment
        notes = [];         % Placeholder to allow text to be saved
        
        % Parameters
        binWidth = [];      % Bin width used to estimate trajectories
    end  % End or properties
    
    methods
        function GP = GPFAData(); end;  % Class constructor
        
        % Access methods
        GP = setOnsetIdx(GP,tOnset,varargin);
        
        % Analysis methods
        GP = average(GP,varargin)   % Average GPFA trajectories
        
        % Plot methods
        F = plot(GP,varargin)           % Plot latent state vs. time
        F = trajectory(GP,varargin)     % Plot neural trajectory
        
    end  % End of methods
end  % End of classdef
        