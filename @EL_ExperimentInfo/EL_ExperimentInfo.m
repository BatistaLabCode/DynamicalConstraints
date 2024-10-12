% Energy Landscape Experiment Info class (@EL_ExperimentInfo)
%
% This is a s
classdef EL_ExperimentInfo
    
    properties
        % Dataset properties
        subject = [];
        dataset = [];
        dataDir = [];       % Directory where data is located
        saveNameBase = [];
        dataPath = [];
        analysisPath = [];
        
        % GPFA properties
        exclCh = [];        % Channels excluded from GPFA
        jaggedDim = [];     % Identified jagged dimensions in GPFA
        onsetMode = [];     % Mode used to identify onset of trajectories
        gpfaRunID = [];     % Unique GPFA run identifier
        xSpec = [];         % Type of latents used ('xsm' or 'xorth')
        causal = [];        % Use causal GPFA
        causalMode = []     % Causal mode used ('full' or 'rt')
        causal_Tmax = [];   % Duration used to estimate filter coeffs
        causal_nBins= [];   % Number of time bins to truncate smoothing coeffs to
        
        % Projection/Decoder properties
        targPair = [];      % Target pair for identified asymmetry
        invPri = [];        % Invert primary axis
        invSec = [];        % Invert secondary axis
        isotropicWVal = []; % Smoothing constant
        A = [];
        x_c = [];
        p_c = [];
        P = [];             % Selected projection
        OptProj = [];       % All available optimized projections

    end
    
    methods
        function E = EL_ExperimentInfo(); end
            
    end
end