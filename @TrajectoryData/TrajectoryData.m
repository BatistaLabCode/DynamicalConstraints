classdef TrajectoryData
    properties
        % Meta info (keep at top)
        subject = [];       % Suject ID
        date = [];          % Date
        trialName = [];     % Trial name
        trialID = [];       % Unique trial identifier
        successful = false; % Specify if trajectory is successful or not
        states = [];        % States present in given trajectory
        stateOnset = [];    % Onset time for each state
        stateOffset = [];   % Offset time for each state
        
        % Hand Kinematics
        goCue = [];             % Time of go cue
        trajOnset = [];         % Onset of trajectory relative to the start of the trial
        trajOffset = [];        % Offset of trajectory relative to the start of the trial
        targetAcquired = [];    % Time target acquired
        trialLen = [];          % Total trial length (not just the selected states)
        
        % Brain control kinematics
        brainKin = KinematicData();
                
        % Target data
        startPosDefined = [];   % Specifies if the start position has been defined
        startPos = [];          % Trajectory start
        targPos = [];           % Target position
        targSize = [];          % Target size
        targetCode = [];        % Unique target code
        intTargPos = [];        % Target position for intermediate states
        intTargSz = [];         % Target size for intermediate states
        tube = TubeObject();    % Tube
        
        % Neural Data
        GPFA = GPFAData();      % GPFA neural trajectory
        
        % Decoder Data
        decodeSpikeCounts = []; % Binned spike times used during online decoding
        decodeTime = [];        % RT time (for online binned spikes)
        decodeState = [];       % Raw decoder output
        decoderName = [];       % Name of decoder used for given trial
        decoderNum = [];        % Number of decoder used for given trial
        decoderBinWidth = [];   % Bin width used for decodeSpikeCounts
        
        % Flags
        averaged = false;   % Specify if data is averaged (across targets)
        normalized = false; % Normalized with respect to the center target
    end % End of properties
    
    methods
        function TD = TrajectoryData(); end
        
        % Save/load methods
        save(TD,pathName,fileBase)      % Save data
        [TD,directoryName] = load(TD,directoryName)
        
        % Processing methods
        TD = normalize(TD,explicitStartPos);
        TDavg = average(TD,varargin);

        % Data access methods
        KD = getKinematicData(TD,varargin);
        TD = setKinematicData(TD,KD);
        TD = setTargetCode(TD,tC);
        TD = checkTrialSucces(TD);
        
        % Plotting methods
        plot(TD,varargin);      % Plot trajectories
        
    end % End of methods
end % End of classdef
