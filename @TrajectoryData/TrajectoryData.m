classdef TrajectoryData < HSTData
    properties
        % Meta info (keep at top)
        subject = [];       % Suject ID
        date = [];          % Date
        trialName = [];     % Trial name
        datasetName = [];   % Unique identifier for dataset 
        trialID = [];       % Unique trial identifier
        tag = [];           % Tag that can be used to disambiguate among data sets
        successful = false; % Specify if trajectory is successful or not
        states = [];        % States present in given trajectory
        stateOnset = [];    % Onset time for each state
        stateOffset = [];   % Offset time for each state
        
        % Specify source of control signal
        controlSource = '';     % Type of control (hand control/brain control)
        
        % Hand Kinematics
        kinTime = [];           % Time (kinematics)
        pos = [];               % Position
        vel = [];               % Velocity
        acc = [];               % Acceleration
        goCue = [];             % Time of go cue
        trajOnset = [];         % Onset of trajectory relative to the start of the trial
        trajOffset = [];        % Offset of trajectory relative to the start of the trial
        moveOnset = [];         % Movement onset
        moveOffset = [];        % Movement offset
        targetAcquired = [];    % Time target acquired
        trialLen = [];          % Total trial length (not just the selected states)
        handKin = KinematicData();
        forceKin = ForceData();
        perturbKin = PerturbationData();
        
        % Brain control kinematics
        brainControlPos = [];
        brainControlVel = [];
        brainKin = KinematicData();
                
        % Target data
        startPosDefined = [];   % Specifies if the start position has been defined
        startPos = [];          % Trajectory start
        targPos = [];           % Target position
        targSize = [];          % Target size
        targetCode = [];        % Unique target code
        intTargPos = [];        % Target position for intermediate states
        intTargSz = [];         % Target size for intermediate states
        cursorSize = [];        % Cursor size
        tube = TubeObject();    % Tube
        
        % Neural Data
        spikes = [];            % Spike data
        binTime = [];           % Binned spike time
        binnedSpikes = [];      % Binned spike data
        binnedChannel = [];     % Channel numbers for binned data
        binnedSort = [];        % Sort numbers for binned data
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

        % Data access methods
        KD = getKinematicData(TD,varargin);
        TD = setKinematicData(TD,KD);
        TD = setTargetCode(TD,tC);
        TD = setTag(TD,tagVal);
        
        % Plotting methods
        plot(TD,varargin);      % Plot trajectories
        S = plotProgress(TD);   % Plot success rate and acquisition time
        
        % Analysis methods
        [lag,F] = calcSpikeKinLag(TD,varargin);   % Calculate lag between spiking and kinematic data
    end % End of methods
end % End of classdef
