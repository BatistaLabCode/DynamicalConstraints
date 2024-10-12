% Class definition for intermediate target experiment
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com
classdef IntTargExp
    properties
        % Experiment info
        subject = [];
        date = [];
        intermediate_target_num = [];
        
        % Information about specific blocks
        unconstrainedBlockDir = [];
        constrainedBlockDir = [];
        unconstrainedWashoutBlockDir = [];
        ttIntBlockDir = [];
        ttRotBlockDir = [];
        
        % Trajectory data objects
        TDunconstrained = [];
        TDconstrained = [];
        TDunconstrainedWashout = [];
        TD_tube_best = [];
        TD_unconst_tube = [];
        TD_unconst_tube_washout = [];
        TD_tt_int = [];
        TD_tt_rot = [];
        
        % Tube information
        startTargPos = [];
        constrainedTubeRadius = [];
        tubeObject = [];
        
        % Performance metrics
        block_size = [];
        performanceHit = [];
        performanceRecovery = [];
        successRate = [];
        expectedSuccessRate = [];
        expectedSuccessRate_washout
        
    end
    
    methods
        function EL_IT = IntTargExp(); end
    end
end