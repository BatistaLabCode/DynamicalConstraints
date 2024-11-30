% Class definition for intermediate target experiment
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
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
        ttRotBlockDir = [];
        
        % Trajectory data objects
        TDunconstrained = [];
        TDconstrained = [];
        TD_tt_rot = [];
        
        % Tube information
        startTargPos = [];
        constrainedTubeRadius = [];
        tubeObject = [];
        
        % Performance metrics
        block_size = []; 
        successRate = [];
        predictedSuccessRate = [];
    end
    
    methods
        function EL_IT = IntTargExp(); end
    end
end