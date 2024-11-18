classdef KinematicData 
% Kinematic Data class
%
% This class contains time-varying kinematic data, including position,
% velocity, and acceleration data.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com
    properties
        % Kinematics
        time = [];  % Time
        pos = [];   % Position
        vel = [];   % Velocity
        acc = [];   % Acceleration
        
        % Data type/flags
        source = '';    % Source of kinematic data (typ. brain or hand)
        moveOnset = [];
        moveOffset = [];
    end % End of properties
    
    methods
        function K = KinematicData(); end   % Class constructor
        [k,binTimes] = binData(K,dataType,binWidth,t0,window,tOffset);
        [K] = calcMoveOnset(K,varargin);
        
        [s,t,K] = speedProfile(K,varargin); % Plot speed profile
    end % End of methods
end % End of class definition