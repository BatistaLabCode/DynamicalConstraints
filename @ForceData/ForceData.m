classdef ForceData < KinematicData
    % Force data class
    %
    % This class inherits from the 'KinematicData' class.  Force
    % technically isn't a type of kinematic data, but is left here for
    % consistency (instead of changing the name of the parent class)
    %
    % Copyright (C) by Alan Degenhart and Erinn Grigsby
    % Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com
    properties
        % In addition to the parent class properties (time, pos, etc.), we
        % need to keep track of the actual force and torque.  There is a
        % bit of ambiguity here, as this code doesn't really distinguish
        % between force/kinematic data that is used for control or not
        force = [];
        torque = [];
        
    end
    
    methods
        function F = ForceData()
            F.source = 'force';
        end
        
        % Note -- might need to overload the movement onset/speed profile
        % code for this class at some point
    end
end