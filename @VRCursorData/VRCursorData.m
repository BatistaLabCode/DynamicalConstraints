classdef VRCursorData < KinematicData
    properties
        markerID = [];
        
        % Transformation parameters
        rotation = [];
        rotationCenter = [];
        scaling = [];
    end % End of properties
    methods
        function CD = VRCursorData()
            % Class constructor
        end
    end % End of methods
end % End of classdef