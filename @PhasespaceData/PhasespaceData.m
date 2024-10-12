classdef PhasespaceData < KinematicData
    properties
        markerID = [];
        rawPositions = [];
        
        % Transformation parameters
        rotation = [];
        rotationCenter = [];
        scaling = [];
    end % End of properties
    methods
        function PD = PhasespaceData()
            % Class constructor
        end
    end % End of methods
end % End of classdef