classdef TubeObject < HSTData
   properties
       % Defined
       name = [];           % Name of tube
       path = [];           % Tube path
       radius = [];         % Tube radius
       displayRadius = [];  % Displayed radius
       rotation = [];       % Rotation
       tubeCode = [];       % Index of unique tube type (rotation/name)
       
       % Unused properties
       type = [];
       window = [];
       pseudorandomize = [];
       pseudorandomTarget = [];
       size = [];
       segment = 0;
       maxSegment = 0;
       
       % Private
       boundary = [];
   end  % End of properties
   methods
       function TU = TubeObject(); end    % Class constructor
       TU = calcBoundary(TU);   % Calculate tube boundary
       TU = plot(TU,varargin);  % Plot tube boundary
       TU = setTubeCode(TU);
       TU = setCondCode(TU,tC);
       containsTD = calcTubeContainment(TU,TD);
   end
end