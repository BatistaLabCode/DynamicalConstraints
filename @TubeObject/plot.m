function plot(TU,varargin)
% Plot method for Tube class

color = 'k';
lineSpec = '-';
lineWidth = 1;
tC = [];
plotSegment = [];

assignopts(who,varargin);

% Unique tubes should be defined for each target.  

% Get names and rotations for unique 

% Determine if we are just plotting a segment of the tube.
if ~isempty(plotSegment)
    pathIdx = round(linspace(1,size(TU.path,1),TU.maxSegment+1));
    TU.path = TU.path(pathIdx(plotSegment):pathIdx(plotSegment+1),:);
    TU = TU.calcBoundary;
end

% Get boundary and plot
x = TU.boundary(:,1);
y = TU.boundary(:,2);
plot(x,y,lineSpec,'color',color,'LineWidth',lineWidth)
axis equal