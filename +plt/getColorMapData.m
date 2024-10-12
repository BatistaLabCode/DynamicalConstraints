function cData = getColorMapData(x,cLim,map)
% getColorMapData       Get RGB color for desired colormap
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

if nargin > 2
    if ischar(map) % allow user to specify a colormap or provide one directly
        cMap = colormap(map);
    else
        cMap = map;
    end
else
    cMap = colormap; % Use default colormap if not specified
end

% Create data axis and find closest value
nPts = size(cMap,1);
dataMap = linspace(cLim(1),cLim(end),nPts);
[~,I] = min(abs(dataMap - x));

% Get color corresponding to index
cData = cMap(I,:);