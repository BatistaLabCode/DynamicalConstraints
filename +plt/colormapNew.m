function [M] = colormapNew(Ax,mapName)
% [M] = colormapNew(mapName)
%
% colormap function with updated functionality.
%
% This function behaves in the same way as the standard 'colormap'
% function, but allows for user-defined custom colors.
%
% Author:       Alan D. Degenhart
% Created:      2016.11.21
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Handle case where only the map name is specified.  In this case, there is
% only one input
if nargin == 1
    if ischar(Ax)
        mapName = Ax;
        Ax = gca;
    else
        error('Must specify color map name as a string.')
    end
end

mapName = lower(mapName);

switch mapName
    case 'blue-orange'
        O = [linspace(0.2,1,32)' linspace(.1,.5,32)' linspace(0,0,32)'];
        B = [linspace(0,0,32)' linspace(.1,.5,32)' linspace(.2,1,32)'];
        M = [flipud(B);O];
        colormap(Ax,M)
    case 'red'
        R = linspace(1,0,64)';
        M = [ones(64,1) R R];
        colormap(Ax,M)
    case 'green'
        cRGB = [0 .75 0]; % A darker green
        cHSV = rgb2hsv(cRGB);
        mapHSV = repmat(cHSV,64,1);
        S = mapHSV(1,2);
        V = mapHSV(1,3);
        %S = linspace(S,0,64);
        V = linspace(V,0,64);
        mapHSV(:,2) = S;
        mapHSV(:,3) = V;
        M = flipud(hsv2rgb(mapHSV));
        colormap(Ax,M)
    case 'red-black'
        cRGB = [0.75 0 0]; % A darker green
        cHSV = rgb2hsv(cRGB);
        mapHSV = repmat(cHSV,64,1);
        S = mapHSV(1,2);
        V = mapHSV(1,3);
        %S = linspace(S,0,64);
        V = linspace(V,0,64);
        mapHSV(:,2) = S;
        mapHSV(:,3) = V;
        M = flipud(hsv2rgb(mapHSV));
        colormap(Ax,M)
    case 'blue-red'
        M = [linspace(1,0,64)' zeros(64,1) linspace(0,1,64)'];
        M = flipud(M);
        colormap(Ax,M)
    otherwise
        colormap(Ax,mapName)
end