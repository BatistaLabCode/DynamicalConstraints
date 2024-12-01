function c = interpColor(n,sF,map,wrapOption,varargin)
% [c] = interpColor(n,sF,map,wrapOption)
%
% Return interpolated colors for multiple objects.
%
% This function returns a color array of N colors, such that the colors are
% distributed evenly over the specified colomap.  This function is intended
% for conditions where object colors are intended to vary gradually.
%
% Inputs:
%   n              Number of points.
%   sF             Color scale factor
%   map            Defined colormap
%   wrapOption     Define how to interpolate the data.
%   
% Optional Inputs:
%   valueScale     Sets the value range.
% 
% Outputs:
%   c               Colormap matrix
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional arguments
valueScale = 1;         % Sets the value range

assignopts(who,varargin);

% Get desired colormap.  As the 'colormap' function needs an active figure,
% one is created, the map is set/returned, and the figure is closed.
F = figure;
m = colormap(map);
close(F)
nCol = size(m,1);

% Define positions to evaluate interpolation at.
if wrapOption
    nTemp = n + 1;
else
    nTemp = n;
end

% Determine points to evaluate interpolation at
xi = linspace(1,nCol,nTemp);
x = 1:nCol;

% Interpolate to get new colormap
c = nan(nTemp,3);
for i = 1:3
    c(:,i) = interp1(x,m(:,i),xi,'spline');
end

% Get rid of extra entry (used for 'wrapOption').  Also set range to [0 1],
% as some values may be slightly below 0 or above 1 due to interpolation.
c = c(1:n,:);
c(c<0) = 0;
c(c>1) = 1;

% Convert to HSV
cHSV = rgb2hsv(c);

% Set scale factor
cHSV(:,2) = cHSV(:,2) * sF;
cHSV(:,3) = cHSV(:,3) * valueScale;

% Convert back to RGB
c = hsv2rgb(cHSV);