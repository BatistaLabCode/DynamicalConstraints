% drawArrowMarker       Draw arrow-style marker at the endpoint of a line.
%
% This function draws an arrow-style marker at the end of a line.  The
% direction of the marker is along the direction defined by the last two
% data points of the provided data.  The 
%
% Usage:
%   plt.drawArrowMarker(x,y,c,sz)
%
% Inputs:
%   x   x-axis data
%   y   y-axis data
%   c   Color
%   sz  Desired size (length of the arrow)
%
% Optional Inputs:
%   edge_color   Set color
%
% Outputs:
%   H   Handle to the drawn marker
%
% Author:   Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function H = drawArrowMarker(x,y,c,sz,varargin)

% Optional arguments
edge_color = 'k';

assignopts(who, varargin);

% Truncate x and y to make sure they only have 2 elements
x = x(end-1:end);
y = y(end-1:end);

% Get vector
xc = x-x(1);
yc = y-y(1);
d = norm([xc(2) yc(2)]);
xc = (xc/d) * sz;
yc = (yc/d) * sz;

% Determine angle and rotation matrix
a = atan2(yc(end),xc(end));
R = [cos(a) -sin(a);sin(a) cos(a)];

% Define unscaled marker positions
xLim = [-0.25 1];
yLim = 0.5;
xM = [xLim(2) xLim(1) 0 xLim(1) xLim(2)];
yM = [0 yLim 0 -yLim 0];

P = R*[xM;yM];

% Determine offset
xc = xc-xc(2);
yc = yc-yc(2);
o = [x(2) + xc(1);y(2) + yc(1)];
O = repmat(o,1,length(xM));
Pscaled = P*sz+O;

% Draw marker
H = patch(Pscaled(1,:),Pscaled(2,:),c);
set(H,'EdgeColor', edge_color)

% Set rotation and offset data (this could be used in the future to
% re-scale the marker size if the plot is re-drawn, but currently is not
% used)
Data.P = P;
Data.O = O;
set(H,'UserData',Data);