% quiver_new   Create quiver plot (with added functionality)
%
% Usage:
%    plt.quiver_new(x, y, dx, dy)
%
% Inputs:
%   x   Vector of x-positions
%   y   Vector of y-positions
%   dx  Vector of x-axis velocities (1 per position)
%   dy  Vector of y-axis velocities (1 per position)
%
% Optional inputs:
%   scale_factor   Scale factor to apply to vectors (default: 1)
%   color          Color for vectors (one, or one per data point)
%
% Outputs:
%   F   Figure handle
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function quiver_new(x, y, dx, dy, varargin)

% Parse optional arguments
scale_factor = 1;
color = [0, 0, 0];
line_width = 1;
arrow_size = 7;

assignopts(who, varargin);
hold on;

% Make sure the 'color' input is the appropriate size
n_pts = length(x);
if size(color, 1) == 1
    color = repmat(color, n_pts, 1);
end

% Loop over points and plot
for i = 1:n_pts
    % Define x and y vectors to plot
    x_temp = [x(i), x(i) + dx(i)*scale_factor];
    y_temp = [y(i), y(i) + dy(i)*scale_factor];
    
    % Plot line
    plot(x_temp, y_temp, 'color', color(i, :), 'LineWidth', line_width)
    
    % Plot marker -- don't show the edges
    if ~isempty(arrow_size)
        plt.drawArrowMarker(x_temp, y_temp, color(i, :), arrow_size, ...
            'edge_color', 'none');
    end
end