function plotTitle(titleStr,subtitleStr,varargin)
% Generate title for multi-plot figure
%
% Inputs:
%  titleStr         Main figure title.
%  subtitleStr      Sub-title.
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

interpreter = 'none';
FontSize = 14;
assignopts(who,varargin);

% Determine plot size.  Do this based on the figure size.  This ensures
% that the location of the title text on the figure is the same regardless
% of the figure size.
sz = get(gcf, 'Position');
h = 34 / sz(4);

axes('Position',[0 1-h 1 h])
% Plot title
A = text(0.5,1,titleStr, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontWeight','normal','Interpreter',interpreter,'FontSize',FontSize);
% Plot sub-title if desired
if (nargin == 2) && (~isempty(subtitleStr))
    B = text(0.5,0.5,subtitleStr, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'Interpreter',interpreter, ...
        'FontSize', FontSize-2); % Make subtitle font size slightly smaller
end
axis off