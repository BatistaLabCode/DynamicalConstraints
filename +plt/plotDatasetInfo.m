function plotDatasetInfo(subject,date,trialID,varargin)
% Plot dataset info on current figure
%
% Useage:
%   plotDatasetInfo(subject,date,trialID)
%
% Inputs:
%   subject    Subject ID
%   date       Session date
%   trialId    Trial number
%
% Optional Inputs:
%   plotMode        Determine if the information is plotted on a figure or
%                       on a axis.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

plotMode = 'figure';

assignopts(who,varargin);

% Generate dataset string
trialRange = [min(trialID) max(trialID)];
s = sprintf('%s%strl%d:%d',subject,datestr(date,'yyyymmdd'), ...
    trialRange(1),trialRange(2));

switch plotMode
    case 'figure'
        % Generate plot axes
        axes('Position',[.5 0 .5 .05]);

        % Plot dataset string
        h = text(.95,.5,s,'HorizontalAlignment','right');
        axis off
    case 'axis'
        % Get axis limits
        xLim = get(gca,'XLim');
        yLim = get(gca,'YLim');
        
        h = text(xLim(1),yLim(1),s,'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom');
end