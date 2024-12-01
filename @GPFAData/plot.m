function F = plot(GP,plotDim,varargin)
% Plot GPFA latent state vs time
%
% GP.plot(plotDim)
%
% Inputs:
%   GP             GPFA object
%   plotDim        Defined dimensions to plot
%
% Optional Arguments:
%   targSubset      Subset of targets to plot
%   nPlotMax        Maximum number of trials to plot per condition
%   condCode        Condition code (used to plot different conditions in different colors)
%   Tmin            Minimum timestep to plot
%   Tmax            Maximum timestep to plot
%   colorMode       Mode to use to select plot colors
%   alignMode       Mode used to align
%   xspec           Type of latent to plot ('xsm': smoothed, non-orthornormalized)
%   tStep           Time step for plotting (ms)
%   plotAlign       Plot line for alignment point
%   lineW           Width of trajectories
%   lineWavg        Width of average trajectories
%   color           Color place holder
%   ColMat          Colormap structure place holder
%   colIndMode      Method used to select colors based on provided condCode
%   colInfo         Color Info
%
% Adapted from 'plotEachDimVsTime' by Byron Yu
%
% Author:   Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

nPlotMax  = length(GP); % Maximum number of trials to plot per condition
condCode  = [];         % Condition code (used to plot different conditions in different colors)
Tmin      = [];         % Minimum timestep to plot
Tmax      = [];         % Maximum timestep to plot
colorMode = 'auto';     % Mode to use to select plot colors
alignMode = 'start';    % Mode used to align
xspec     = 'xsm';      % Type of latent to plot ('xsm': smoothed, non-orthornormalized)
tStep     = 500;        % Time step for plotting (ms)
plotAlign = true;       % Plot line for alignment point
lineW     = 0.5;        % Width of trajectories
lineWavg  = 2;          % Width of average trajectories
color     = [];         % Color place holder
ColMat    = [];         % Colormap structure place holder
colIndMode = 'normal';  % Method used to select colors based on provided condition code
colInfo     = [];       % Color Info

assignopts(who, varargin);

% Set colors for conditions.  If a condition code vector is supplied,
% choose colors that are evenly spaced in the HSV map.  Otherwise, plot
% in black.
if ~isempty(condCode)
    colorFlag = true;
    uniCode = unique(condCode);
    nCond = length(uniCode);

    switch colorMode
        case 'auto'
            if nCond > 2
                c = interpColor(nCond,1,'hsv',1);
            elseif nCond == 2
                c = [ones(1,3)*.75;1 0 0];
            else
                c = [ones(1,3)*.75];
            end
        case 'gridTask'
            c = getGridTaskCol;
        case 'manual'
            c = color;
        case 'ColMat'
            if GP(1).averaged
                c = ColMat.col_dark;
            else
                c = ColMat.col_light;
            end
    end
else
    colorFlag = false;
    c = zeros(1,3);
end

% Set line color saturation and line type.  If the data being plotted is
% averaged, plot thick dark lines.  If data are single-trial trajectories,
% use thin/light lines.
if GP(1).averaged
    lw = lineWavg;
    if ~strcmp(colorMode,'ColMat')
        cHSV = rgb2hsv(c);
        s = cHSV(:,end);
        s = 1 - s;
        s = s*1.5;
        s(s>1) = 1;
        cHSV(:,3) = 1-s;
        c = hsv2rgb(cHSV);
    end
else
    lw = lineW;
    if ~strcmp(colorMode,'ColMat')
        cHSV = rgb2hsv(c);
        cHSV(:,2) = cHSV(:,2) * 0.5;
        c = hsv2rgb(cHSV);
    end
end

% Get alignment time
switch alignMode
    case 'index'
        alignIdx = [GP.onsetIdx];
    case 'start'
        alignIdx = 1;
    case 'end'
        alignIdx = [GP.T];
end

% Get bin width
binWidth = unique([GP.binWidth]);
binWidth = binWidth(1); % Handle case where more than 1 bin width is present (this should never happen)

% Get all data and determine x-axis range
Xall = [GP.(xspec)];
xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1

% Loop over trials and plot
Tmin_all = nan(length(GP),1);
Tmax_all = nan(length(GP),1);
hold on;
for n = 1:min(length(GP), nPlotMax)
    % Get data for current trial
    dat = GP(n).(xspec);
    T   = GP(n).T;

    % Define time axis.  This allows the user to specify alignment points.
    t = 1:T;
    if ~isempty(alignIdx)
        t = t - alignIdx(n);
    end
    t = t*binWidth;
    Tmin_all(n) = t(1);
    Tmax_all(n) = t(end);

    if colorFlag
        % Choose appropriate color
        switch colIndMode
            case 'normal'
                % Select color as the index into the list of unique
                % conditions
                cInd = ismember(uniCode,condCode(n));
            case 'targetCode'
                % Condition code is the index into the provided ColMat
                % structure
                cInd = condCode(n);
        end
        col = c(cInd,:);
    else
        col = 0.2 * ones(1,3);
    end

    % If 'colInfo' is specified, use the provided color
    if ~isempty(colInfo)
        if GP(n).averaged
            col = colInfo{n,2};
        else
            col = colInfo{n,1};
        end
    end
    
    % Plot data for specified dimension
    plot(t, dat(plotDim,:), 'linewidth', lw, 'color', col);
end

% Set up axis limits
if isempty(Tmin)
  Tmin = min(Tmin_all);
end

if isempty(Tmax)
  Tmax = max(Tmax_all);
end

xtk     = [fliplr(0:-tStep:Tmin) tStep:tStep:Tmax];
ytk     = [-xMax 0 xMax];

axis([Tmin Tmax 1.1*min(ytk) 1.1*max(ytk)]);

% Plot dashed line for t=0
if plotAlign
    plot([0 0],get(gca,'YLim'),'k--')
end

set(gca,'xtick',xtk,'xticklabel',xtk);
set(gca,'ytick',ytk,'yticklabel',ytk);
set(gca,'TickDir','out')
xlabel('Time (ms)');

% Plot title
if isequal(xspec, 'xorth')
  str = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',plotDim);
else
  str = sprintf('$${\\mathbf x}_{%d,:}$$',plotDim);
end
title(str, 'interpreter', 'latex', 'fontsize', 12);