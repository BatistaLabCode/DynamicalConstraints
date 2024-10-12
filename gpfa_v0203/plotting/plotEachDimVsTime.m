function f = plotEachDimVsTime(seq, xspec, binWidth, varargin)
%
% plotEachDimVsTime(seq, xspec, binWidth, ...)
%
% Plot each state dimension versus time in a separate panel.
%
% INPUTS:
%
% seq       - data structure containing extracted trajectories
% xspec     - field name of trajectories in 'seq' to be plotted 
%             (e.g., 'xorth' or 'xsm')
% binWidth  - spike bin width used when fitting model
%
% OPTIONAL ARGUMENTS:
%
% nPlotMax  - maximum number of trials to plot (default: 20)
% redTrials - vector of trialIds whose trajectories are plotted in red
%             (default: [])
% nCols     - number of subplot columns (default: 4)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

nPlotMax  = 20;
redTrials = [];
nCols     = 4;
condCode  = [];   % Condition code (used to plot different conditions in different colors)
alignTime = [];   % Timepoint for each trial to align (e.g., movement onset)
Tmin      = [];   % Minimum timestep to plot
Tmax      = [];   % Maximum timestep to plot
colorMode = 'auto'; % Mode to use to select plot colors

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
                c = [0 0 0;1 0 0];
            else
                c = [0 0 0];
            end
        case 'gridTask'
            c = getGridTaskCol;
    end
else
    colorFlag = false;
    c = zeros(1,3);
end

f = figure;
pos = get(gcf, 'position');
set(f, 'position', [100 100 1500 800]);

Xall = [seq.(xspec)];
xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1

%   Tmax    = max([seq.T]);  
%   xtkStep = ceil(Tmax/25)*5;
%   xtk     = 1:xtkStep:Tmax;
%   xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;
%   ytk     = [-xMax 0 xMax];

nRows   = ceil(size(Xall, 1) / nCols);

Tmin_all = nan(length(seq),1);
Tmax_all = nan(length(seq),1);

for n = 1:min(length(seq), nPlotMax)
dat = seq(n).(xspec);
T   = seq(n).T;

% Define time axis.  This allows the user to specify alignment points.
t = 1:T;
if ~isempty(alignTime)
    tAlign = round(alignTime(n)/binWidth);
    t = t - tAlign;
end
Tmin_all(n) = t(1);
Tmax_all(n) = t(end);

if colorFlag
    % Choose appropriate color
    cMask = ismember(uniCode,condCode(n));
    col = c(cMask,:);
else
    col = 0.2 * ones(1,3);
end

for k = 1:size(dat,1)
  subplot(nRows, nCols, k);
  hold on;
  lw = 0.05;
%       if ismember(seq(n).trialId, redTrials)
%         col = [1 0 0]; % red
%         lw  = 3;
%       else
%         col = 0.2 * [1 1 1]; % gray
%         lw = 0.05;
%       end      
  plot(t, dat(k,:), 'linewidth', lw, 'color', col);
end
end

% Set up axis limits
if isempty(Tmin)
  Tmin = min(Tmin_all);
else
  Tmin = Tmin/binWidth;
end

if isempty(Tmax)
  Tmax = max(Tmax_all);
else
  Tmax = Tmax/binWidth;
end

xtkStep = ceil((Tmax-Tmin)/25)*5;
xtk     = Tmin:xtkStep:Tmax;
xtkl    = (Tmin*binWidth):(xtkStep*binWidth):(Tmax-1)*binWidth;
ytk     = [-xMax 0 xMax];

for k = 1:size(dat,1)
h = subplot(nRows, nCols, k);
axis([Tmin Tmax 1.1*min(ytk) 1.1*max(ytk)]);

% Plot dashed line for t=0
hold on
plot([0 0],[-xMax xMax],'k--')

if isequal(xspec, 'xorth')
  str = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',k);
else
  str = sprintf('$${\\mathbf x}_{%d,:}$$',k);
end
title(str, 'interpreter', 'latex', 'fontsize', 16);

set(h, 'xtick', xtk, 'xticklabel', xtkl);
set(h, 'ytick', ytk, 'yticklabel', ytk);
xlabel('Time (ms)');
end
