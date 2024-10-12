function trajectory(GP,varargin)
% Plot GPFA trajectory colored by target
%
% GP.plot(tC)
%
% Inputs:
%   GP             GPFA object
%   tC             Vector of target code values
%
% Optional Arguments:
%   targSubset      Subset of targets to plot
%   fracTrials      Fraction of trials to plot at random.
%
% Author:   Alan D. Degenhart
% Created:  2017.02.16
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional argument default values
targSubset = [];            % Subset of targets to plot
fracTrials = 1;             % Fraction of trials to plot at random.
condCode = [];              % Condition code (used for setting colors)
plotPreOnset = false;       % pre-onset data if onset index is specified
xSpec = 'xorth';            % Latent space
plotMode = '3D';            % plotting mode
plotDim = [1 2 3];          % Plot dimension
colorMode = 'auto';         % Determine the color code
col = [];                   % color place holder
colorIndMode = 'normal';    % 'normal': use unique target number for identifying colors
targetID = [];              % How to label/identify the conditions
lineStr = '.-';             % Default line type
markerSize = [];            % End point marker size
markerScale = 0.2;          % Scale for arrow marker
plotLegend = false;         % Don't plot legend by default
plotAxLabels = true;        % Plot axis labels by default
lineWidth     = 0.5;        % Width of trajectories
lineWidthavg  = 2;          % Width of average trajectories
preOnsetSamp = 0;           % Number of pre onset samples to plot
ColMat = [];                % Colormap structure place holder
colInfo = [];               % Color Info

% Parse optional agruments
assignopts(who,varargin);

% Set condition code if not specified
if isempty(condCode)
    condCode = ones(length(GP),1);
end

% Determine plot mode
if length(plotDim) == 2
    plotMode = '2D';
end

% Specify colors
[uniCond,I] = unique(condCode);
nTarg = length(uniCond);

% Get target names
if isempty(targetID)
    targNames = uniCond;
else
    targNames = targetID(I);
end

% Specify colormap based on number of targets.  If only 2 targets are
% provided, use red and black
switch colorMode
    case 'auto'
        if nTarg > 2
            c = interpColor(nTarg,1,'hsv',1);   % Dark
            cL = interpColor(nTarg,1,'hsv',1);  % Light
        elseif nTarg == 2
            c = [.5 .5 .5;1 0 0];
            cL = c;
        else
            c = [.5 .5 .5];
        end
    case 'gridTask'
        c = getGridTaskCol;
        cL = c;
    case 'twoTarget'
        c = [.8 .15 0;0 .6 .85];
        cL = c;
    case 'specified'
        c = col;
        cL = c;
    case 'ColMat'
        if GP(1).averaged
            c = ColMat.col_dark;
            cL = ColMat.col_dark;
        else
            c = ColMat.col_light;
            cL = ColMat.col_light;  
        end
end

% If subset is specified, get subset of trials and remove unnecessary
% colors as well
if ~isempty(targSubset)
    % Get subset of trials
    tCMask = ismember(condCode,targSubset);
    GP = GP(tCMask);
    condCode = condCode(tCMask);
    
    % Get subset of colors (this results in better color choices)
    uniMask = ismember(uniCond,targSubset);
    uniCond = uniCond(uniMask);
    targNames = targNames(uniMask);
    c = c(uniMask,:);
    cL = cL(uniMask,:);
end

% If there is only one target, plot in red
nTarg = length(uniCond);
if nTarg == 1 && ~strcmp(colorMode,'specified')
    c = [1 0 0];
    cL = c;
end

% Set plot line/color options depending on the data type.  This will allow
% the user to override these options using the optional input arguments.
% By default, individual trajectories will be plotted in a lighter color
% than trial-averaged ones

if ~GP(1).averaged
    % Individual trial trajectories
    lineWidth = lineWidth;
    if isempty(markerSize)
        markerSize = 5;
    end
    if ~strcmp(colorMode,'ColMat') && ~strcmp(colorMode,'specified')
        c_HSV = rgb2hsv(c);
        cL_HSV = rgb2hsv(cL);
        c_HSV(:,2) = c_HSV(:,2) * 0.5;
        cL_HSV(:,2) = cL_HSV(:,2) * 0.5;
        c = hsv2rgb(c_HSV);
        cL = hsv2rgb(cL_HSV);
    end
else
    lineWidth = lineWidthavg;
    if isempty(markerSize)
        markerSize = 10;
    end
end

hold on
% Loop over targets, then trials to plot
AX = nan(nTarg,1);
legStr = cell(nTarg,1);
for i = 1:nTarg
    % Get data for specific target
    tcMask = (condCode == uniCond(i));
    GPtemp = GP(tcMask);
    nSeq = sum(tcMask);
    
    if ~isempty(colInfo)
        colInfoTemp = colInfo(tcMask,:);
    end
    
    % Determine color index mode
    switch colorIndMode
        case 'normal'
            cInd = i;
        case 'targetCode'
            cInd = uniCond(i);
    end
    
    % Handle case where only one target exists.  When this is the case,
    % only one color is used
    if nTarg == 1 && ~strcmp(colorMode,'specified')
        cInd = 1;
    end
        
    % Downsampling trajectories if desired
    if fracTrials < 1
        nFrac = round(nSeq*fracTrials); % Number of trials to plot
        trialNum = randperm(nSeq,nFrac);
        GPtemp = GPtemp(trialNum);
        nSeq = nFrac;
        
        if ~isempty(colInfo)
            colInfoTemp = colInfoTemp(trialNum,:);
        end
    end
    
    % Plot each trajectory
    for j = 1:nSeq
        dat = GPtemp(j).(xSpec)(plotDim,:);
        
        % Plot pre-onset data if onset index is specified
        idx = GPtemp(j).onsetIdx;
        if plotPreOnset
            switch plotMode
                case '3D'
                    plot3(dat(1,1:(idx-1)),dat(2,1:(idx-1)), ...
                        dat(3,1:(idx-1)),lineStr,'color',ones(1,3)*.25);
                    plot3(dat(1,1),dat(2,1),dat(3,1),'o','color','k',...
                        'MarkerFaceColor','k','MarkerSize',5)      % Start (black)
                case '2D'
                    plot(dat(1,1:(idx-1)),dat(2,1:(idx-1)), ...
                        lineStr,'color',ones(1,3)*.25);
                    plot(dat(1,1),dat(2,1),'o','color','k',...
                        'MarkerFaceColor','k','MarkerSize',5)      % Start (black)
            end
        end
        
        % Adjust onset index if pre-onset samples are to be plotted
        idx = idx - preOnsetSamp;
        
        % Get color for plotting
        if ~isempty(colInfo)
            if GPtemp(j).averaged
                colTemp = colInfoTemp{j,2};
            else
                colTemp = colInfoTemp{j,1};
            end
        else
            colTemp = c(cInd,:);
        end
        
        % Plot trajectory
        switch plotMode
            case '3D'
                AX(i) = plot3(dat(1,idx:end),dat(2,idx:end),dat(3,idx:end), ...
                    lineStr,'color',cL(cInd,:),'LineWidth',lineWidth,'MarkerSize',markerSize/4);
                % Plot start and end points
                plot3(dat(1,idx),dat(2,idx),dat(3,idx),'o','color','k',...
                    'MarkerFaceColor',colTemp,'MarkerSize',markerSize/2)      % Start (black)
                plot3(dat(1,end),dat(2,end),dat(3,end),'o','color','k',...
                    'MarkerFaceColor',colTemp,'MarkerSize',markerSize)   % End (colored)
            case '2D'
                % Plot trajectory
                AX(i) = plot(dat(1,idx:end),dat(2,idx:end), ...
                    lineStr,'color',colTemp,'LineWidth',lineWidth);
                
                % Plot start and end points
                plot(dat(1,idx),dat(2,idx),'o', ...
                    'color','k', ...
                    'MarkerFaceColor', colTemp, ...
                    'MarkerEdgeColor', 'None', ...
                    'MarkerSize', markerSize/2)
                plt.drawArrowMarker(dat(1,:),dat(2,:),colTemp,markerScale);
        end
    end
    if ~isempty(colInfo)
        legStr{i} = colInfoTemp{1,3};
    elseif ~strcmp(colorMode,'ColMat')
        legStr{i} = sprintf('Target %d',targNames(i));
    else
        legStr{i} = ColMat.condLabel{uniCond(i)};
    end
end

if plotLegend
    Hl = legend(AX,legStr,'Location','SouthWest','AutoUpdate','off');
    set(Hl,'box','off')
end

set(gca,'TickDir','out')

% Generate axis labels
if plotAxLabels
    %axis equal; % This can make scaling problematic
    switch xSpec
        case 'xsm'
            str1 = sprintf('$${\\mathbf x}_{%d,:}$$',plotDim(1));
            str2 = sprintf('$${\\mathbf x}_{%d,:}$$',plotDim(2));
        case 'xorth'
            str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',plotDim(1));
            str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',plotDim(2));
    end

    xlabel(str1,'fontsize', 16);
    ylabel(str2,'fontsize', 16);

    if strcmp(plotMode,'3D')
        switch xSpec
            case 'xsm'
                str3 = sprintf('$${\\mathbf x}_{%d,:}$$',plotDim(3));
            case 'xorth'
                str3 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',plotDim(3));
        end
        zlabel(str3,'fontsize',16);
    end
else
    %set(gca,'XTickLabel',[],'YTickLabel',[])
end