function plot(TD,varargin)
% plot          Plotting function for trajectoryData class
%
% Usage:
%   TD.plot()
%
% There are numerous optional arguments to adjust the plot including color
% settings, subselecting trials, trajectory parameters, target parameters,
% etc. See full list of options below, grouped by function.
%
% Optional Inputs:
%   1) Define trial and timings to plot:
%       trialsPerTarg     Allows for subselection of trials
%       setRandSeed       Allows user to fix the random permutation for trial 
%                            subselection. Default: randomization not fixed.
%       plotStates        Plot only the desired states.
%       plotSubSetTimePt  Plot only the desired number of time points.
%       flipSubSet        Determine whether the desired number of time
%                            points is taking from the start or end of the
%                            trajectory.
%   2) Defines colors:
%       colorMap          Default color map to use
%       color             Specify exact color, this will overwrite the default
%       nCol              Number of colors to use
%       ColMat            Use a color map structure with target definitions
%       trajColorWeight   Saturation of the trajectory color
%       targetFaceCol     Define target face color, overwrites the default
%       sF                Scale factor for target colors
%       intTargCol        Set intermediate target color. Default-Light gray
%       intTargEdgeCol    Set intermediate target edge color.
%       targColMatIdx     Determine whether to plot the targets light (1) 
%                             or dark (2)
%   3) Define trajectory options:
%       plotTrajectories  Determines whether to plot the trajectories
%       lineWidth         Width of trajectories
%       lineStyle         Line type default '-'
%       markerSize        Size of the display marker
%       startMarker       Define the marker for the start of the trajectory.
%       endMarker         Define the marker for the end of the trajectory. 
%                           This includes an arrow though it will just be a
%                           visualization and not an exact match to arrows
%                           defined with illustrator.
%   4) Define target options:
%       plotTargets       Determines whether to display the targets
%       plotStartTarg     Plot start target
%       plotEndTarg       Plot end target
%       plotIntTarg       Plot targets for intermediate states
%       targetScale       Scale factor for the target size
%       targetStyle       Target style (circle or square)
%       targetLineColor   Target outline color
%       targetLineStyle   Define target line style
%       intTargWidth      Set the intermediate target line width
%       plotTargNum       Determine whether to display the target number
%       fontSize          Target Font Size
%       uniTargNum        Determine whether to plot the unique target 
%                           number or the unique condition number.
%   5) Define plotting options:
%       verbose           Display messages
%       plot3d            Determine if the trajectories should plot in 3D
%       axesLimits        Axis limits
%
% Note: do not need to worry about the cursor radius here, as the actual
% target window requirements are such that if the *center* of the cursor is
% inside of the window, the window requirements have been met.  We
% arbitrarily set the size of the cursor and target *objects* to reflect
% the desired feedback to the animal.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

%% Set default values for optional arguments

% Define trial and timings to plot
trialsPerTarg = [];     % Allows for subselection of trials
setRandSeed = [];       % Allows user to fix the random permutation for trial 
                        %   subselection. Default: randomization not fixed.
plotStates = [];        % Plot only the desired states.
plotSubSetTimePt = [];  % Plot only the desired number of time points.
flipSubSet = false;     % Determine whether the desired number of time
                        %   points is taking from the start or end of the
                        %   trajectory.

% Defines colors
colorMap = 'hsv';       % Default color map to use
color = [];             % Specify exact color, this will overwrite the default
nCol = [];              % Number of colors to use
ColMat = [];            % Use a color map structure with target definitions
trajColorWeight = 1;    % Saturation of the trajectory color
targetFaceCol = [];     % Define target face color, overwrites the default
sF = 0.75;              % Scale factor for target colors
intTargCol = .9*[1 1 1];% Set intermediate target color. Default: light gray
intTargEdgeCol = .5*[1 1 1]; % Set intermediate target edge color.
targColMatIdx = 1;      % Determine whether to plot the targets light (1) 
                        %   or dark (2)

% Define trajectory options
plotTrajectories = 1;   % Determines whether to plot the trajectories
lineWidth = [];         % Width of trajectories
lineStyle = [];         % Line type default '-'
markerSize = [];        % Size of the display marker
startMarker = 'o';      % Define the marker for the start of the trajectory.
endMarker = 'o';        % Define the marker for the end of the trajectory. 
                        %   This includes an arrow though it will just be a
                        %   visualization and not an exact match to arrows
                        %   defined with illustrator.

% Define target options
plotTargets = 1;        % Determines whether to display the targets
plotStartTarg = false;  % Plot start target
plotEndTarg = true;     % Plot end target
plotIntTarg = false;    % Plot targets for intermediate states
targetScale = 1;        % Scale factor for the target size
targetStyle = 'circle'; % Target style (circle or square)
targetLineColor = 'k';  % Target outline color
targetLineStyle = '-';  % Define target line style
intTargWidth = 1.5;     % Set the intermediate target line width
plotTargNum = true;     % Determine whether to display the target number
fontSize = 24;          % Target Font Size
uniTargNum = false;     % Determine whether to plot the unique target 
                        %   number or the unique condition number.

% Define plotting options
verbose = true;         % Display messages
plot3d = false;         % Determine if the trajectories should plot in 3D
axesLimits = [];        % Axis limits


%% Parse optional agruments
assignopts(who,varargin);

nTrials = length(TD);
hold on

% Set line width.  If this is specified, use the provided value. If not,
% choose depending on whether the data is averaged
if isempty(lineWidth)
    if TD(1).averaged
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end
end

% Plot targets
tC = [TD.targetCode];
[uniTC,IA] = unique(tC);
nTarg = length(uniTC);

% Get target codes for start and end targets
startPos = [TD.startPos]';
endPos = [TD.targPos]';
[tcStart] = util.pos2targetCode([startPos; endPos]);
tcEnd = tcStart(size(startPos,1)+1:end);
tcStart = tcStart(1:size(startPos,1));

% Get unique target combinations
tcComb = [tcStart tcEnd];

% Allow for user to speficy the number of colors to interpolate over.  This
% handles the case where fewer than the full number of unique targets in a
% larger dataset are present in the dataset being plotted, and ensures that
% target and trajectory colors are consistent across mulitple plots if
% desired.
if isempty(nCol)
    nCol = nTarg;
    colIndMode = 'targNumber'; % Use target number as index into color matrix
else
    % If specifying the number of plot colors to use, check to make sure
    % that the number of colors specified is equal to or greater than the
    % maximum target code value.
    assert(nCol >= max(tC),'Number of unique targets exceeds the number of plot colors specified.')
    colIndMode = 'targetCode'; % Use target code as index into color matrix
end

% Handle case where a single color is provided.
if ~iscell(color) && size(color,1) == 1
    color = repmat(color,nTarg,1);
end

% Get target colors
if isempty(color)
    cTarg = plt.interpColor(nCol,sF,colorMap,1);
elseif iscell(color)
    cTarg = color;
else
    % Convert to HSV
    cHSV = rgb2hsv(color);
    % Set scale factor
    cHSV(:,2) = cHSV(:,2) * sF;
    % Convert back to RGB
    cTarg = hsv2rgb(cHSV);
end

% Define target style.  MonkeyHost supports square of circular targets.
% This information should be pulled from the data file, but for now it has
% been added as an optional argument.
switch targetStyle
    case 'circle'
        curv = [1 1];
    case 'square'
        curv = [0 0];
end

if plotTargets
    for i = 1:nTarg
        % Get start and end target positions
        startPos = TD(IA(i)).startPos; % Center of target
        targPos = TD(IA(i)).targPos; % Center of target

        % If colors are automatically provided, then it is assumed that the
        % plot colors correspond to the unique target values.  Future
        % updates to this function could ensure that custom color maps
        % could be provided which specify the color code as target-code
        % dependent.
        switch colIndMode
            case 'targNumber'
                colInd = i;
            case 'targetCode'
                colInd = uniTC(i);
        end

        % Get color
        if iscell(color)
            cTemp = cTarg{colInd};
        else
            cTemp = cTarg(colInd,:);
        end

        % If ColMat is not empty, use it to determine the target color
        % automatically
        if ~isempty(ColMat)
            col_info = util.getColorInfo(ColMat, startPos',targPos');
            cTemp = col_info{targColMatIdx};
        end

        % If target face color was specified explicitly, override
        if ~isempty(targetFaceCol)
            cTemp = targetFaceCol;
        end

        % Determine the target text
        if uniTargNum
            targText = num2str(tcComb(IA(i),2));
        else
            targText = num2str(uniTC(i));
        end

        % Plot end target position if desired
        if plotEndTarg
            % Plot target
            r = targetScale*TD(IA(i)).targSize;      % Target radius
            targPosCorner = targPos - r;
            rectangle('Position',[targPosCorner(1),targPosCorner(2),2*r,2*r], ...
                'Curvature',curv,'FaceColor',cTemp,'EdgeColor',targetLineColor, ...
                'LineStyle',targetLineStyle)

            % Plot target number text
            if plotTargNum
                text(targPos(1),targPos(2),targText, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'FontSize',fontSize,'FontWeight','bold','color',targetLineColor)
            end
        end

        % Plot starting position if desired.  Note that currently the
        % trajectory data structure does not have a field for the start
        % position/target size.  Usually this will be the target position,
        % but eventually the class should be updated to contain this
        % information.
        if plotStartTarg
            r = targetScale*TD(IA(i)).targSize;      % Target radius
            startPosCorner = startPos - r;
            rectangle('Position',[startPosCorner(1),startPosCorner(2),2*r,2*r], ...
                'Curvature',curv,'FaceColor',cTemp,'EdgeColor',targetLineColor)
        end
    end
end

% Plot intermediate targets
if plotIntTarg && plotTargets
    % Get intermediate target position and size
    intTargPos = [TD.intTargPos];
    intTargSz = [TD.intTargSz];
    intTargInfo = [intTargPos' intTargSz'];
    uniIntTarg = unique(intTargInfo,'rows');

    % Get end target position.  Any intermediate target positions that are
    % the same as the end target positions will not be plotted.
    endTargPos = [TD.targPos]';
    uniEndTargPos = unique(endTargPos,'rows');

    mask = ismember(uniIntTarg(:,1:end-1),uniEndTargPos,'rows');
    uniIntTarg = uniIntTarg(~mask,:);

    % Loop over intermediate targets and plot
    for i = 1:size(uniIntTarg,1)
        targPos = uniIntTarg(i,1:end-1)';   % Center of target
        r = uniIntTarg(i,end);              % Target radius
        targPosCorner = targPos - r;
        rectangle('Position',[targPosCorner(1),targPosCorner(2),2*r,2*r], ...
            'Curvature',curv,'FaceColor',intTargCol,'EdgeColor',intTargEdgeCol,...
            'LineWidth',intTargWidth)
    end
end

% Get trajectory colors
if isempty(color)
    cTraj = plt.interpColor(nCol,trajColorWeight,colorMap,1);
elseif iscell(color)
    cTraj = color;
else
    % Convert to HSV
    cHSV = rgb2hsv(color);
    % Set scale factor
    cHSV(:,2) = cHSV(:,2) * trajColorWeight;
    % Convert back to RGB
    cTraj = hsv2rgb(cHSV);
end

% Subselect trajectories if desired
if ~isempty(trialsPerTarg)
    trialSubset = cell(1,nTarg);
    for i = 1:nTarg
        targMask = tC == uniTC(i);
        trl = find(targMask);
        if ~isempty(setRandSeed)
            rng(setRandSeed)
            trl = trl(randperm(length(trl)));
        else
            trl = trl(randperm(length(trl)));
        end
        trialSubset{i} = trl(1:trialsPerTarg);
    end

    trlList = sort(cell2mat(trialSubset));
    TD = TD(trlList);
    nTrials = length(TD);
end


if plotTrajectories
    % Plot individual trajectories
    for i = 1:nTrials
        tC = TD(i).targetCode;
        % Get start and end target positions
        startPos = TD(i).startPos; % Center of target
        targPos = TD(i).targPos; % Center of target

        % Determine marker
        if isempty(markerSize)
            if TD(i).averaged
                markerSize = 10;
            else
                markerSize = 5;
            end
        end

        % Get color index
        switch colIndMode
            case 'targNumber'
                colInd = find(uniTC == tC);
            case 'targetCode'
                colInd = tC;
        end

        % Get color
        if iscell(color)
            cTemp = cTraj{colInd};
        else
            cTemp = cTraj(colInd,:);
        end

        % Get trajectory color if 'ColMat' is specified.  Currently only
        % set the color based on the end target.
        if ~isempty(ColMat)
            col_info = util.getColorInfo(ColMat, startPos',targPos');
            cTemp = col_info{2};

            % Convert to HSV
            cHSV = rgb2hsv(cTemp);
            % Set scale factor
            cHSV(:,2) = cHSV(:,2) * trajColorWeight;
            % Convert back to RGB
            cTemp = hsv2rgb(cHSV);
        end

        % Get position data.  
        srcStr = 'brainKin';
        pos = TD(i).(srcStr).pos;

        % Select line style based on trial success if line Style is empty
        if isempty(lineStyle)
            if logical(TD(i).successful)
                lineStyleTemp = '-';
            else
                lineStyleTemp = '--';
            end
        else
            lineStyleTemp = lineStyle;
        end

        % Check if the trajectory is being truncated to specific states
        if ~isempty(plotStates)
            % Find the desired states
            states = find(ismember(TD(i).states,plotStates));

            if ~isempty(states)
                % Find the times of the desired states
                stateOn = TD(i).stateOnset(states);
                stateOff = TD(i).stateOffset(states);
                trunTime = [min(stateOn) max(stateOff)];

                % Get the timing data within the truncation time period and
                % shorten the position data.
                timeKin = TD(i).(srcStr).time;
                trunPts = timeKin >= trunTime(1) & timeKin < trunTime(2);
                pos = pos(trunPts,:);
            else
                if verbose
                    fprintf(['\n Desired state not found for trial '...
                        num2str(TD(i).trialID)])
                end
                pos = [];
            end
        end

        % Check if the trajectory is being truncated
        if ~isempty(plotSubSetTimePt)
            if size(pos,1)>=plotSubSetTimePt(end)
                if flipSubSet
                    pos = pos(end-(fliplr(plotSubSetTimePt)-1),:);
                else
                    pos = pos(plotSubSetTimePt,:);
                end
            elseif size(pos,1)>=plotSubSetTimePt(1)
                if verbose
                    fprintf(['\n Position data for Trial ' num2str(TD(i).trialID)...
                        ' is shorter than expected. Plotting available data'])
                end
                if flipSubSet
                    idx = find(size(pos,1)-(plotSubSetTimePt-1)>0);
                    pos = pos(end - (fliplr(plotSubSetTimePt(idx))-1),:);
                else
                    pos = pos(plotSubSetTimePt(1):end,:);
                end
            else
                if verbose
                    fprintf(['\n No position data in the desired range for '...
                        'Trial ' num2str(TD(i).trialID) ' Skipping'])
                end
                continue
            end
        end

        % Plot the trajectory
        if size(pos,1) > 0
            if plot3d
                % Plot entire trajectory
                plot3(pos(:,1),pos(:,2),pos(:,3),'color',cTemp,...
                    'LineWidth',lineWidth,'LineStyle',lineStyleTemp)
                % Plot start position
                plot3(pos(1,1),pos(1,2),pos(1,3),startMarker,'MarkerSize',...
                    markerSize/2,'MarkerFaceColor',cTemp, 'MarkerEdgeColor','None')
                % Plot end position
                plot3(pos(end,1),pos(end,2),pos(end,3),endMarker,...
                    'MarkerSize',markerSize,'MarkerFaceColor',cTemp,...
                    'MarkerEdgeColor', 'k')
            else
                % Plot entire trajectory
                plot(pos(:,1),pos(:,2),'color',cTemp,...
                    'LineWidth',lineWidth,'LineStyle',lineStyleTemp)
                % Plot start position
                if ~isempty(startMarker)
                    plot(pos(1,1),pos(1,2),startMarker,'MarkerSize',...
                        markerSize/2,'MarkerFaceColor',cTemp,...
                        'MarkerEdgeColor', 'None')
                end
                % Plot end position
                if ~isempty(endMarker)
                    if strcmp(endMarker,'arrow')
                        plt.drawArrowMarker(pos(:,1),pos(:,2),cTemp,2*markerSize);
                    else
                        plot(pos(end,1),pos(end,2),endMarker,'MarkerSize',...
                            markerSize,'MarkerFaceColor',cTemp,'MarkerEdgeColor','k')
                    end
                end
            end
        else
            if verbose
                fprintf('No position data for Trial %d\n',TD(i).trialID)
            end
        end
    end
end
set(gca,'Box','on')

if ~isempty(axesLimits)
    set(gca,'XLim',axesLimits(1:2),'YLim',axesLimits(3:4))    
end
axis square
axis off