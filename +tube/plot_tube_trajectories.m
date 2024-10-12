% tube.plot_tube_trajectories  Plot cursor trajectories for constrained
% path experiment.
%


function [fh] = plot_tube_trajectories(IT, varargin)

r_factor = 0.15;  % Factor used for spatial averaging
center_pos = [0, 0, 0];
plot_scale = 220;
plotStates = {'Step 2'};
avgMode = 'spatial';
plotSegment = [];
C = [];
trialsPerCondition = [];     % Allow for subselection of trials
rmvFail = 1;

assignopts(who, varargin);

n_tube = length(IT.constrainedTubeRadius);
tube_radius = IT.constrainedTubeRadius;

if isempty(C)
    C = util.defineTaskColormap('bc_rot');
    C.targPos = C.targPos.*(norm(IT.startTargPos)/100)/.9;
end

% Initialize figure handles
fh = nan(1, 2);

% Setup figure -- all trajectories
nRow = 3;
nCol = n_tube;
axSp = 5;
axW = 300;
[fW, fH, Ax1] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
fh(1) = figure('Position',[10 10 fW fH]);
figName = sprintf('%s%s_IntTargConstAll_%0.2d', IT.subject, IT.date, ...
    IT.intermediate_target_num);
set(fh(1), 'Name', figName)
titleStr = sprintf('%s %s : Intermediate constrained trajectories (all)', ...
    IT.subject,IT.date);
plt.plotTitle(titleStr)

% Iterate over presented tubes and plot successful trajectories
TD_avg = cell(1, n_tube);
TD_best_avg = cell(1, n_tube);
TD_unconst_avg = cell(1, n_tube);
n_avg = nan(3, n_tube);
for i = 1:n_tube
    % Plot all trajectories
    plotIdx1 = i;
    plt.subplotSimple(nRow,nCol,plotIdx1,'Ax',Ax1);
    
    % Plot trials
    TD = IT.TDconstrained{i};
    util.plotIntermediateTargetTaskTrajectories(TD, C, '', ...
        'r_factor', r_factor, ...
        'avgMode', avgMode, ...
        'centerPos', center_pos, ...
        'plotTubes', true, ...
        'plotScale', plot_scale, ...
        'createFigure', false, ...
        'plotStates',plotStates,...
        'plotSegment',plotSegment,...
        'trialsPerCondition',trialsPerCondition,...
        'rmvFail',rmvFail);
    plt.plotDatasetInfo(TD(1).subject,TD(1).date,[TD.trialID],'plotMode','axis')
    if i == 1; ylabel('All trials','FontSize', 16, 'FontWeight', 'Bold'); end
    title(sprintf('Tube r = %d',tube_radius(i)),'FontWeight','bold', ...
        'FontSize',16)
    
    % Plot best block of trials
    plotIdx1 = n_tube + i;
    plt.subplotSimple(nRow,nCol,plotIdx1,'Ax',Ax1);
    
    % Plot trials
    TD_best = IT.TD_tube_best{i};
    util.plotIntermediateTargetTaskTrajectories(TD_best,C,'', ...
        'r_factor', r_factor, ...
        'avgMode', avgMode, ...
        'centerPos', center_pos, ...
        'plotTubes', true, ...
        'plotScale', plot_scale, ...
        'createFigure', false, ...
        'plotStates',plotStates,...
        'plotSegment',plotSegment,...
        'trialsPerCondition',trialsPerCondition,...
        'rmvFail',rmvFail);
    plt.plotDatasetInfo(TD_best(1).subject, TD_best(1).date, ...
        [TD_best.trialID],'plotMode','axis')
    if i == 1
        label_str = sprintf('Best %d trials', IT.block_size);
        ylabel(label_str,'FontSize',16,'FontWeight','Bold');
    end
    
    % Plot unconstrained trials with tube constraint applied
    plotIdx1 = n_tube * 2 + i;
    plt.subplotSimple(nRow,nCol,plotIdx1,'Ax',Ax1);
    
    % Plot trials
    TD_unconst = IT.TD_unconst_tube{i};
    util.plotIntermediateTargetTaskTrajectories(TD_unconst,C,'', ...
        'r_factor', r_factor, ...
        'avgMode', avgMode, ...
        'avgMode', avgMode, ...
        'centerPos', center_pos, ...
        'plotTubes', true, ...
        'plotScale', plot_scale, ...
        'createFigure', false, ...
        'plotStates',plotStates,...
        'plotSegment',plotSegment,...
        'trialsPerCondition',trialsPerCondition,...
        'rmvFail',rmvFail);
    plt.plotDatasetInfo(TD_unconst(1).subject, TD_unconst(1).date, ...
        [TD_unconst.trialID],'plotMode','axis')
    if i == 1
        label_str = 'Unconstrained trials';
        ylabel(label_str,'FontSize',16,'FontWeight','Bold');
    end
    
    % Average all trials
    sc_mask = logical([TD.successful]);
    TD = TD(sc_mask);
    TD_avg{i} = TD.average('avgMode', avgMode, 'r_factor', r_factor);
    n_avg(1, i) = sum(sc_mask);
    
    % Average best trials
    sc_mask = logical([TD_best.successful]);
    TD_best = TD_best(sc_mask);
    TD_best_avg{i} = TD_best.average('avgMode', avgMode, 'r_factor', r_factor);
    n_avg(2, i) = sum(sc_mask);
    
    % Average unconstrained trials
    sc_mask = logical([TD_unconst.successful]);
    TD_unconst = TD_unconst(sc_mask);
    TD_unconst_avg{i} = TD_unconst.average('avgMode', avgMode, 'r_factor', r_factor);
    n_avg(3, i) = sum(sc_mask);
end

%--------------------------------------------------------------------------
% Plot average trajectories for each tube size
% Setup figure
nRow = 1;
nCol = 3;
axSp = 5;
axW = 300;
[fW, fH, Ax1] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
fh(2) = figure('Position',[10 10 fW fH]);
figName = sprintf('%s%s_IntTargConstAvg_%0.2d', IT.subject, IT.date, ...
    IT.intermediate_target_num);
set(fh(2), 'Name', figName)
titleStr = sprintf('%s %s : Intermediate constrained trajectories (average)', ...
    IT.subject,IT.date);
plt.plotTitle(titleStr)

% Get tube colors
col = tube.get_tube_col(IT.startTargPos, n_tube,'C',C);

% Combine all and best trajectories into one cell array
TD_avg_combined = {TD_avg, TD_best_avg, TD_unconst_avg};
titleStr = cell(1, 3);
titleStr{1} = 'All trials';
titleStr{2} = sprintf('Best %d trials', IT.block_size);
titleStr{3} = 'Unconstrained trials';
n_panels = length(TD_avg_combined);

% Loop over conditions (all trials/best trials)
for i = 1:n_panels
    plt.subplotSimple(nRow,nCol,i,'Ax',Ax1);
    TDavg_temp = TD_avg_combined{i};

    % Initialize arrays for legend
    LH = nan(1, n_tube);
    leg_str = cell(1, n_tube);

    % Iterate over tubes and plot
    for j = 1:length(TDavg_temp)
        if j == 1
            % Plot targets
            TDavg_temp{j}.plot( ...
                'color', col.gray, ...
                'plotTrajectories', 0, ...
                'plotStates', {'Step 2'}, ...
                'plotIntTarg', false, ...
                'plotTargNum', false)
            TDavg_temp{j}.plot( ...
                'color', col.tube(1,:), ...
                'plotStartTarg', true, ...
                'plotEndTarg', false)
        end
        TDavg_temp{j}.plot( ...
            'color', col.tube(j+1,:), ...
            'plotTargets', 0, ...
            'lineWidth', 3)
        IT.tubeObject(j).plot('color', col.tube(j+1,:),...
            'lineSpec', '-','plotSegment',plotSegment);

        % Get line object handle and legend text
        LH(j) = plot(nan, nan, '-', 'color', col.tube(j+1,:), 'lineWidth',3);
        leg_str{j} = sprintf('r = %d mm (%d trials)', tube_radius(j), ...
            n_avg(i, j));
    end

    % Format average trajectory figure
    if ~isempty(plot_scale)
        axis([-1 1 -1 1]*plot_scale)
    end
    plt.scaleBar(gca,20,'mm')

    % Plot legend
    legend(LH, leg_str, 'box', 'off', ...
        'Location', 'SouthWest', ...
        'FontSize', 12)

    axis on
    set(gca,'XTick',[],'YTick',[])
    title(titleStr{i}, 'FontSize', 16, 'FontWeight','Bold');
end

