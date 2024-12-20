% flow.plot_all_session_summary  Creates a summary scatter plot of all
% sessions comparing A_int_rot and A_pred_rot. Significant markers will be
% filled and the NS markers will be hallow. It is possible to highlight an
% example session as well.
%
% Usage:
%   h = flow.batch_create_session_figs()
%
% Optional Inputs:
%   highlight_exp      Highlight example sessions
%   data_save_loc      Location to save flow field session data
%   fig_save_loc       Location to save flow field session figures
%   pltLegend          Determine whether to plot the legend
%   plotOverlap        Determine whether to plot a scatter of the number of
%                           voxel overlap
%   xTextOffset        Distance between the text printed on the screen
%   measMet            Plot the median and mean metric values
%   metric_str         Metric name in the structure
%   metric_plot_title  Metric name for printing text 
%
% Outputs:
%   h               Summary figure
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function h = plot_all_session_summary(varargin)

% Optional agruments
highlight_exp = [];
data_save_loc = [];
fig_save_loc = [];
pltLegend = 1;
plotOverlap = 0;
xTextOffset = 0.6;
measMet = 'median'; % Will either plot the median and mean metric values

% Define metrics to analyze
metric_str = {'mse'};
metric_plot_title = {'Mean squared error'};

assignopts(who, varargin);

% Load all data
FR = flow.load_session_results('data_save_loc',data_save_loc);
n_metric = length(metric_str);

% Determine index to highlight
if ~isempty(highlight_exp)
    hlt_subject = highlight_exp(1:(end-8));
    hlt_dataset = highlight_exp((end-7):end);
    all_subject = {FR.subject};
    all_dataset = {FR.dataset};
    subject_mask = ismember(all_subject, hlt_subject);
    dataset_mask = ismember(all_dataset, hlt_dataset);
    hlt_idx = find(subject_mask & dataset_mask);
end

% Set up figure
n_row = 1;
if plotOverlap
    n_col = n_metric + 1;
else
    n_col = n_metric;
end
axSz = 400;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(n_row,n_col,axSz,axSz,axSp);

h = figure('Position',[10 10 fW fH]);
set(h,'Name',sprintf('FlowAnalysis_SessionSummary_%s',measMet))

% Get vectorized data
A_int_rot = [FR.A_int_rot];
A_pred_rot = [FR.A_pred_rot];
stats = [FR.stats];
p_vals = [stats.pvals];

% Get subject information
subAll = {FR.subject};
uniSubj = {'Monkey E','Monkey D','Monkey Q'};
marker = {'ko','md','cs'};
plotSub = 1;
% Iterate over error metrics
for i = 1:n_metric
    % Get data for all sessions
    metric_int_rot = [A_int_rot.(metric_str{i})];
    metric_pred_rot = [A_pred_rot.(metric_str{i})];
    err_int_rot = [metric_int_rot.(measMet)];
    err_pred_rot = [metric_pred_rot.(measMet)];
    p_metric = [p_vals.(metric_str{i})];
    p = [p_metric.p];

    % Plot
    plt.subplotSimple(n_row, n_col, i, 'Ax', Ax); hold on;
    if plotSub
        for n = 1:length(unique(subAll))
            mask = find(ismember(subAll,uniSubj(n)));
            if ismember(hlt_idx,mask)
                hlt_idxTmp = find(ismember(mask,hlt_idx));
                [ptsLab(n),LEGEND{n}] = plot_error(err_int_rot(mask), err_pred_rot(mask),...
                    p(mask),metric_plot_title{i},hlt_idxTmp,marker{n},n,~pltLegend);
            else
                [ptsLab(n),LEGEND{n}] = plot_error(err_int_rot(mask), err_pred_rot(mask),...
                    p(mask),metric_plot_title{i},[],marker{n},n,~pltLegend);
            end
            % Add the monkey ID to the legend
            LEGEND{n} = [uniSubj{n} ', ' LEGEND{n}];
        end
        
        % Get axis limits
        x_lim = get(gca, 'XLim');
        y_lim = get(gca, 'YLim');
        ax_lim = nan(1, 2);
        ax_lim(1) = 0;
        ax_lim(2) = max([x_lim, y_lim]);

        % Add the N value information here
        font_size = 14;
        text(ax_lim(1) + diff(ax_lim) * xTextOffset, ax_lim(2), ...
            sprintf('N = %d sessions', length(p)), ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top', ...
            'FontSize', font_size)
    else
        [ptsLab(n),LEGEND{n}] = plot_error(err_int_rot, err_pred_rot, p, metric_plot_title{i}, hlt_idx);
    end

end
legend(ptsLab,LEGEND,'location','best')

% Plot number of overlap points
if plotOverlap
    plt.subplotSimple(n_row, n_col, n_col, 'Ax', Ax); hold on;

    % Get data to plot
    n_valid_int_rot = [A_int_rot.n_valid];
    n_valid_pred_rot = [A_pred_rot.n_valid];
    p = zeros(1, length(n_valid_int_rot));

    % Plot -- use the error plotting function
    if plotSub
        for n = 1:length(unique(subAll))
            mask = find(ismember(subAll,uniSubj(n)));
            if ismember(hlt_idx,mask)
                hlt_idxTmp = find(ismember(mask,hlt_idx));
                [ptsLab(n),LEGEND{n}] = plot_error(n_valid_int_rot(mask), n_valid_pred_rot(mask),...
                    p(mask),'Number of overlap points', hlt_idxTmp,marker{n},n)
            else
                [ptsLab(n),LEGEND{n}] = plot_error(n_valid_int_rot(mask), n_valid_pred_rot(mask),...
                    p(mask),'Number of overlap points', [],marker{n},n)
            end
        end

        % Get axis limits
        x_lim = get(gca, 'XLim');
        y_lim = get(gca, 'YLim');
        ax_lim = nan(1, 2);
        ax_lim(1) = 0;
        ax_lim(2) = max([x_lim, y_lim]);

        % Add the N value information here
        text(ax_lim(1) + diff(ax_lim) * xTextOffset, ax_lim(2), ...
            sprintf('N = %d sessions', length(p)), ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top', ...
            'FontSize', font_size)
    else
        plot_error(n_valid_int_rot, n_valid_pred_rot, p, ...
            'Number of overlap points', hlt_idx)
    end

    %plot_error(n_valid_int_rot, n_valid_pred_rot, p, ...
    %   'Number of overlap points', hlt_idx)
    xlabel('Intuitive vs SM')
    ylabel('Predicted SM flow vs actual SM flow')
end

% Plot title
plt.plotTitle(sprintf('Energy landscape :: Flow analysis :: Session summary :: %s',measMet))

% Save figure
if ~isempty(fig_save_loc)
    saveFigurePDF(h, fig_save_loc)
end
end


% Function to generate error plot
function [ptsLab,LEGEND]  = plot_error(err_int_rot, err_pred_rot, p, plot_title, hlt_idx,marker,shift,pltText)

if nargin == 5
    marker = 'k';
    shift = 0;
    pltText = 1;
elseif nargin == 6
    shift = 0;
    pltText = 1;
elseif nargin == 7
    pltText = 1;
end
font_size = 14;

% Create p-value mask
p_mask = p < 0.05;

% Plot significant and non-significant points
scatter(err_int_rot(~p_mask), err_pred_rot(~p_mask), 10, marker, ...
    'MarkerFaceColor', 'none')
ptsLab = scatter(err_int_rot(p_mask), err_pred_rot(p_mask), 10, marker, ...
    'MarkerFaceColor', 'flat');

% Get axis limits
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
ax_lim = nan(1, 2);
ax_lim(1) = 0;
ax_lim(2) = max([x_lim, y_lim]);

% Format plot
set(gca, 'TickDir', 'out', 'XLim', ax_lim, 'YLim', ax_lim)
tick = get(gca, 'XTick');
set(gca, 'XTick', tick, 'YTick', tick)
y_ax = get(gca, 'YAxis');
x_ax = get(gca, 'XAxis');
y_ax.FontSize = font_size;
x_ax.FontSize = font_size;

% Plot highlighed point
if ~isempty(hlt_idx)
    % Get data to plot
    x_hlt = err_int_rot(hlt_idx);
    y_hlt = err_pred_rot(hlt_idx);
    p_hlt = p_mask(hlt_idx);

    % Determine marker face color
    if p_hlt  % Point was significant
        marker_face_col = 'r';
    else
        marker_face_col = 'none';
    end

    % Plot
    plot(x_hlt, y_hlt, 'o', ...
        'MarkerFaceColor', marker_face_col, ...
        'MarkerEdgeColor', 'r', ...
        'MarkerSize', 8)
end

% Plot diagonal and number of session
plot(ax_lim, ax_lim, 'k--')
n = length(p);
if pltText
    text(ax_lim(1) + diff(ax_lim) * 0.6, ...
        ax_lim(2) - diff(ax_lim) * 0.05*(1+shift), ...
        sprintf('N_{sig} = %d/%d', sum(p_mask),n), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', font_size);
end
LEGEND = sprintf('N_{sig} = %d/%d', sum(p_mask),n);

title(plot_title, 'FontSize', font_size)
xlabel('Intuitive vs rotated error', 'FontSize', font_size)
ylabel('Predicted rotated vs actual rotated error', 'FontSize', font_size)
end