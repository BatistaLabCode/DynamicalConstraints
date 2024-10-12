% flow.plot_all_session_histogram

function h = plot_all_session_histogram(varargin)

% Optional arguments
highlight_exp = [];
data_save_loc = '/Volumes/Samsung_T5/Analysis/EnergyLandscape/FlowComparison/mat';
fig_save_loc = '/Volumes/Samsung_T5/Analysis/EnergyLandscape/FlowComparison/fig';
measMet = 'median'; % Will either plot the median and mean metric values
sepSig = 0;
assignopts(who, varargin);

% Load all data
FR = flow.load_session_results('data_save_loc',data_save_loc);

% Define metrics to analyze
metric_str = {'ang', 'mag', 'mse'};
metric_plot_title = {...
    'Angular error', ...
    'Magnitude error', ...
    'Mean squared error'};
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
n_col = n_metric + 1;
axSz = 300;
axSp = 75;
[fW,fH,Ax] = plt.calcFigureSize(n_row,n_col,axSz,axSz,axSp);

h = figure('Position',[10 10 fW fH]);
set(h,'Name',sprintf('FlowAnalysis_SessionSummaryHistogram_%s',measMet))

% Get vectorized data
A_int_rot = [FR.A_int_rot];
A_pred_rot = [FR.A_pred_rot];
stats = [FR.stats];
p_vals = [stats.pvals];

% Collect and plot the histograms of the different conditions
for n = 1:length(metric_str)
    % Get data for all sessions
    metric_int_rot = [A_int_rot.(metric_str{n})];
    metric_pred_rot = [A_pred_rot.(metric_str{n})];

    % Determine the binWidth size
    range = [metric_pred_rot.mean metric_int_rot.mean];
    stepSz = 3*std(range)./15;
    stepSz = min([stepSz 30]);

    % Calculate the paired t-test on the data
    [hVal,pVal] = ttest([metric_int_rot.(measMet)],[metric_pred_rot.(measMet)]);

    % Plot the histograms
    figure(h)
    plt.subplotSimple(n_row,n_col,n,'Ax',Ax); hold on
    histogram([metric_int_rot.(measMet)],'BinWidth',stepSz,'DisplayStyle','stairs')
    plot(metric_int_rot(hlt_idx).(measMet),0,'V')
    histogram([metric_pred_rot.(measMet)],'BinWidth',stepSz,'DisplayStyle','stairs')
    plot(metric_pred_rot(hlt_idx).(measMet),0,'^')
    title(sprintf('%s \n h = %1.0f, p = %1.3d',metric_plot_title{n},hVal,pVal))
    xlabel('Error'),ylabel('Session Count')
    legend('Int vs SM','Int vs SM Example','Pred SM vs Actual SM','Pred SM vs SM Example')

    
end

% Get data to plot
n_valid_int_rot = [A_int_rot.n_valid];
n_valid_pred_rot = [A_pred_rot.n_valid];

% Calculate the paired t-test on the data
[hVal,pVal] = ttest(n_valid_int_rot,n_valid_pred_rot);

% Determine the binWidth size
range = [n_valid_int_rot n_valid_pred_rot];
stepSz = 3*std(range)./15;
stepSz = min([stepSz 30]);

% Plot the number of overlap points
plt.subplotSimple(n_row, n_col, n_col, 'Ax', Ax); hold on
histogram(n_valid_int_rot,'BinWidth',stepSz,'DisplayStyle','stairs')
plot(A_int_rot(hlt_idx).n_valid,0,'V')
histogram(n_valid_pred_rot,'BinWidth',stepSz,'DisplayStyle','stairs')
plot(A_pred_rot(hlt_idx).n_valid,0,'^')
xlabel('Overlap points'),ylabel('Session Count')
legend('Int vs Rot','Int vs Rot Example','Pred vs Rot','Pred vs Rot Example')
title(sprintf('Number of overlap points \n h = %1.0f, p = %1.3d',hVal,pVal))

% Plot title
plt.plotTitle(sprintf('Energy landscape :: Flow analysis :: Session histogram summary :: %s',measMet))

% Save figure
saveFigurePDF(h, fig_save_loc)
end