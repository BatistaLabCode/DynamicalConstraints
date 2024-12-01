function [F_ff, h2] = fig_5_plot_flow_field_analysis(dataLoc,varargin)
% [F_ff, h2] = fig_5_plot_flow_field_analysis(dataLoc) Plots the flow field
% the two-target task of example sessions in four different conditions: 
%      1) MoveInt projection feedback trials, MoveInt projection (Left-top)
%      2) MoveInt projection feedback trials, SepMax projection (Right-top)
%      3) SepMax projection feedback trials, MoveInt projection (Left-bot.)
%      4) SepMax projection feedback trials, SepMax projection (Right-bot.)
% 
% Additionally, it plots the summary scatter plot showing the difference in
% flow fields when comparing 1-2 and 2-4 from the list above.
%
% Inputs:
%   dataLoc    Paths for the main data folder
%
% Optional Inputs:
%   exampleSess         The example session used in the paper (fig6/7).
%   saveFig             Determine whether or not to save the data
%   savePathBase        Where to save the figures.
%   arrowSize           Size of the arrow heads to plot on the flowfield
%
% Outputs:
%   F_ff               Figure of the flow fields plotted in 2x2 grid.
%   h2                 Figure scatter plot for all sesssions. Example
%                       sessions are highlighted
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Load in the D structure and the data
exampleSess = {'20190719'}; % The example session used in the paper.
saveFig = 0;                % Determine whether or not to save the data,
                            % default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
arrowSize = 7;              % Size of the arrow heads to plot on the flowfield

% Assign the optional inputs
assignopts(who,varargin);

data_save_loc = fullfile(dataLoc,'flowAnalysis','mat'); % Where the flow field data is saved.

% Load the flow fields
FlowResults = flow.load_session_results('data_save_loc',data_save_loc);

% Determine the sessions with the correct data
mask = ismember({FlowResults.dataset},exampleSess);
tFR = FlowResults(mask == 1);

% Iterate through the different examples.
for k = 1:size(tFR,2)
    
    % Create a summary trajectory figures for the session.
    axSp = 75;
    axW = 400;
    [fW, fH, Ax] = plt.calcFigureSize(2,2,axW,axW,axSp);

    F_ff(k) = figure('Position',[10 10 fW fH]); % Flow field plot
    F_ff(k).Name = sprintf('%s%s_flowFieldSummary',...
        tFR(k).subject, tFR(k).dataset);
    plt.plotTitle(sprintf("Flow Field Summary Plot: %s %s",...
        tFR(k).subject, tFR(k).dataset))

    % Make the subplots for the flow field
    for n = 1:4
        figure(F_ff(k))
        axFF_Fin(n) = plt.subplotSimple(2,2,n,'Ax',Ax); % Set the new axes positions
        axis off
    end

    % Create a grid of condition information and flow fields.
    condStr = {'MoveInt projection feedback trials, MoveInt projection',...
        'MoveInt projection feedback trials, SepMax projection';...
        'SepMax projection feedback trials, MoveInt projection',...
        'SepMax projection feedback trials, SepMax projection'};
    FF = [tFR(k).FF_int tFR(k).FF_pred; tFR(k).FF_pred_Int tFR(k).FF_rot];

    for n = 1:2 % Loop through feedback condition (1=MovInt,2=SepMax)
        for i = 1:2 % Loop through projection condition (1=MovInt,2=SepMax)
            if n == 1
                C = util.defineTaskColormap('bc_int');
            else
                C = util.defineTaskColormap('bc_rot');
            end
            % Adjust the color target scale to match the targets
            scale = norm(FF(n,i).kin.start_pos(1,1:2))./100;
            C.targPos = C.targPos.*(scale/.9);

            % Plot the flow field session
            [fh] = flow.plot_flow_field(FF(n,i),'C',C,...
                'arrow_size',arrowSize);

            % Add the target pair plot to the flow field summary figure
            for m = 1:2
                axFF(n+2*(i-1),m) = copyobj(fh(1).Children(3+2*(m-1)), F_ff(k));
                axFF(n+2*(i-1),m).Position = axFF_Fin(i+2*(n-1)).Position;
                axes(axFF(n+2*(i-1),m))
                h_cb = colorbar;
                h_cb.Ticks = linspace(0,1,10);
                h_cb.TickLabels = linspace(0,90,10);
                h_cb.Label.String = 'Number of time points';
                h_cb.Label.Rotation = 270;
                h_cb.Label.FontSize = 14;
                if m == 2 && i == 2
                    h_cb.Visible = 'off';
                end
                if m == 2
                    axFF(n+2*(i-1),m).Color = 'none';
                    title(condStr{n,i},'Interpreter','none')
                else
                    title('')
                end 
                axis square
            end
            close([fh(:)])
        end
    end

    % Plot the summary stats
    h2(k) = flow.plot_all_session_summary('data_save_loc',data_save_loc,...
        'metric_str',{'mse'},'metric_plot_title',{'Mean squared error'},...
        'highlight_exp',[tFR(k).subject, tFR(k).dataset]);   
end

% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF([F_ff(:); h2(:)],savePathBase)
end
end