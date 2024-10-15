function [F_ff, h2] = fig_5_plot_flow_field_analysis(dataLoc,varargin)
% [Fhist,FintAng,Ftube,FintAngTube,F_EL] = fig_6_fig_7(dataLoc) Plots the 
% SepMax trajectories for the intermediate target task and the constrained
% target task. Additionally, it plots the summary histograms showing the
% change in initial angel for different conditions.
%
% Inputs:
%   dataLoc    Paths for the main data folder
%
% Optional Inputs:
%   exampleSess         The example session used in the paper (fig6/7).
%   saveFig             Determine whether or not to save the data
%   savePathBase        Where to save the figures.
%   C                   Colormap struct to use
%   alpha               Alpha for statistical tests (early vs late comparison)
%   trlCnt              Number of trials to use for early vs late comparison
%   plotExamples        Plots the example sessions on the histograms
%   avgMode             Average method for the trajectories
%   plotScale           Axis Limits
%   trialsPerCondition  Number of trials to subselect for plotting
%   setRandSeed         Fix the random number generator for a consistenet
%                           permutation for trial subselection.
%   useExample          Use the hardcode example trajectories (what is in
%                           the figure)
%
% Outputs:
%   Fhist               Figure histogram of initial angle summaries for all
%                           sessions and animals.
%   FintAng             Figure initial angle intermediate target task
%   Ftube               Figure trial examples of trajectories for the tube
%                           task.
%   FintAngTube         Figure Initial angle tube example
%   F_EL                Figure Early vs Late trajectory example
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Load in the D structure and the data
exampleSess = {'20190719'}; % The example session used in the paper.
saveFig = 0;                % Determine whether or not to save the data,
                            % default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
plotScale = 200;               % Axis Limits plotScal*[-1 1 -1 1]
centerPos = [0 0 0];           % Defined center of the workspace
avgMode = 'samp';              % Average method for the trajectories.
plotStates = {'Step 2','BC Step 2'}; % Which task states to include in flow field calculation
d_grid = 20;                % Allow the flow field grid to be pre-specified
arrowSize = 7;              % Size of the arrow heads to plot on the flowfield
min_pts_per_voxel = 2;      % Minimum number of points in voxel to calculate flow field.

% Assign the optional inputs
assignopts(who,varargin);

F_ffcb =[];

% Determine the sessions with the correct data
load(fullfile(dataLoc,'publicationQualitySessions.mat'));
D = D(ismember({D.dataset},exampleSess));
dir_list = db.get_task_datasets(D, {'tt_int','tt_rot'});

data_save_loc = fullfile(dataLoc,'flowAnalysis','mat'); % Where the flow field data is saved.

% Iterate through the different conditions.
for k = 1:size(dir_list,1)
    [TD_int, P_int, result_int] = util.loadSessionData(dir_list(k,1));
    [TD_rot, P_rot, result_rot] = util.loadSessionData(dir_list(k,2));

    % Normalize the data
    if ~ismember(mean(unique([TD_rot.startPos]')),centerPos)
        TD_int = TD_int.normalize(mean(unique([TD_rot.startPos]','rows')));
        TD_rot = TD_rot.normalize(mean(unique([TD_rot.startPos]','rows')));
    else
        TD_int = TD_int.normalize(centerPos);
        TD_rot = TD_rot.normalize(centerPos);
    end

    % Exclude any trials that were not part of the two target pair.
    TD_int = TD_int(ismember([TD_int.targPos]',unique([TD_int.startPos]','rows'),'rows'));
    TD_rot = TD_rot(ismember([TD_rot.targPos]',unique([TD_rot.startPos]','rows'),'rows'));

    % Create a summary trajectory figures for the session.
    axSp = 75;
    axW = 400;
    [fW, fH, Ax] = plt.calcFigureSize(2,2,axW,axW,axSp);

    F_ff(k) = figure('Position',[10 10 fW fH]); % Flow field plot
    F_ff(k).Name = sprintf('%s%s_flowFieldSummary',...
        dir_list(k,1).subject, dir_list(k,1).dataset);
    plt.plotTitle(sprintf("Flow Field Summary Plot: %s %s",dir_list(k,1).subject,...
        dir_list(k,1).dataset))

    % Make the subplots for the flow field
    for n = 1:4
        figure(F_ff(k))
        axFF_Fin(n) = plt.subplotSimple(2,2,n,'Ax',Ax); % Set the new axes positions
        axis off
    end

    % Create the projected data
    dir_names = {dir_list(k,1).trajectory{1}, dir_list(k,2).trajectory{1}};
    TD_cell = cell(0);
    for n = 1:2 % Result and decoder loop
        % Define decoder
        if n == 1
            result = result_int;
            P = P_int;
            decode_name_file = 'IntDecoder';
        else
            result = result_rot;
            P = P_rot;
            decode_name_file = 'SMDecoder';
        end

        % Find orthonormal
        [~,~,TT] = orthogonalize(zeros(P.xDim,1),result.estParams.C);

        for i = 1:2 % Trial loop
            [~, f_name] = fileparts(dir_names{i});
            if i == 1
                TD = TD_int;
                C = util.defineTaskColormap('bc_int');
                if n == 1
                    startMarker = 'o';
                    endMarker = 'o';
                    plotTargets = 1;
                    condStr = 'MoveInt projection feedback trials, MoveInt projection';
                else
                    startMarker = [];
                    endMarker = [];
                    plotTargets = 0;
                    condStr = 'MoveInt projection feedback trials, SepMax projection';
                end
            else
                TD = TD_rot;
                C = util.defineTaskColormap('bc_rot');
                if n == 2
                    startMarker = 'o';
                    endMarker = 'o';
                    plotTargets = 1;
                    condStr = 'SepMax projection feedback trials, SepMax projection';
                else
                    startMarker = [];
                    endMarker = [];
                    plotTargets = 0;
                    condStr = 'SepMax projection feedback trials, MoveInt projection';
                end
            end

            % File name
            save_name_base = [f_name(1:end-18) '_' decode_name_file];

            % Add the GPFA data for intuitive trials
            TD = util.predictDecodeState_GPFA(TD,P,'spikeCountSrc','decodeSpikeCounts',...
                'useT',1,'TT',TT,'predictPos',1);
            TD_cell{n,i} = TD;

            % Normalized and set the target center
            TDnorm = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

            % Adjust the color target scale to match the targets
            scale = norm(TDnorm(1).startPos(1:2))./100;
            C.targPos = C.targPos.*(scale/.9);

            % Run flow field analysis
            [FF(n,i)] = flow.calc_flow_field(TD, d_grid, save_name_base, ...
                'min_pts_per_voxel', min_pts_per_voxel);

            % Plot the flow field session
            [fh] = flow.plot_flow_field(FF(n,i),'C',C,...
                'arrow_size',arrowSize);

            % Add the target pair plot to the flow field summary figure
            for m = 1:2
                axFF(n+2*(i-1),m) = copyobj(fh(1).Children(3+2*(m-1)), F_ff(k));
                axFF(n+2*(i-1),m).Position = axFF_Fin(n+2*(i-1)).Position;
                axes(axFF(n+2*(i-1),m))
                h_cb = colorbar;
                h_cb.Ticks = linspace(0,1,10);
                h_cb.TickLabels = linspace(0,90,10);
                h_cb.Label.String = 'Number of time points';
                h_cb.Label.Rotation = 270;
                h_cb.Label.FontSize = 14;
                if m == 2 && n == 2
                    h_cb.Visible = 'off';
                end
                if m == 2
                    axFF(n+2*(i-1),m).Color = 'none';
                    title(condStr,'Interpreter','none')
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
        'highlight_exp',[dir_list(k,1).subject, dir_list(k,1).dataset]);   
end

% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF([F_ff(:); h2(:)],savePathBase)
end
end