function [F] = fig_4_plot_view_vs_projections(dataLoc,varargin)
% ADD the header
%
%
%
%
%
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Load and choose the data
exampleSess = {'20190719'}; % Figure 4 example session
saveFig = 0;                % Determine whether or not to save the data,
                            % default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
plotScale = 180;            % Axis scale for the plot panels
centerPos = [0 0 0];
plotStates = {'Step 2','BC Step 2'}; % Trial states to plot
avgMode = 'samp';                    % Average mode, default is by sample 
                                     % and not with a spatial temporal warp
                                     % ('spatial')
endMarker = 'arrow';                 % Plots an arrow at the end of the 
                                     % trajectories

% Assign the optional inputs
assignopts(who,varargin);

% Determine the sessions with the correct data
load(fullfile(dataLoc,'publicationQualitySessions.mat'))
D = D(ismember({D.dataset},exampleSess));
dir_list = db.get_task_datasets(D, {'tt_int','tt_rot'});

for k = 1:size(dir_list,1)
    % Load the example trajectory data
    [TD_int, P_int, result_int] = util.loadSessionData(dir_list(k,1));
    [TD_rot, P_rot, result_rot] = util.loadSessionData(dir_list(k,2));

    % Normalize the data
    uniTD_intPos = unique([TD_int.startPos]','rows');
    uniTD_rotPos = unique([TD_rot.startPos]','rows');
    if ~ismember(mean(unique([TD_rot.startPos]','rows')),centerPos)
        TD_int = TD_int.normalize(mean(uniTD_rotPos));
        TD_rot = TD_rot.normalize(mean(uniTD_rotPos));
    else
        TD_int = TD_int.normalize(centerPos);
        TD_rot = TD_rot.normalize(centerPos);
    end

    % Exclude any trials that were not part of the two target pair.
    TD_int = TD_int(ismember([TD_int.targPos]',uniTD_intPos,'rows'));
    TD_rot = TD_rot(ismember([TD_rot.targPos]',uniTD_rotPos,'rows'));

    % Create a summary trajectory figures for the session.
    axSp = 75; % Space between the axis
    axW = 300; % Size of the axis
    [fW, fH, Ax] = plt.calcFigureSize(2,2,axW,axW,axSp);

    F(k) = figure('Position',[10 10 fW fH]); % Trajectory plot
    F(k).Name = sprintf('%s%s_TwoTargetSummary',...
        dir_list(k,1).subject, dir_list(k,1).dataset);
    plt.plotTitle(sprintf("Summary Plot: %s %s",dir_list(k,1).subject,...
        dir_list(k,1).dataset));
    axFin = repmat(axes,4,1);
    for n = 1:4
        figure(F(k))
        axFin(n) = plt.subplotSimple(2,2,n,'Ax',Ax); % Set the new axes positions
        axis off
    end

    % Create the projected data
    dir_names = {dir_list(k,1).trajectory{1}, dir_list(k,2).trajectory{1}};
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
                    markerSize = [1 10];
                    plotTargets = 1;
                    condStr = 'MI projection feedback trials, MI projection';
                else
                    startMarker = [];
                    markerSize = [0.1 10];
                    plotTargets = 0;
                    condStr = 'MI projection feedback trials, SM projection';
                end
            else
                TD = TD_rot;
                C = util.defineTaskColormap('bc_rot');
                if n == 2
                    startMarker = 'o';
                    markerSize = [1 10];
                    plotTargets = 1;
                    condStr = 'SM projection feedback trials, SM projection';
                else
                    startMarker = [];
                    markerSize = [0.1 10];
                    plotTargets = 0;
                    condStr = 'SM projection feedback trials, IM projection';
                end
            end

            % File name
            save_name_base = [f_name(1:end-18) '_' decode_name_file];

            % Add the GPFA data for intuitive trials
            TD = util.predictDecodeState_GPFA(TD,P,'spikeCountSrc','decodeSpikeCounts',...
                'useT',1,'TT',TT,'predictPos',1);

            % Normalized and set the target center
            [TDnorm,targInfo] = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

            % Adjust the color target scale to match the targets
            scale = norm(TDnorm(1).startPos(1:2))./100;
            C.targPos = C.targPos.*(scale/.9);

            % Plot the cursor trajectories
            f_cursor = util.plotGridTaskTrajectories(TDnorm,save_name_base, ...
                'ColMat', C, ...
                'plotScale', plotScale,...
                'plotStates',plotStates,...
                'avgMode',avgMode,...
                'startMarker',startMarker,...
                'endMarker',endMarker,...
                'plotTargets',plotTargets,...
                'trlMarkerSz',markerSize(1),...
                'avgMarkerSz',markerSize(2));

            % Add the target pair plot to the two target summary figure
            ax(n+2*(i-1)) = copyobj(f_cursor(1).Children(1), F(k));
            ax(n+2*(i-1)).Position = axFin(n+2*(i-1)).Position; 
            axes(ax(n+2*(i-1)))
            if i ~= n
                axis off
            end
            title(condStr)

            % Close f_cursor
            close(f_cursor)
        end
    end 
    plt.matchAxis(ax)
end
% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF([F(:)],savePathBase)
end
end