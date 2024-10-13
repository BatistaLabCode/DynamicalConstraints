function [Fhist,FintAng,Ftube,FintAngTube,F_EL] = fig_6_fig_7(dataLoc,varargin)
%% ADD the header
%
%
%
%
%
%
% Created by Erinn Grigsby (erinn.grigsby@gmail.com)

%% Load in the D structure and the data
exampleSess = {'20190719'}; % The example session used in the paper (fig2).
saveFig = 0;                % Determine whether or not to save the data,
% default is to not save the data (0).
savePathBase = [];          % Where to save the figures.
plotExamples = 0; % Plots the example sessions on the histograms
alpha = 0.05; % Alpha for statistical tests (early vs late comparison)
trlCnt = 20; % Number of trials to use for early vs late comparison

plotScale = 220;               % Axis Limits plotScal*[-1 1 -1 1]
avgMode = 'samp';              % Average method for the trajectories.
C = util.defineTaskColormap('bc_rot');  % Colormap struct. Default bc_rot
trialsPerCondition = 10;     % Allow for subselection of trials
setRandSeed = 9; % Allow to fix the random permutation for trial subselection.
useExample = 1;              % Use the hardcode example trajectories

% Assign the optional inputs
assignopts(who,varargin);

% Load the initial angle analysis data
% AD_compare , ang_antFlow, ang_unCon, ang_flow, and AD_co
load(fullfile(dataLoc,'initialAng_tt_vs_uncon+con_Simplified.mat'))
load(fullfile(dataLoc,'centerOut_AngleAnalysisData_Signed.mat'))

% Load the location where the intermediate target data is saved.
data_save_loc = fullfile(dataLoc,'ConstrainedPath\mat\int_targ_data');

%% Creates the histogram plots for figure 6 and figure 7
% Create normIntAng variable
data = nan(size(AD_compare,1),5);
for n = 1:size(AD_compare,1)
    data(n,1) = AD_compare(n,1).compAng(2);
    data(n,2) = AD_compare(n,1).avgUncon(2);
    data(n,3) = AD_compare(n,1).avgCon(2,end);
    data(n,4) = std(AD_compare(n,1).rotTT(2,:));
    data(n,5) = std(AD_compare(n,1).uncon(2,:));
end
avgCon = data(:,3);

% Calculate the early vs late control angle differencess
datNoChange = repmat(struct('subject',[],'date',[],...
    'idx_E',[],'avg_E',[],'ci_E',[],'idx_L',[],'avg_L',[],'ci_L',[],...
    'stat_EL',[],'cmpNorm_EL',[],'pEL',[],'hEL',[]),1,size(AD_compare,1));
for n = 1:size(AD_compare,1)
    % Add the standard information
    datNoChange(n).subject = AD_compare(n,1).subject;
    datNoChange(n).date = AD_compare(n,1).dataset;

    % Collect the two target data in the rotated plan
    ttAng = AD_compare(n,1).rotTT(2,:);

    % Collect the early and late data
    idxE = 1:trlCnt;
    idxL = (length(ttAng)-(trlCnt-1)):length(ttAng);
    datNoChange(n).idx_E = idxE;
    avg_E = mean(ttAng(idxE),'omitnan');
    datNoChange(n).avg_E = avg_E;
    datNoChange(n).ci_E = quantile(ttAng(idxE), [alpha/2, 1-(alpha/2)]);
    datNoChange(n).idx_L = idxL;
    avg_L = mean(ttAng(idxL),'omitnan');
    datNoChange(n).avg_L = avg_L;
    datNoChange(n).ci_L = quantile(ttAng(idxL), [alpha/2, 1-(alpha/2)]);
    datNoChange(n).cmpNorm_EL = (avg_E-avg_L)./avg_E;
    [pEL,hEL] = ranksum(ttAng(idxE),ttAng(idxL),alpha);
    datNoChange(n).pEL = pEL;
    datNoChange(n).hEL = hEL;
end

%% Create the normalize center out full change comparison.
ref = [8 1:8 1];

[tempCO,tempP45,tempN45] = deal(nan(size(AD_co))); % Collect the average trajectories
for n = 2:9
    % collect the center out data
    tmpCO = [AD_co(ref(n),:).avgAng_ZIA];
    tempCO(ref(n),:) = tmpCO(2,:);

    % Collect the negative offset degree target,remember negative for CO so it
    % will be the positive offset for the other target
    tmpPOff = [AD_co(ref(n-1),:).avgAng_ZIA];
    tempN45(ref(n),:) = tmpPOff(1,:);

    % Collect the positive offset degree target,remember positive for CO so it
    % will be the negative offset for the other target
    tmpPOff = [AD_co(ref(n+1),:).avgAng_ZIA];
    tempP45(ref(n),:) = tmpPOff(3,:);
end

g= (tempP45-tempCO)./tempP45;
g2= (tempN45-tempCO)./tempN45;
tempCO45 = .5*(g+g2);

%% Create the histograms for no change, full change, unconstrained, and
% constrained data. Include the mean and standard deviation in each plot,
% show the example session and list the p-values.
infoTitle ={'Unconstrained','No Change','Full Change','Constrained'};
collDat = {100*(ang_flow2intTarg_r-ang_unCon2intTarg_r)./ang_flow2intTarg_r,... % Unconstrained
    100*[datNoChange.cmpNorm_EL],... % No change (early vs late tt)
    100*median(tempCO45),... % Full Change (center out)
    (100*([datNoChange.avg_E]' - avgCon)./[datNoChange.avg_E]')'}; % Constrained

% Setup figure
nRow = 2;
nCol = 2;
axSp = 110;
axH = 250;
[fW, fH, Ax1] = plt.calcFigureSize(nRow,nCol,2*axH,axH,axSp);
Fhist = figure('Position',[10 10 fW fH]);
Fhist.Name = 'Updated_histograms_figure6+7';

% Iterate through the four conditions.
for n = 1:4
    % Plot all the animala data together
    figure(Fhist)
    plt.subplotSimple(nRow,nCol,n,'Ax',Ax1); hold on
    histogram(collDat{n},'BinWidth',5)
    xline([0 100],'k','LineStyle','--','LineWidth',2)
    xline(mean(collDat{n},'omitnan'),'LineStyle','-')
    xline(mean(collDat{n},'omitnan') + std(collDat{n},'omitnan')*[-1 1],'LineStyle',':')
    LEGEND = {'initial angle','o','100',...
        ['mean =' num2str(mean(collDat{n},'omitnan'))],...
        ['std =' num2str(std(collDat{n}))],''};

    % Plot the example sessions
    if plotExamples
        posSes = find(ismember({datNoChange.date},exampleSess));
        plot(collDat{n}(posSes),0,'o')
        LEGEND = [LEGEND,exampleSess];
    end

    % The run two t-tests, one comparing the data to the 0 value and the
    % other to the 100 value. This was to determine if there was a
    % significant difference between the real data distribution and the ideal
    % no change and full change conditions.
    [h,p] = ttest(collDat{n});
    [h100,p100] = ttest(100-collDat{n});
    if n == 1 || n == 4
        [hNC,pNC] = ttest(collDat{n},collDat{2});
        [hFC,pFC] = ttest(collDat{n},collDat{3});
        title(sprintf(['%s, To zero comparison p= %0.3d, h = %0.0f \n'...
            'To 100 comparison p= %0.3d, h=%0.0f \n'...
            'To NoChange p= %0.3d, h=%0.0f, To FullChange p= %0.3d, h=%0.0f'],...
            infoTitle{n},p,h,p100,h100,pNC,hNC,pFC,hFC));
    else
        title(sprintf(['%s, To zero comparison p= %0.3d, h = %0.0f \n'...
            'To 100 comparison p= %0.3d, h=%0.0f'],infoTitle{n},p,h,p100,h100));
    end
    legend(LEGEND)
    %legend('initial angle','o','100',['mean =' num2str(mean(collDat{n},'omitnan'))],...
    %    ['std =' num2str(std(collDat{n}))])
    xlabel('% difference in initial angle (relative to direct path)')
    ylabel('# Experiments')
end
plt.matchAxis(Fhist);

%% Plot the example session trajectory data
for n = 1:length(exampleSess)
    % Find and load the data
    valid_files = util.findDirContents(data_save_loc, [exampleSess{n} '_int_targ.mat']);
    load(fullfile(data_save_loc, valid_files{1}));

    scale = norm(IT.startTargPos(1:2))./100;
    C.targPos = C.targPos.*(scale/.9);
    %% Plot the trajectory examples for the angle difference
    TD_tt_rot = IT.TD_tt_rot{1};
    startPos = IT.TDconstrained{1}(1).startPos';

    % 6D) Calculate the angle difference for the unconstrained trials
    TD = IT.TDunconstrained;
    TD = TD(ismember([TD.startPos]',startPos,'rows'));
    [FintAng] = plot_intAngComparisons(TD,TD_tt_rot,startPos,C);
    title(sprintf('%s%s Average trajectories and initial angles \n Unconstrained',...
        IT.subject,IT.date))
    axis([-1 1 -1 1]*plotScale)
    axis on,set(gca,'XTick',[],'YTick',[])
    plt.scaleBar(gca,20,'mm')
    FintAng.Name = sprintf('%s%s_initialAngle_unconstrained',IT.subject,IT.date);

    %% Plot trajectories for different tube tube
    % Setup figure
    nRow = 1;
    nCol = 4;
    axSp = 5;
    axW = 300;
    [fW, fH, Ax1] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
    Ftube = figure('Position',[10 10 fW fH]);
    figName = sprintf('%s%s_IntTargConstAvg_%0.2d', IT.subject, IT.date, ...
        IT.intermediate_target_num);
    Ftube.Name = figName;
    titleStr = {'6C) Unconstrained trials','7B) Largest Tube',...
        '7B) Smallest Tube','7D) All Tubes (averages)'};
    n_tube = length(IT.constrainedTubeRadius);
    tube_radius = IT.constrainedTubeRadius;

    % 6A) Create figures for intermediate target - no tube
    plt.subplotSimple(nRow,nCol,1,'Ax',Ax1);
    TD(1).plot('ColMat',C,'plotTrajectories',0,'plotStartTarg',1,...
        'plotEndTarg',0,'plotTargNum',0)
    TD.plot('color',[0 0 0],'plotTargNum',0,'endMarker', [],'startMarker',[],...
        'trialsPerTarg',trialsPerCondition,'setRandSeed',setRandSeed)
    TDallAvg = TD.average('avgMode', avgMode);
    TDallAvg.plot('color',[0 0 0],'lineWidth',3,'plotTargets',0,'endMarker','arrow');

    % 7B) Iterate over presented tubes and plot successful trajectories
    TDavg = cell(1, n_tube);
    n_avg = nan(1, n_tube);
    for i = 1:n_tube
        % Average all successful trials
        TD = IT.TDconstrained{i};
        sc_mask = logical([TD.successful]);
        TDsuc = TD(sc_mask);
        TDfail = TD(~sc_mask);
        TDavg{i} = TDsuc.average('avgMode', avgMode);
        n_avg(1, i) = sum(sc_mask);

        % Plot the smallest tube trajectories
        if i == 1 || i == n_tube
            if i == 1
                plt.subplotSimple(nRow,nCol,2,'Ax',Ax1);
            else
                plt.subplotSimple(nRow,nCol,3,'Ax',Ax1);
            end
            % Plot start target
            TDavg{i}.plot('ColMat',C,'plotTrajectories',0, ...
                'plotStartTarg',1,'plotEndTarg',0,'plotTargNum', 0)
            % Plot end target
            TDavg{i}.plot('color',[0.3 0.3 0.3],...
                'plotTrajectories', 0,'plotTargNum', 0)

            % Plot individual trajectories
            if i == 1 % Plot the largest tube trajectories
                if ~isempty(trialsPerCondition)
                    if useExample
                        pos = [1 15 10 5 6 12 7 9 11 23];
                    elseif ~isempty(setRandSeed)
                        rng(setRandSeed)
                        pos = randperm(length(TD),trialsPerCondition);
                    else
                        pos = randperm(length(TD),trialsPerCondition);
                    end
                else
                    pos = 1:length(TD);
                end
                TD(pos).plot('color',[0 0 0],'plotTargets',0,'endMarker',[])
            else % Plot the smallest tube trajectories
                if ~isempty(trialsPerCondition)
                    if useExample
                        rng(3), posSuc = randperm(length(TDsuc));
                        posSuc = posSuc(1:trialsPerCondition);
                        posFail = [20 31 65 99 109 114 131 147 161 174];
                    elseif ~isempty(setRandSeed)
                        rng(setRandSeed)
                        posSuc = randperm(length(TDsuc),trialsPerCondition);
                        posFail = randperm(length(TDfail),trialsPerCondition);
                    else
                        posSuc = randperm(length(TDsuc),trialsPerCondition);
                        posFail = randperm(length(TDfail),trialsPerCondition);
                    end

                else
                    posSuc = 1:length(TDsuc);
                    posFail = 1:length(TDfail);
                end
                TDsuc(posSuc).plot('color',[0 0 0],'plotTargets',0, ...
                    'endMarker', [],'startMarker',[])
                TDfail(posFail).plot('color',[1 0 0],'plotTargets',0,...
                    'endMarker', [],'startMarker',[],'lineStyle','-')
            end
            % Plot average
            TDallAvg = TD.average('avgMode', avgMode);
            TDallAvg.plot('color',[0 0 0],'lineWidth',3,...
                'plotTargets',0,'endMarker','arrow');

            % Plot the tube
            IT.tubeObject(i).plot('color', [0 0 0],'lineSpec', '-');
        end
    end

    % 7C) Plot average trajectories for each tube size
    % Initialize figure setup including arrays for legend and tube colormap
    plt.subplotSimple(nRow,nCol,4,'Ax',Ax1);
    LH = nan(1, n_tube);
    leg_str = cell(1, n_tube);
    col = tube.get_tube_col(IT.startTargPos, n_tube,'colStart',[204 204 204]/255);

    % Iterate over tubes and plot
    for i = 1:length(TDavg)
        if i == 1
            % Plot targets
            TDavg{i}.plot('color', col.gray,'plotTrajectories',0, ...
                'plotIntTarg',0,'plotTargNum',0)
            TDavg{i}.plot('color',col.tube(1,:),'plotStartTarg',0, ...
                'plotEndTarg',0)
        end
        TDavg{i}.plot('color',col.tube(i+1,:),'plotTargets',0,'lineWidth',3)
        IT.tubeObject(i).plot('color', col.tube(i+1,:),'lineSpec', '-');

        % Get line object handle and legend text
        LH(i) = plot(nan, nan, '-', 'color', col.tube(i+1,:), 'lineWidth',3);
        leg_str{i} = sprintf('r = %d mm (%d trials)', tube_radius(i), ...
            n_avg(1, i));
    end

    % Format average each subpanel of trajectory figure
    for i = 1:4
        plt.subplotSimple(nRow,nCol,i,'Ax',Ax1)
        if ~isempty(plotScale)
            axis([-1 1 -1 1]*plotScale)
        end
        plt.scaleBar(gca,20,'mm')

        if i == 4 % Plot legend
            legend(LH, leg_str,'box','off','Location','SouthWest','FontSize',12)
        end

        axis on
        set(gca,'XTick',[],'YTick',[])
        title(titleStr{i}, 'FontSize', 16, 'FontWeight','Bold');
    end
    %% 7D) Calculate the angle difference for the first constrained trials
    TD = IT.TDconstrained{end};
    TD = TD(ismember([TD.startPos]',startPos,'rows'));

    [FintAngTube] = plot_intAngComparisons(TD,TD_tt_rot,startPos,C);
    TU = [TD.tube];
    TU(end).plot('lineSpec','-') % Plot the smallest tube

    title(sprintf('%s%s Average trajectories and initial angles \n Constrained',...
        IT.subject,IT.date))
    axis([-1 1 -1 1]*plotScale)
    axis on,set(gca,'XTick',[],'YTick',[])
    plt.scaleBar(gca,20,'mm')
    FintAngTube.Name = sprintf('%s%s_initialAngle_constrained',IT.subject,IT.date);

    %% 6F) The average trajectories of early vs late
    posSes = find(ismember({datNoChange.date},exampleSess(n)));

    % Plot the ortho/opposite traces
    TD = IT.TD_tt_rot{1};
    avgTD = TD.average();
    startPos = IT.TDconstrained{1}(1).startPos';
    pos = ismember([TD.startPos]',startPos,'rows');
    TD = TD(pos==1);

    avgTD_e = TD(datNoChange(posSes).idx_E).average('avgMode','samp');
    avgTD_l = TD(datNoChange(posSes).idx_L).average('avgMode','samp');
    F_EL = figure;
    avgTD.plot('ColMat',C,'plotTrajectories',0,'plotTargNum',0,...
        'plotEndTarg',0,'plotStartTarg',1)
    avgTD_e.plot('ColMat',C,'endMarker','arrow','plotTarget',0,'lineWidth',3)
    avgTD_l.plot('ColMat',C,'endMarker','arrow','plotTarget',0,...
        'lineWidth',3,'lineStyle','--')
    % Plot intermediate target
    IT.TDconstrained{end}(1).plot('color',[0 0 0],'plotTrajectories',0, ...
        'plotTargNum',0)
    axis([-1 1 -1 1]*plotScale)
    axis on,set(gca,'XTick',[],'YTick',[])
    plt.scaleBar(gca,20,'mm')
    title(sprintf('%s%s Average trajectories \n Early vs Late',...
        IT.subject,IT.date))
    F_EL.Name = sprintf('%s%s_Early_Vs_LateTraj',IT.subject,IT.date);
end

%% Save the figures
if saveFig
    if isempty(savePathBase)
        savePathBase = uigetdir;
    end

    saveFigurePDF([Fhist,FintAng,Ftube,FintAngTube,F_EL],savePathBase)
end
end

%% Internal function
function [F] = plot_intAngComparisons(TD,TD_tt_rot,startPos,C)

rotTD = TD_tt_rot.average;
rotTDplt = TD_tt_rot.average('avgMode','samp');

mask = ismembertol([rotTD.startPos]',startPos,1e-6,'ByRows',true);

% Calculate the flow and againstFlow trajectories using the time-warped
% average
flowTraj = diff(rotTD(mask).brainKin.pos([1 25],:))';
flowTraj = flowTraj./norm(flowTraj);

antFlowTraj = diff(rotTD(~mask).brainKin.pos([end end-19],:))';
antFlowTraj = antFlowTraj./norm(antFlowTraj);

% Calculate the average trajectory
TD = TD(ismember([TD.startPos]',startPos,'rows'));
aTD = TD.average('avgMode','samp');

F = figure; hold on
rotTDplt(mask).plot('trajColorWeight',.5,'ColMat',C,'endMarker','arrow',...
    'lineWidth',4,'plotTargNum',0,'plotStartTarg',1,'plotEndTarg',0);
aTD.plot('color',[.7 .7 .7],'plotTargNum',0,'endMarker','arrow','lineWidth',4)
line(aTD.startPos(1:2),aTD.targPos(1:2),'Color','k', ...
    'LineStyle','--','LineWidth',2)
end