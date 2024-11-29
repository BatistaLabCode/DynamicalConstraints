% [P] = opt.get_data(x, cond_cond) get_data  Get data for optimization
%
% This function is a standalone version of the 'calcOptimizationData'
% function used online to get the data used during optimization. This
% version was created in order to simplify the process of finding optimized
% projections offline without needing to use the full 'GPFA' data
% structure.
%
% Usage:
%   [P] = opt.get_data(x, cond_cond)
%
% Inputs:
%   x           Cell array of latent states for each trial (1 x n_trials)
%   cond_code   Code indicating the condition for each element of x
%
% Optional Inputs:
%   plotFigs    Determine if the data should be plotted
%   nAvg        Number of datapoints to average over
% 
% Outputs:
%   P           Structure of parameterized data passed to the optimization.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [P,D] = get_data(X, cond_code, varargin)

% Parse optional arguments
plotFigs = false;
nAvg = 1;               % Number of datapoints to average over

assignopts(who,varargin);

% Get dimensionality of inputs
n_trials = length(X);
x_dim = size(X{1}, 1);

% --- Calculate mean neural state at the start and end of trajectories ---
nAB = sum(cond_code == 1);
nBA = sum(cond_code == 2);
xAMat_start = nan(x_dim, nAvg*nAB);
xAMat_end = nan(x_dim, nAvg*nAB);
xBMat_start = nan(x_dim, nAvg*nBA);
xBMat_end = nan(x_dim, nAvg*nBA);

% Get data for start of A->B trajectory
X_AB = X(cond_code == 1);
onset = 1;
for i = 1:nAB
    offset = onset + nAvg - 1;
    xAMat_start(:, onset:offset) = X_AB{i}(:,1:nAvg);
    xAMat_end(:, onset:offset) = X_AB{i}(:, (end - nAvg + 1):end);
    onset = offset + 1;
end

% Get data for start of B-A trajectory
X_BA = X(cond_code == 2);
onset = 1;
for i = 1:nBA
    offset = onset + nAvg - 1;
    xBMat_start(:, onset:offset) = X_BA{i}(:, 1:nAvg);
    xBMat_end(:, onset:offset) = X_BA{i}(:, (end - nAvg + 1):end);
    onset = offset + 1;
end

% Calculate start and end points
xA_start = mean(xAMat_start,2,'omitnan');  % Start point for A->B
xA_end = mean(xAMat_end,2,'omitnan');  % End point for A->B
xB_start = mean(xBMat_start,2,'omitnan');  % Start point for B->A
xB_end = mean(xBMat_end,2,'omitnan');  % End point for B->A
xMid = (xA_start + xA_end)/2;  % Midpoint between start of A->B and B->A

% --- Calculate various metrics potentially used during optimization ---
% NOTE: some of these are legacy metrics that are no longer used but are
% being left here for consistency with previous versions of the function
% this is based on.

% Calculate projection vector
p = (xA_start - xB_start)/norm(xA_start - xB_start);

% Calculate the following:
% - Distance from the start of A->B
% - Distance from the start of B->A
% - Distance from the midpoint
% - Projection onto line between A and B

D = repmat(struct(), n_trials, 1);
nSamp = nan(n_trials, 1);
for i = 1:n_trials
    % Get neural trajectories
    x = X{i};
    nSamp(i) = size(x,2);
    
    % Calculate distance from A
    xD = x - xA_start*ones(1,nSamp(i));
    dA = nan(nSamp(i),1);
    for j = 1:nSamp(i)
        dA(j) = norm(xD(:,j));
    end
    D(i).dA = dA;
    
    % Calculate distance from B
    xD = x - xB_start*ones(1,nSamp(i));
    dB = nan(nSamp(i),1);
    for j = 1:nSamp(i)
        dB(j) = norm(xD(:,j));
    end
    D(i).dB = dB;
    
    % Calculate distance from midpoint
    xD = x - xMid*ones(1,nSamp(i));
    dM = nan(nSamp(i),1);
    for j = 1:nSamp(i)
        dM(j) = norm(xD(:,j));
    end
    D(i).dM = dM;
    
    % Calculate projection onto line between two start points
    dP = nan(nSamp(i),1);
    for j = 1:nSamp(i)
        dP(j) = x(:,j)' * p;
    end
    D(i).dP = dP;
    D(i).pStart = mean(dP(1:2),'omitnan');
    D(i).pEnd = mean(dP(end-1:end),'omitnan');
    
    % Get TC info
    D(i).tC = cond_code(i);
end

% Identify closest data point to projection

% Find center of projection
pStart = [D.pStart];
pEnd = [D.pEnd];
tcMask = cond_code == 1;
pAB_Start = mean(pStart(tcMask),'omitnan');
tcMask = cond_code == 2;
pBA_Start = mean(pStart(tcMask),'omitnan');
pCenter = mean([pAB_Start pBA_Start]);

% Find closest point to dP = 0
for i = 1:n_trials
    dP = D(i).dP - pCenter;
    [~,I] = min(abs(dP));
    D(i).midIdx = I;
end

% Calculate midpoints and covariance at midpoint
xDim = size(X{1}, 1);
x_mid = nan(xDim, n_trials);

% Iterate over all trials and get data at the midpoint
for i = 1:n_trials
    idx = D(i).midIdx;
    x_mid(:, i) = X{i}(:, idx);
end

mu_A = xA_start;
mu_A_end = xA_end;
mu_B = xB_start;
mu_B_end = xB_end;

% Get xAB (tC == 1)
tcMask = cond_code == 1;

X_temp = X(tcMask);
D_temp = D(tcMask);
nTrials_temp = sum(tcMask);
xMat = nan(xDim, nTrials_temp);

for i = 1:length(D_temp)
    idx = D_temp(i).midIdx;
    xMat(:,i) = X_temp{i}(:, idx);
end
mu_AB = mean(xMat,2);
sig_AB = cov(xMat');
xMat_AB = xMat;

% Get xBA (tC == 2)
tcMask = cond_code == 2;

X_temp = X(tcMask);
D_temp = D(tcMask);
nTrials_temp = sum(tcMask);
xMat = nan(xDim,nTrials_temp);

for i = 1:length(D_temp)
    idx = D_temp(i).midIdx;
    xMat(:,i) = X_temp{i}(:, idx);
end
mu_BA = mean(xMat,2);
sig_BA = cov(xMat');
xMat_BA = xMat;
% Define data structure used in optimization
P.D = xDim;
P.N = [];       % Not needed
P.mu_A = mu_A;
P.sig_A = [];   % Not needed
P.mu_B = mu_B;
P.sig_B = [];   % Not needed
P.mu_AB = mu_AB;
P.sig_AB = sig_AB;
P.mu_BA = mu_BA;
P.sig_BA = sig_BA;
P.x_mid = x_mid;
P.cond_code = cond_code;
P.xMat_AB = xMat_AB;
P.xMat_BA = xMat_BA;

% End points of trajectories
P.mu_A_end = mu_A_end;
P.mu_B_end = mu_B_end;

% --- Plot distance metric data if desired ---
if plotFigs
    % Setup figure
    nCol = 1;
    nRow = 4;
    axW = 800;
    axH = 175;
    axSp = 15;
    [fW,fH,Ax] = calcFigureSize(nRow,nCol,axW,axH,axSp);
    F1 = figure('Position',[100 100 fW fH]);

    col = [64 129 183;181 82 39]/255;
    plotIdx = 1;
    H = subplotSimple(nRow,nCol,plotIdx,'Ax',Ax); hold on;
    for i = 1:n_trials
        plot(D(i).dA,'color',col(D(i).tC,:))
    end

    % Plot lines for legend (these aren't shown)
    L(1) = plot([-1 -1],[0 1],'color',col(1,:));
    L(2) = plot([-1 -1],[0 1],'color',col(2,:));

    ylabel('Distance from $\mathbf{x}_A$','Interpreter','latex','FontSize',14)
    set(H.XAxis,'Visible','off');
    set(H,'TickDir','out')
    set(gca,'XLim',[1 max(nSamp)])
    legend(L,{'$A \rightarrow B$','$B \rightarrow A$'}, ...
        'Interpreter','latex','FontSize',14)

    plotIdx = plotIdx + 1;
    H = subplotSimple(nRow,nCol,plotIdx,'Ax',Ax); hold on;
    for i = 1:n_trials
        plot(D(i).dB,'color',col(D(i).tC,:))
    end
    set(gca,'XLim',[1 max(nSamp)])
    ylabel('Distance from $\mathbf{x}_B$','Interpreter','latex','FontSize',14)
    set(H.XAxis,'Visible','off');
    set(H,'TickDir','out')

    plotIdx = plotIdx + 1;
    H = subplotSimple(nRow,nCol,plotIdx,'Ax',Ax); hold on;
    for i = 1:n_trials
        plot(D(i).dM,'color',col(D(i).tC,:))
    end
    set(gca,'XLim',[1 max(nSamp)])
    ylabel('Distance from midpoint','Interpreter','latex','FontSize',14)
    set(H.XAxis,'Visible','off');
    set(H,'TickDir','out')

    plotIdx = plotIdx + 1;
    H = subplotSimple(nRow,nCol,plotIdx,'Ax',Ax); hold on;
    for i = 1:n_trials
        plot(D(i).dP,'color',col(D(i).tC,:))
    end
    set(gca,'XLim',[1 max(nSamp)])
    plot(get(gca,'XLim'),[0 0],'k--')
    ylabel('Projection onto $\mathbf{x}_A - \mathbf{x}_B$ axis', ...
        'Interpreter','latex','FontSize',14)
    set(H.XAxis,'Visible','off');
    set(H,'TickDir','out')

    plotTitle('Neural trajectory distance metrics')
end

end
