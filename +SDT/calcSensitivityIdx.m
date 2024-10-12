function [dP,datROC] = calcSensitivityIdx(xA,xB,varargin)
% Calculate the sensitivity index for your data. This code currently
% assumes that the data is 1-dimensional for all comparisions. 
%
% Inputs: vectors xA and xB, the two conditions that are being compared. 
%
% Outputs: 
%       - dP a single scalar value. 
%       - datROC a structure that returns the data necessary to remake an
%           ROC curve for the data, ie P, N, sensitivity (TPR = TP/P),
%           specificity (TN/N = 1-FPR), criterion, and t-test results
%           comparing the two distributions.
%       * Note: A will always be defined as the positive trials and B will
%           always be the negative trials.
%
% Optional Inputs: 
%       - muA: Mean of A, if empty calculated from the given data
%       - muB: Mean of B, if empty calculated from the given data
%       - SIGMA: standard deviation of data, if empty calculated as the 
%               average of the std for A and the std for B.
%       - sig2A: variance of A, if empty calculated
%       - sig2B: variance of B, if empty calculated
%       - crit: Criterion to compare TPR and FPR over, if empty will simply
%                  test over the range of unique values for A and B
%                  combined.
%       - plotFig: If the ROC should plotted
%       - Ax:       Axis to plot the data, if empty a new figure is created
%
% Author:       Erinn Grigsby
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional inputs
muA = []; % Default empty
muB = []; % Default empty
meanSig = []; % Default empty
sig2A = []; % Default empty
sig2B = []; % Default empty
crit = []; % Default empty
plotFig = false; % Default false
Ax = []; % Default empty

optArg = assignopts(who,varargin);

% Calculate the SIGMA value
if isempty(meanSig)
    % Check if sigA is empty, otherwise calculate
    if isempty(sig2A)
        sig2A = var(xA);
    end
    % Check if sigB is empty, otherwise calculate
    if isempty(sig2B)
        sig2B = var(xB);
    end
    
    SIGMA = sqrt(mean([sig2A sig2B]));
else
    SIGMA = meanSig;
end

% Calculate the mean values
if isempty(muA)
    muA = mean(xA);
end
if isempty(muB)
    muB = mean(xB);
end

% Calculate d'
dP = (muA - muB)./SIGMA;

% Calculate the ROC curve data
if isempty(crit)
    crit(1,:) = unique([xA; xB]);
end

P = length(xA);
N = length(xB);
sensitivity = sum(xA>=crit)./P;
specificity = sum(xB<crit)./N;

% Run the unpair ttest on the distributions.
[t_test.reject,t_test.P,t_test.CI,t_test.STATS] = ttest2(xA,xB);

% Save the ROC data
datROC.dPrime = dP;
datROC.xA = xA;
datROC.P = P;
datROC.TPR = sensitivity;
datROC.xB = xB;
datROC.N = N;
datROC.TNR = specificity;
datROC.criterion = crit;
datROC.auc = (trapz(flipud(1-specificity'),flipud(sensitivity'))-0.5)./0.5;

% Save the t-test in the structure
datROC.t_test = t_test; 

% Plot the ROC 
if plotFig
    if isempty(Ax)
        F = figure;
    else
        axes(Ax)
        F = gcf;
    end
    
    hold on
    plot(1-specificity, sensitivity, optArg{:})
    line([0 1],[0 1],'LineStyle','--')
    axis square
    title('ROC')
    ylabel('Sensitivity'),xlabel('1-Specificity')
end
end