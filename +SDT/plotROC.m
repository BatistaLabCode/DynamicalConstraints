function [F] = plotROC(datROC,varargin)
% Plots the ROC curve data

% Optional Inputs
plotHist = 1; %Default plot
plotROC = 1; %Default plot
plotSS = 1; %Default plot
AxHist = []; %Default empty
AxROC = []; %Default empty
AxSS = []; %Default empty
F = [];
nBins = 25;

optArg = assignopts(who,varargin);

% Check if any axis are defined, if all three are not defined then create a
% shared figure
if isempty(AxHist) | isempty(AxROC) | isempty(AxSS)
    %Determine how many subplots to make
    subplt = sum(plotHist+plotROC+plotSS);
    
    F = figure('Position',[25 15 subplt*350 350]);
    
    % Plot the data
    validPlts = find([plotHist plotSS plotROC]);
    for n = 1:subplt
        subplot(1,subplt,n), hold on
        if validPlts(n) == 1
            [~,bins] = histcounts([datROC.xA; datROC.xB],nBins);
            histogram(datROC.xA,bins,'Normalization','probability')
            histogram(datROC.xB,bins,'Normalization','probability')
            xlabel('Distributions'),ylabel('Probability')
            legend('xA','xB')
            axis square
        elseif validPlts(n) == 2
            plot(datROC.criterion,datROC.TPR)
            plot(datROC.criterion,datROC.TNR)
            xlabel('Criterion'),ylabel('Probability')
            legend('Specificity','Sensitivity')
            axis square
        elseif validPlts(n) == 3
            plot(1-datROC.TNR,datROC.TPR)
            line([0 1],[0 1],'LineStyle','--')
            axis square
            title('ROC')
            ylabel('TPR'),xlabel('FPR')
        end
    end
end

if ~isempty(AxHist)
    axes(AxHist)
    [~,bins] = histcounts([datROC.xA; datROC.xB],nBins);
    hold on
    histogram(datROC.xA,bins,'Normalization','probability')
    histogram(datROC.xB,bins,'Normalization','probability')
    xlabel('Distributions'),ylabel('Probability')
    legend('xA','xB')
    axis square
end

if ~isempty(AxSS)
    axes(AxSS),hold on
    plot(datROC.criterion,datROC.TPR)
    plot(datROC.criterion,datROC.TNR)
    xlabel('Criterion'),ylabel('Probability')
    legend('Specificity','Sensitivity')
    axis square
end

if ~isempty(AxROC)
    axes(AxROC),hold on
    plot(1-datROC.TNR,datROC.TPR)
    line([0 1],[0 1],'LineStyle','--')
    axis square
    title('ROC')
    ylabel('TPR'),xlabel('FPR')
end
end