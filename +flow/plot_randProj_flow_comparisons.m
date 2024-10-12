function [F,FlowDataRaw,FlowData] = plot_randProj_flow_comparisons(varargin)
%% Load the flow data data information
filePath = [];
savePath = [];
exampleSes = {'20190719'};
condition = {'ii','iv','ia','ir','rr','rv','ra'};

assignopts(who, varargin);

if isempty(filePath)
    filePath = uigetdir('','Load flow field main folder');
end
if isempty(savePath)
    savePath = uigetdir('','Determine where to save the figures');
end
D = dir(fullfile(filePath,'*flowRandProject_results.mat'));

FlowDataRaw = repmat(struct('subject',[],'dataset',[],'targPair',[],...
    'angMean',[],'angMedian',[],'magMean',[],'magMedian',[],...
    'mseMean',[],'mseMedian',[],'n_valid',[]),length(D),1);

FlowData = repmat(struct('subject',[],'dataset',[],'targPair',[],...
    'angMean',[],'angMedian',[],'magMean',[],'magMedian',[],...
    'mseMean',[],'mseMedian',[],'n_valid',[]),length(D),1);

measVal = {'angMean','angMedian','magMean','magMedian','mseMean',...
    'mseMedian','n_valid'};
%%
for n = 1:size(D,1)
    tFF = load(fullfile(D(n).folder,D(n).name));
    tFF = tFF.FF;
    
    % Save subject and session information
    FlowDataRaw(n).subject = tFF(1).angMean.subject;
    FlowDataRaw(n).dataset = tFF(1).angMean.dataset;
    %FlowDataRaw(n).targPair = uniTarg;
    FlowData(n).subject = tFF(1).angMean.subject;
    FlowData(n).dataset = tFF(1).angMean.dataset;
%    FlowData(n).targPair = uniTarg;
    
    % Iterate through the measurements
    for k = 1:length(measVal)
        tDat = [tFF.(measVal{k})];
        
        % Identify the normalization values
        if ismember(measVal{k},'n_valid') % Use m_same and m_alt
            m_top = mean([tDat.ia]); % m_alt
            m_bot = mean([tDat.ii]); % m_same
        else
            m_top = mean([tDat.ii]); % m_same
            m_bot = mean([tDat.iv]); % m_rev
        end
        
        % Identify the mean int vs int data for each condition
        t_ii = [tDat.ii];
        m_ii(1) = mean(t_ii);%mean(t_ii(1:2:end));
        m_ii(2) = mean(t_ii);%mean(t_ii(2:2:end));
        
        % Iterate through the conditions
        for m = 1:length(condition)
            tMeas = [tDat.(condition{m})];
            tVal = (tMeas - m_top)./(m_bot - m_top);
            tNorm = tMeas;
            tNorm(1:2:end) = tNorm(1:2:end)./m_ii(1);
            tNorm(2:2:end) = tNorm(2:2:end)./m_ii(2);
            
            % Save the data
            FlowDataRaw(n).(measVal{k}).(condition{m}) = mean(tMeas); 
            FlowDataRaw(n).(measVal{k}).([condition{m} 'all']) = tMeas;
            FlowData(n).(measVal{k}).(condition{m}) = mean(tVal);
            FlowData(n).([measVal{k} '_new']).(condition{m}) = mean(tNorm);
        end
    end  
    
end

% Identify the example sessions locations
mask = ismember([{FlowData.dataset}],exampleSes);

%% Set up the figure
nCol = 4;
nRow = 3;
axSp = 65;
axW = 350;
[fW, fH, Ax] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
F(1) = figure('Position',[10 10 fW fH]);
F(1).Name = 'FlowField_summary_allAnimals';

% Create the colormaps
cII = [64,64,64]/255;
cIR = [186,186,186]/255;
cIV = [202,0,32]/255;
cIA = [244,165,130]/255;

% Create the animal masks and symbols
maskE = ismember([{FlowData.subject}],'Earl');
maskD = ismember([{FlowData.subject}],'Dwight');
maskQ = ismember([{FlowData.subject}],'Quincy');

%% Plot the voxel overlap data
OL = [FlowData.n_valid];
OLraw = [FlowDataRaw.n_valid];
ii = [OL.ii];
ir = [OL.ir];
iv = [OL.iv];
ia = [OL.ia];

% Plot the histograms for voxel overlap
dt = .025;
[~,bins] = histcounts([ii ir iv ia],-0.5+.5*dt:dt:1.5);

plt.subplotSimple(nRow,nCol,1,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized flow field overlap')
ylabel('Num. Experiments')
axis square

% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

%Plot the scatter data
cond = {'ii','iv','ia'};
irr = [OLraw.ir];
olLim = [0 100];
for n = 1:3
    tScat = [OLraw.(cond{n})];
    plt.subplotSimple(nRow,nCol,1+n,'Ax',Ax), hold on
    line(olLim,olLim,'LineStyle','--')
    plot(irr(maskE),tScat(maskE),'ko')
    plot(irr(maskD),tScat(maskD),'kd')
    plot(irr(maskQ),tScat(maskQ),'ks')
    xlabel('int vs rot'),ylabel(cond{n})
    axis square
    
    % Add the example session mark
    plot(irr(mask==1),tScat(mask==1),'o','MarkerFaceColor','r','MarkerSize',8)
end

%% Plot the histograms for flow field difference -- Match with what we did 
% for the original figure 4 (ie the python code)
%(ismember([{FlowData.subject}],'Quincy')==1)
flowDiff = [FlowData.magMedian];
flowDiffraw = [FlowDataRaw.magMedian];
ii = [flowDiff.ii];
ir = [flowDiff.ir];
iv = [flowDiff.iv];
ia = [flowDiff.ia];

dt = 0.25;
[~,bins] = histcounts([ii ir iv ia],-2.5+.5*dt:dt:5);%'binWidth',.2);
plt.subplotSimple(nRow,nCol,1+nCol,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized flow field difference')
ylabel('Num. Experiments')
axis square

% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

%Plot the scatter data
cond = {'ii','iv','ia'};
irr = [flowDiffraw.ir];
fdLim = [0 2];
for n = 1:3
    tScat = [flowDiffraw.(cond{n})];
    plt.subplotSimple(nRow,nCol,1+n+nCol,'Ax',Ax), hold on
    line(fdLim,fdLim,'LineStyle','--')
    plot(irr(maskE),tScat(maskE),'ko')
    plot(irr(maskD),tScat(maskD),'kd')
    plot(irr(maskQ),tScat(maskQ),'ks')
    xlabel('int vs rot'),ylabel(cond{n})
    axis square
    
    % Add the example session mark
    plot(irr(mask==1),tScat(mask==1),'o','MarkerFaceColor','r','MarkerSize',8)
end

%% Plot the MSE data
OL = [FlowData.mseMean];
OLraw = [FlowDataRaw.mseMean];
ii = [OL.ii];
ir = [OL.ir];
iv = [OL.iv];
ia = [OL.ia];

% Plot the histograms for voxel overlap
dt = 0.05
[~,bins] = histcounts([ii ir iv ia],-0.5+.5*dt:dt:1.5);

plt.subplotSimple(nRow,nCol,1+2*nCol,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized MSE')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

%Plot the scatter data
cond = {'ii','iv','ia'};
irr = [OLraw.ir];
olLim = [0 20];
for n = 1:3
    tScat = [OLraw.(cond{n})];
    plt.subplotSimple(nRow,nCol,1+n+2*nCol,'Ax',Ax), hold on
    line(olLim,olLim,'LineStyle','--')
    plot(irr(maskE),tScat(maskE),'ko')
    plot(irr(maskD),tScat(maskD),'kd')
    plot(irr(maskQ),tScat(maskQ),'ks')
    xlabel('int vs rot'),ylabel(cond{n})
    axis square
    
    % Add the example session mark
    plot(irr(mask==1),tScat(mask==1),'o','MarkerFaceColor','r','MarkerSize',8)
end

%% Plot an alternative metric for the flow field difference - Normalized

% Set up the figure
nCol = 4;
nRow = 1;
axSp = 65;
axW = 350;
[fW, fH, Ax] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
F(2) = figure('Position',[10 10 fW fH]);
F(2).Name = 'FlowField_summary_allAnimals_new';

% Overlap voxel
ffDiff = [FlowData.n_valid_new];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,1,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
gx = gca; gx.XLim(1) = 0
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized overlap difference')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

% Flow field magnitude difference
ffDiff = [FlowData.magMedian_new];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,2,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized flow field difference')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

% Flow field angle difference
ffDiff = [FlowData.angMean_new];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,3,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized flow field angle difference')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

% Flow field magnitude difference
ffDiff = [FlowData.mseMean_new];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,4,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Normalized MSE')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'v','Color',cII)
plot(ia(mask==1),0,'v','Color',cIA)
plot(iv(mask==1),0,'v','Color',cIV)
plot(ir(mask==1),0,'v','Color',cIR)

%% Plot an alternative metric for the flow field difference - Raw data

% Set up the figure
nCol = 4;
nRow = 1;
axSp = 65;
axW = 350;
[fW, fH, Ax] = plt.calcFigureSize(nRow,nCol,axW,axW,axSp);
F(3) = figure('Position',[10 10 fW fH]);
F(3).Name = 'FlowField_summary_allAnimals_raw';

% Overlap voxel
ffDiff = [FlowDataRaw.n_valid];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,1,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Overlap difference')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'o','Color',cII)
plot(ia(mask==1),0,'o','Color',cIA)
plot(iv(mask==1),0,'o','Color',cIV)
plot(ir(mask==1),0,'o','Color',cIR)

% Flow field magnitude difference
ffDiff = [FlowDataRaw.magMedian];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,2,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Flow field difference')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'o','Color',cII)
plot(ia(mask==1),0,'o','Color',cIA)
plot(iv(mask==1),0,'o','Color',cIV)
plot(ir(mask==1),0,'o','Color',cIR)

% Flow field angle difference
ffDiff = [FlowDataRaw.angMean];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,3,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('Flow field angle difference')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'o','Color',cII)
plot(ia(mask==1),0,'o','Color',cIA)
plot(iv(mask==1),0,'o','Color',cIV)
plot(ir(mask==1),0,'o','Color',cIR)

% Flow field magnitude difference
ffDiff = [FlowDataRaw.mseMean];
ii = [ffDiff.ii];
ir = [ffDiff.ir];
iv = [ffDiff.iv];
ia = [ffDiff.ia];

[~,bins] = histcounts([ii ir iv ia],50);%'binWidth',.2);
bins = bins + 0.5*mean(diff(bins));
plt.subplotSimple(nRow,nCol,4,'Ax',Ax), hold on
histogram(ii,bins,'EdgeColor',cII,'DisplayStyle','stairs')
histogram(ir,bins,'EdgeColor',cIR,'DisplayStyle','stairs')
histogram(iv,bins,'EdgeColor',cIV,'DisplayStyle','stairs')
histogram(ia,bins,'EdgeColor',cIA,'DisplayStyle','stairs')
%xlim([-0.5 1.5])
legend('int vs int','int vs rot','int vs int reverse','int vs int alt')
xlabel('MSE')
ylabel('Num. Experiments')
axis square
% Add the example session mark
plot(ii(mask==1),0,'o','Color',cII)
plot(ia(mask==1),0,'o','Color',cIA)
plot(iv(mask==1),0,'o','Color',cIV)
plot(ir(mask==1),0,'o','Color',cIR)

%% Save the figure and the data
%save(fullfile(filePath,'FlowComparison'),'FlowData','FlowDataRaw')
%saveFigurePDF(F,savePath)
end