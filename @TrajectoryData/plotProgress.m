function [S] = plotProgress(TD)
% plotProgress
%
% Plot progress summary for TrajectoryData object
%
% This function will generate a performance progress summary plot for the
% input TrajectoryData object.
%
% Author:  Alan D. Degenhart
% Date Created: 2016/06/21
% Last Updated: 2016/06/21
% Last Update:  Initial version of function

winSize = 50; % Window size for calculating running average
showInitTrials = false;

% Get data size info and initialize matrices
nTrials = length(TD);
srSmoothed = nan(nTrials,1);    % Smoothed success rate
atSmoothed = nan(nTrials,1);    % Smoothed acquisition time

% Get success data.
sc = logical([TD.successful]);

% Get acquisition times.
trajOnset = [TD.trajOnset];
trajOffset = [TD.trajOffset];
at = trajOffset - trajOnset;

% Get flags indicating transitions between states.  Typically this will be
% the decoder number, but this function could be written to be generic
% enough to take arbitrary flags.
decNum = [TD.decoderNum];
decChgInd = find(diff(decNum));
decChgOnset = [1 (decChgInd + 1)];  % Add first index, account for diff
decChgOffset = [decChgInd nTrials]; % Add last index

ax(1) = subplotSimple(2,1,1,'xMarg',[.1 .1],'yMarg',[.1 .1]); hold on;
ax(2) = subplotSimple(2,1,2,'xMarg',[.1 .1],'yMarg',[.1 .1]); hold on;

% Plot block transitions
yLim = [.4 1.05;300 2500]; % Plot limits (currently hard-coded)

% Loop over blocks and smooth success rate data
nBlocks = length(decChgOnset);
decStr = cell(nBlocks,1);
for i = 1:nBlocks
   % Get data for current block
   blkInd = decChgOnset(i):decChgOffset(i);     % Indices of trials for block
   blkTrials = length(blkInd);                  % Number of trials in the current block
   scBlk = sc(blkInd);
   atBlk = at(blkInd);
   winBlk = min([blkTrials winSize]);
   
   % Set acquisition times for failed trials to NaN
   atBlk(scBlk == 0) = NaN;
   
   % Calculate running average using convolution
   srBlk = movmean(scBlk,[winSize-1 0],'omitnan');
   atBlk = movmean(atBlk,[winSize-1 0],'omitnan');
   
   % Plot success rate, acquisition time, and block boundaries
   invalidInd = 1:(winBlk-1);   % "Invalid" section of block window
   validInd = winBlk:blkTrials; % "Valid" section of block window
   axes(ax(1))
   plot(blkInd(validInd),srBlk(validInd),'k')
   axes(ax(2))
   plot(blkInd(validInd),atBlk(validInd),'k')
   
   % Plot invalid trials if desired
   if showInitTrials  
       axes(ax(1))
       plot(blkInd(invalidInd),srBlk(invalidInd),'r')
       axes(ax(2))
       plot(blkInd(invalidInd),atBlk(invalidInd),'r')
   else
       % Shade out region for first 50 trials
       pInd = blkInd(invalidInd);
       x = [pInd(1) pInd(end) pInd(end) pInd(1)];
       ySR = [yLim(1,1) yLim(1,1) yLim(1,2) yLim(1,2)];
       yAT = [yLim(2,1) yLim(2,1) yLim(2,2) yLim(2,2)];
       
       col = ones(1,3) * .9;
       axes(ax(1))
       patch(x,ySR,col,'EdgeColor','k')
       axes(ax(2))
       patch(x,yAT,col,'EdgeColor','k')
   end
   
   % Add smoothed data for the current block into the trial-wide vector
   srSmoothed(blkInd) = srBlk;
   atSmoothed(blkInd) = atBlk;
   
   % Set text string for current decoder
   decStr{i} = sprintf('Decoder %d',TD(blkInd(1)).decoderNum);
end

% Plot block transitions
%yLim = [0 1.05;300 2500]; % Plot limits
xLim = [0 nTrials];
x = repmat(decChgOnset - 0.5,2,1);

axes(ax(1))
y = repmat(yLim(1,:)',1,nBlocks);
plot(x,y,'k--')
text(x(1,:),y(1,:),decStr,'color','k','Rotation',90, ...
    'VerticalAlignment','top')
set(gca,'Box','on','YLim',yLim(1,:),'XLim',xLim,'XTickLabel',[],'yTick',yLim(1,1):.2:yLim(1,2))
ylabel('Success Rate')

axes(ax(2))
y = repmat(yLim(2,:)',1,nBlocks);
plot(x,y,'k--')
text(x(1,:),y(1,:),decStr,'color','k','Rotation',90, ...
    'VerticalAlignment','top')
set(gca,'Box','on','YLim',yLim(2,:),'XLim',xLim)
ylabel('Acquisition Time (ms)')
xlabel('Trial')

% Pack progress data into output structure
S.subject = TD(1).subject;
S.date = TD(1).date;
S.dataset = datestr(datenum(S.date),'yyyymmdd');
S.successCode = sc;
S.acquireTime = at;
S.successRate_smoothed = srSmoothed;
S.acquireTime_smoothed = atSmoothed;
S.nBlocks = nBlocks;
S.blockOnset = decChgOnset;
S.blockOffset = decChgOffset;
S.blockDecoder = decNum(decChgOnset);