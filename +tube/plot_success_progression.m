% plot_success_progression  Plots the success progression of a single
% session of intermediate target task
%
% Usage:
%   tube.plot_success_progression(IT)
%
% 
% This will automatically collect and fill in the intermediate target 
% experiment data and success rates based on the filepath locations. 
%
% Input:
%   IT          IntTargExp object
%
% Optional Input:
%   binSize     Largest number of trials between sucess rate evaluation
%   startPos    The start target of the experiment
%   inclPred    Plot the success rate of the tubes applied to the
%                   unconstrained trials.
%   C           Default ColMat structure
%
% Output:
%   IT          IntTargExp object
%
% Created by: Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [F] = plot_success_progression(IT,varargin)
binSize = 50;
startPos = [];
inclPred = 1;
C = util.defineTaskColormap('bc_rot');

% Parse optional agruments
assignopts(who,varargin);

if isempty(startPos)
    startPos = IT.TDconstrained{1}(1).startPos';
end

% Set the color maps
scale = norm(startPos)./100;
C.targPos = C.targPos.*(scale/.9);
h = tube.get_tube_col(startPos,size(IT.TDconstrained,1),'C',C);
    
% Collect the success rates for each block
startPos = IT.TDconstrained{1}(1).startPos';

% Collect the unconstrained blocks, but only for the tested start target
unconSuc = [IT.TDunconstrained.successful];
mask = ismember([IT.TDunconstrained.startPos]',startPos,'rows');
unconSuc = unconSuc(mask==1);

% Collect the constrained blocks
for n = 1:size(IT.TDconstrained,1)
    con{n} = [IT.TDconstrained{n}.successful];
end

% Combine the session
sesDat = [unconSuc con];

% Calculate the moving average for each block
for n = 1:size(sesDat,2)
    temp = conv(sesDat{n},ones(1,binSize));
    temp = temp(1:size(sesDat{n},2));
    
    % Adjust the convolution for the first binSize number of trials
    if size(temp,2) <= binSize
        sesDat{n} = 100*temp./(1:size(temp,2));
    else
        rL = size(temp,2) - binSize; % remaining length to account for
        sesDat{n} = 100*temp./[1:binSize binSize*ones(1,rL)];
    end
    
    blockL(n,:) = size(sesDat{n},2);
end

% Plot the success rate trend
allData = [sesDat{:}];

F = figure; hold on
ylim([0 100])
thres = line([0 size([sesDat{:}],2)],[75 75],'Color','k','LineStyle','--');

for n = 1:size(sesDat,2)
    % Set the x-axis positions and the color
    if n == 1
        sPos = 1;
        endPos = blockL(1);
        col = 'k';
        lineStyle = ':';
    else
        sPos = sum(blockL(1:n-1))+1;
        endPos = sum(blockL(1:n));
        col = h.tube(end-(n-1),:);
        lineStyle = '-';
    end
    
    tXval = sPos:endPos;
    
    % Add the block separation and the uncertainty shading
    if blockL(n)> binSize
        rectangle('Position',[sPos 0 binSize 100],...
            'FaceColor',[.9 .9 .9],'EdgeColor','none');
    else
        rectangle('Position',[sPos 0 blockL(n) 100],...
            'FaceColor',[.9 .9 .9],'EdgeColor','none');
    end
    r = plot(nan,nan,'s','MarkerSize',15,...
        'MarkerFaceColor',[.9 .9 .9],'Color',[.9 .9 .9]);
    line([endPos endPos],[0 100],'Color',.6*[1 1 1],'LineStyle','--')
    
    plot(tXval,sesDat{n},'Color',col,'LineStyle',lineStyle)
    
    % Find the threshold crossings
    if n >1 && blockL(n)<binSize
        pos = blockL(n);
    elseif n>1 && n == size(sesDat,2)
        pos = 50:50:blockL(n);
    elseif n>1 && n ~= size(sesDat,2)
        pos = blockL(n):-50:50;
        pos = fliplr(pos);
    else 
        pos = [];
    end

    % Add a point at 25th trial, if this is first time that the block
    % length is greater or equal to 50.
    if n == (find(blockL(2:end)>=binSize,1,'first')+1)
        pos = [25 pos];
    end
    
    % Plot the success threshold crossings
    mask = sesDat{n}(pos)>=75;
    if sum(mask==1)>0
        leg1 = plot(tXval(pos(mask==1)),sesDat{n}(pos(mask==1)),'go','MarkerFaceColor','g');
    end
    if sum(mask==0)>0
        leg2 = plot(tXval(pos(mask==0)),sesDat{n}(pos(mask==0)),'ro','MarkerFaceColor','r');
    end
    
    % Plot the predicted success rate
    if inclPred
        if n == 1
            predSuc = sum(unconSuc)./length(unconSuc);
        else
            predSuc = IT.predictedSuccessRate(n-1);
        end
        predR = line([sPos endPos],100*[predSuc predSuc],'Color','r','LineStyle','--');
    end
end

if inclPred
    legend([leg1,leg2,r,thres,predR],{'Succeeded threshold',...
        'Failed threshold','Flux success rate','75% Threshold',['Predicted' ...
        'success rate']})
else
    legend([leg1,leg2,r,thres],{'Succeeded threshold','Failed threshold',...
        'Flux success rate','75% Threshold'})
end
xlim([0 endPos])
xlabel('Trial Count')
ylabel('Success rate (%)')
plt.plotTitle(sprintf('%s%s success progression constrained blocks',...
    IT.subject,IT.date));
F.Name = sprintf('%s%s_successProgression_constrainedBlocks',...
    IT.subject,IT.date);

end