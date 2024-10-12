function [F] = plot_success_progression(IT,varargin)
binSize = 50;
startPos = [];
inclWO = 0;
inclPred = 1;
C = util.defineTaskColormap('bc_rot');

% Parse optional agruments
assignopts(who,varargin);

if isempty(startPos)
    startPos = IT.TDconstrained{1}(1).startPos';
end

wo_pres = 0;

% Set the color maps
scale = norm(startPos)./100;
C.targPos = C.targPos.*(scale/.9);
h = tube.get_tube_col(startPos,size(IT.TDconstrained,1),'C',C);
    
%% Collect the success rates for each block
startPos = IT.TDconstrained{1}(1).startPos';

% Collect the unconstrained blocksm, but only for the tested start target
unconSuc = [IT.TDunconstrained.successful];
mask = ismember([IT.TDunconstrained.startPos]',startPos,'rows');
unconSuc = unconSuc(mask==1);

if ~isempty(IT.TDunconstrainedWashout) & inclWO
    unconWOSuc = [IT.TDunconstrainedWashout.successful];
    mask = ismember([IT.TDunconstrainedWashout.startPos]',startPos,'rows');
    unconWOSuc = unconWOSuc(mask==1);
    wo_pres = 1;
end

% Collect the constrained blocks
for n = 1:size(IT.TDconstrained,1)
    con{n} = [IT.TDconstrained{n}.successful];
end

% Combine the session
if wo_pres
    sesDat = [unconSuc con unconWOSuc];
else
    sesDat = [unconSuc con];
end

% Calculate the moving average for each block
for n = 1:size(sesDat,2)
    temp = conv(sesDat{n},ones(1,binSize));
    temp = temp(1:size(sesDat{n},2));
    
    % Adjust the convolution for the first binSize number of trials
    if size(temp,2) <= binSize;
        sesDat{n} = 100*temp./[1:size(temp,2)];
    else
        rL = size(temp,2) - binSize; % remaining length to account for
        sesDat{n} = 100*temp./[1:binSize binSize*ones(1,rL)];
    end
    
    blockL(n,:) = size(sesDat{n},2);
end

%% Plot the success rate trend
j = 1;
allData = [sesDat{:}];

F = figure, hold on
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
        sPos = sum(blockL([1:n-1]))+1;
        endPos = sum(blockL([1:n]));
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
    if n >1 & blockL(n)<binSize;
        pos = blockL(n);
    elseif n>1
        pos = blockL(n):-50:50;
        pos = fliplr(pos);
    else 
        pos = [];
    end
    
    % Plot the success threshold crossings
    if sum(sesDat{n}(pos)>=75)
        tPos = pos(sesDat{n}(pos)>=75);
        leg1 = plot(tXval(tPos),sesDat{n}(tPos),'go','MarkerFaceColor','g');
    else
        tPos = pos(sesDat{n}(pos)<75);
        leg2 = plot(tXval(tPos),sesDat{n}(tPos),'ro','MarkerFaceColor','r');
    end
    
    % Plot the predicted success rate
    if inclPred
        if n == 1
            predSuc = sum(unconSuc)./length(unconSuc);
        elseif ~inclWO
            predTrl = IT.TD_unconst_tube{n-1};
            predSuc = sum([predTrl.successful]==1)./length(predTrl)
        end
        line([sPos endPos],100*[predSuc predSuc],'Color','r','LineStyle','--');
    end
end

legend([leg1,leg2,r,thres],{'Succeeded threshold','Failed threshold',...
    'Flux success rate','75% Threshold'})
xlim([0 endPos])
xlabel('Trial Count')
ylabel('Success rate (%)')
plt.plotTitle(sprintf('%s%s success progression constrained blocks',...
    IT.subject,IT.date));
F.Name = sprintf('%s%s_successProgression_constrainedBlocks',...
    IT.subject,IT.date);

end