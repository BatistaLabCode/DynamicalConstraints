function [angle,dProd,baseAng,F,flowTraj,antFlowTraj,signAng] = calc_angle_diff(TD,baseTD,varargin)
%%  Function calculate the difference of angle between two vectors
srtPt = 1;
endPt = 4;
kinSrc = 'brain';
intTraj = [];
flowTraj = [];
antFlowTraj = [];
dim = [1:2];
F = [];
plotAng = false;
plotScale = 1;
plotBase = true;
debugMode = false;
offset = [0 0];
compInt = 0;
colMat = [];%jet(1:size(TD,1));

extraOpts = assignopts(who, varargin);

% Define the ideal path trajectory
if isempty(intTraj)
    intTraj = TD(1).targPos(dim) - TD(1).startPos(dim);
end
intTraj = intTraj./sqrt(sum(intTraj.^2));

% Find the with-flow trajectory
mask = ismember([baseTD.startPos]',TD(1).startPos','rows');
if isempty(flowTraj)
    flowTraj = diff(baseTD(mask).([kinSrc 'Kin']).pos([srtPt endPt],dim))';
end
flowTraj = flowTraj./sqrt(sum(flowTraj.^2));

% Find the against-flow trajectory
if isempty(antFlowTraj)
    antFlowTraj = diff(baseTD(~mask).([kinSrc 'Kin']).pos(...
        [end-(srtPt-1) end-(endPt-1)],dim))';
end
antFlowTraj = antFlowTraj./sqrt(sum(antFlowTraj.^2));

% Interate through aTD and calculate the angle difference
angle = nan(3,length(TD));

if plotAng & isempty(F) 
    F = figure;
end

for n = 1:length(TD)
    % Skip the trial if the endPt is larger than the size of the trajectory
    if endPt>size(TD(n).([kinSrc 'Kin']).pos,1)
        continue
    end
    
    % Define the initial vector of the comparison trajectory
    startVec = diff(TD(n).([kinSrc 'Kin']).pos([srtPt endPt],dim))';
    startVec = startVec./sqrt(sum(startVec.^2));
    
    % Calculate the dot product and angle between the vectors
    dProd(:,n) = startVec'*[flowTraj intTraj antFlowTraj];
    signAng(:,n) = [det([flowTraj startVec]) det([intTraj startVec]) det([antFlowTraj startVec])];
    signAng(signAng(:,n)<0,n) = -1;
    signAng(signAng(:,n)>=0,n) = 1;
    angle(:,n) = signAng(:,n).*acosd(dProd(:,n));
    %angle(:,n) = acosd(dProd(:,n));
    if compInt
        %baseAng(:,n) = acosd(round(intTraj'*[flowTraj intTraj antFlowTraj],5));
        baseAng(:,n) = acosd(round(intTraj'*[flowTraj intTraj antFlowTraj],10));
        tmp = [det([flowTraj intTraj]) det([intTraj intTraj]) det([antFlowTraj intTraj])];
    else
        baseAng(:,n) = acosd(round(flowTraj'*[flowTraj intTraj antFlowTraj],5));
        tmp = [det([flowTraj flowTraj]) det([intTraj flowTraj]) det([antFlowTraj flowTraj])];
    end
    tmp = sign(tmp);
    tmp(tmp==0) = 1;
    baseAng(:,n) = -tmp'.*baseAng(:,n); % Flipping here because we want it centered around the comparison point.
    % Sanity check the angles
    if plotAng
        figure(F), hold on
        if plotBase
            fP = quiver(offset(1),offset(2),plotScale*flowTraj(1),...
                plotScale*flowTraj(2),'r--','LineWidth',1);
            aP = quiver(offset(1),offset(2),plotScale*antFlowTraj(1),...
                plotScale*antFlowTraj(2),'c--','LineWidth',1);
            iP = quiver(offset(1),offset(2),plotScale*intTraj(1),...
                plotScale*intTraj(2),'--','LineWidth',1,'Color',[.7 .7 .7]);
            
        end
        if ~isempty(colMat)
            traj = quiver(offset(1),offset(2),plotScale*startVec(1),...
            plotScale*startVec(2),'Color',colMat(n,:),'LineWidth',2);
        else
        traj = quiver(offset(1),offset(2),plotScale*startVec(1),...
            plotScale*startVec(2),'k','LineWidth',2);
        end
        if debugMode 
            pause
            cla
            plotBase = true;
        else
            plotBase = false;
        end
    end
end
if plotAng
lgd = legend([fP iP aP traj],...
    {'With Flow','Direct Path','Against Flow','Trajectory Angle'});
lgd.Location = 'eastoutside';
end
end