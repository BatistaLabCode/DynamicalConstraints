% calc_angle_diff      Calculates the difference of angle between two vectors
%
% [angle,dProd,baseAng,F,flowTraj,antFlowTraj,signAng] = calc_angle_diff(TD,baseTD,varargin)
%
% Usage:
%   util.calc_angle_diff(TD,baseTD)
% 
% Inputs:
%   baseDir         Base directory to identified all dependent content
%   itemName        File identifier to subselect for specific file content.              
% 
% Optional Inputs:
%   strtPt          Start index of the initial angle vector
%   endPt           End index of the initial angle vector
%   kinSrc          kinSource for the data
%   intTraj         Vector for the intermediate target
%   flowTraj        Vector for the with flow trajectory
%   antFlowTraj     Vector for the against flow trajectory
%   dim             Dimensions of the workspace
%   F               Figure handle
%   plotAng         Determine whether to plot angle comparison.
%   plotScale       Plotting scale
%   plotBase        Determine whether to plot lines for intTraj, flowTraj 
%                       and antFlowTraj
%   debugMode       Pause the run for debugging and sanity checks
%   offset          Apply an offset to the trajectory plotting
%   compInt         Calculate the angles using intTraj (default) or 
%                       flowTraj as base comparison
%   colMat          Define the colormap
%
% Outputs:
%   angle           Angle comparison to [flowTraj intTraj antFlowTraj]
%   dProd           Dot product between the trajectory vector and 
%                       [flowTraj intTraj antFlowTraj]
%   baseAng         Base angles of [flowTraj intTraj antFlowTraj]
%   F               Figure handle
%   flowTraj        Vector for the with flow trajectory
%   antFlowTraj     Vector for the against flow trajectory
%   signAng         Sign of the calculated angles
%
% Author:   Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [angle,dProd,baseAng,F,flowTraj,antFlowTraj,signAng] = calc_angle_diff(TD,baseTD,varargin)
%  Function calculates the difference of angle between two vectors
% Optional Inputs
strtPt = 1;         % Start index of the initial angle vector
endPt = 4;          % End index of the initial angle vector
kinSrc = 'brain';   % kinSource for the data
intTraj = [];       % Vector for the intermediate target
flowTraj = [];      % Vector for the with flow trajectory
antFlowTraj = [];   % Vector for the against flow trajectory
dim = 1:2;          % Dimensions of the workspace
F = [];             % Figure handle
plotAng = false;    % Determine whether to plot angle comparison.
plotScale = 1;      % Plotting scale
plotBase = true;    % Determine whether to plot lines for intTraj, flowTraj 
                    %   and antFlowTraj
debugMode = false;  % Pause the run for debugging and sanity checks
offset = [0 0];     % Apply an offset to the trajectory plotting
compInt = 0;        % Calculate the angles using intTraj (default) or 
                    %   flowTraj as base comparison
colMat = [];        % Define the colormap

assignopts(who, varargin);

% Define the ideal path trajectory
if isempty(intTraj)
    intTraj = TD(1).targPos(dim) - TD(1).startPos(dim);
end
intTraj = intTraj./sqrt(sum(intTraj.^2));

% Find the with-flow trajectory
mask = ismember([baseTD.startPos]',TD(1).startPos','rows');
if isempty(flowTraj)
    flowTraj = diff(baseTD(mask).([kinSrc 'Kin']).pos([strtPt endPt],dim))';
end
flowTraj = flowTraj./sqrt(sum(flowTraj.^2));

% Find the against-flow trajectory
if isempty(antFlowTraj)
    antFlowTraj = diff(baseTD(~mask).([kinSrc 'Kin']).pos(...
        [end-(strtPt-1) end-(endPt-1)],dim))';
end
antFlowTraj = antFlowTraj./sqrt(sum(antFlowTraj.^2));

% Create a figure for plotting
if plotAng & isempty(F)
    F = figure;
end

% Interate through aTD and calculate the angle difference
[angle,baseAng,dProd,signAng] = deal(nan(3,length(TD)));

for n = 1:length(TD)
    % Skip the trial if the endPt is larger than the size of the trajectory
    if endPt>size(TD(n).([kinSrc 'Kin']).pos,1)
        continue
    end

    % Define the initial vector of the comparison trajectory
    startVec = diff(TD(n).([kinSrc 'Kin']).pos([strtPt endPt],dim))';
    startVec = startVec./sqrt(sum(startVec.^2));

    % Calculate the dot product and angle between the vectors
    dProd(:,n) = startVec'*[flowTraj intTraj antFlowTraj];

    % Calculate the sign for the angle
    signAng(:,n) = [det([flowTraj startVec]) det([intTraj startVec]) det([antFlowTraj startVec])];
    signAng(signAng(:,n)<0,n) = -1;
    signAng(signAng(:,n)>=0,n) = 1;
    angle(:,n) = signAng(:,n).*acosd(dProd(:,n));

    if compInt
        %baseAng(:,n) = acosd(round(intTraj'*[flowTraj intTraj antFlowTraj],5));
        baseAng(:,n) = acosd(round(intTraj'*[flowTraj intTraj antFlowTraj],10));
        tmp = [det([flowTraj intTraj]) det([intTraj intTraj]) det([antFlowTraj intTraj])];
    else
        baseAng(:,n) = acosd(round(flowTraj'*[flowTraj intTraj antFlowTraj],10));
        tmp = [det([flowTraj flowTraj]) det([intTraj flowTraj]) det([antFlowTraj flowTraj])];
    end
    tmp = sign(tmp);
    tmp(tmp==0) = 1;

    % Flip the sign here because we want it centered around the comparison point.
    baseAng(:,n) = -tmp'.*baseAng(:,n);

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
        % Set the color
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

% Add the legend
if plotAng
    lgd = legend([fP iP aP traj],...
        {'With Flow','Direct Path','Against Flow','Trajectory Angle'});
    lgd.Location = 'eastoutside';
end
end