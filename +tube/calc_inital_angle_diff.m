% calc_inital_angle_diff       Creates the initial angle structure for the
% unconstrained and constrained conditions of a session.
%
% [AD,F] = calc_inital_angle_diff(IT)
%
% This function returns the valid items in the directory specified by
% BASEDIR with the name specified by ITEMNAME.
%
% Usage:
%   tube.calc_inital_angle_diff(IT)
% 
% Inputs:
%   baseDir         Base directory to identified all dependent content
%   itemName        File identifier to subselect for specific file content.              
%
% Optional Inputs:
%   baseDir         Base directory to identified all dependent content
%   itemName        File identifier to subselect for specific file content.
% 
% Outputs:
%   dList           List of files with valid elements.
%
% Author:   Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com
function [AD,F] = calc_inital_angle_diff(IT,varargin)
successOnly = true;
F = [];
plotFig = 1;
useRot = 1;
decoderVisual = 'rot';
extraOpts = assignopts(who, varargin);

% Set up the structure
AD = struct('subject',IT.subject,'dataset',IT.date,'compAng',[],...
    'avgUncon',[],'avgCon',[],'uncon',[],'con',[],'rotTT',[],...
    'tubes',IT.constrainedTubeRadius);

startPos = IT.TDconstrained{1}(1).startPos';
if strcmp(decoderVisual,'int')
    C = util.defineTaskColormap('bc_int');
else
    C = util.defineTaskColormap('bc_rot');
end
scale = norm(startPos(1:2))./100;
C.targPos = C.targPos.*(scale/.9);

%Determine the which rotated condition is in the direction of the
%constrained data
if useRot
    rotTD = IT.TD_tt_rot{1}.average;
    mask = ismembertol([rotTD.startPos]',startPos,1e-6,'ByRows',true);
    
    % Calculate the flow and againstFlow trajectories using the time-warped
    % average
    flowTraj = diff(rotTD(mask).brainKin.pos([1 25],:))';
    flowTraj = flowTraj./norm(flowTraj);
    
    antFlowTraj = diff(rotTD(~mask).brainKin.pos([end end-19],:))';
    antFlowTraj = antFlowTraj./norm(antFlowTraj);
else
    flowTraj = [];
    antFlowTraj = [];
    rotTD = IT.TDunconstrained.average;
end
% Calculate the angle difference for the unconstrained trials
TD = IT.TDunconstrained;
if successOnly
    TD = TD([TD.successful]==1);
end
TD = TD(ismember([TD.startPos]',startPos,'rows'));
AD.avgUncon = util.calc_angle_diff(TD.average('avgMode','samp'),...
    rotTD,'flowTraj',flowTraj,'antFlowTraj',antFlowTraj,extraOpts{:});
AD.uncon = util.calc_angle_diff(TD, rotTD,'flowTraj',flowTraj,...
    'antFlowTraj',antFlowTraj,extraOpts{:});
AD.uncon(end+1,:) = [TD.successful];

% Calculate the angle difference for the constrained trials
for n = 1:size(IT.TDconstrained,1)
    TD = IT.TDconstrained{n};
    if successOnly
        TD = TD([TD.successful]==1);
    end
    TD = TD(ismember([TD.startPos]',startPos,'rows'));
    aTD = TD.average('avgMode','samp');
    if n == 1 && plotFig
        F = figure; hold on
        for k = 1:size(rotTD,1)
            rotTD(k).targPos = rotTD(k).startPos;
            rotTD(k).startPos = -rotTD(k).startPos;
        end
        rotTDplt = IT.TD_tt_rot{1}.average('avgMode','samp');
        rotTDplt.plot('trajColorWeight',.5,'ColMat',C,'plotTargNum',0);
        aTD.plot('color',[.7 .7 .7],'plotTargNum',0)
        [AD.avgCon(:,n),~,AD.compAng] = util.calc_angle_diff(aTD,rotTD,...
            'flowTraj',flowTraj,'antFlowTraj',antFlowTraj,'F',F,...
            'plotAng',1,extraOpts{:},'plotScale',50,'offset',startPos);
        title(sprintf('%s%s Average trajectories and initial angles',...
            IT.subject,IT.date))
        axis(175*[-1 1 -1 1])
        F.Name = sprintf('%s%s_initialAngle',IT.subject,IT.date);
    elseif n == 1
        [AD.avgCon(:,n),~,AD.compAng] = util.calc_angle_diff(aTD,rotTD,...
            'flowTraj',flowTraj,'antFlowTraj',antFlowTraj,extraOpts{:});
    else
        AD.avgCon(:,n) = util.calc_angle_diff(aTD,rotTD,...
            'flowTraj',flowTraj,'antFlowTraj',antFlowTraj,extraOpts{:});
    end
    AD.con{n} = util.calc_angle_diff(TD, rotTD,'flowTraj',flowTraj,...
        'antFlowTraj',antFlowTraj,extraOpts{:});
    AD.con{n}(end+1,:) = [TD.successful];
end

% Calculate the null distribution, i.e. the initial angle for the intuitive
% mapping two target task.
if useRot
    % Calculate intTraj
    intTrl = IT.TDconstrained{1}(1);
    intTraj = intTrl.targPos(1:2) - intTrl.startPos(1:2);
    intTraj = intTraj./sqrt(sum(intTraj.^2));

    TD = IT.TD_tt_rot{1};
    TD = TD(ismembertol([TD.startPos]',startPos,1e-6,'ByRows',true));
    AD.rotTT = util.calc_angle_diff(TD, rotTD,'flowTraj',flowTraj,...
        'antFlowTraj',antFlowTraj,'intTraj',intTraj,'colMat',jet(size(TD,1)),...
        'plotAng',plotFig,extraOpts{:});
end
end