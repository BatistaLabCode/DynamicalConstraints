% tube.batch_initial_angle_centerOutControl  Run the initial angle control
% analysis on all datasets
%
% Usage:
%   tube.batch_initial_angle_centerOutControl
%
% This function iterates over every final block of center out datasets
% and runs util.calc_angle_diff function for each target. This dataset can
% then be used as a control of full change trajecotries.
%
% Optional Inputs:
%   data_loc    Location of the data folder
%   centerPos   Center of the workspace
%   save_file   Determine whether or not to save the AD_compare struct
%   save_path   Determine where to save the AD_compare data
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [AD_co] = batch_initial_angle_centerOutControl(varargin)
%% Calculated the initial angles for the gradial training/center out task.
data_loc = [];
centerPos = [0 0 0];
save_file = 0;
save_path = [];

% Parse optional agruments
assignopts(who,varargin);

% Define the data location
if isempty(data_loc)
    data_loc = serverPath();
end

% Determine the valid datasets
int_targ_loc = fullfile(data_loc,'ConstrainedPath\mat\int_targ_data');
valid_files = util.findDirContents(int_targ_loc, '_int_targ.mat');

% Create the AD file based on the number of sessions
AD_co = repmat(struct('subject',[],'date',[],'targAng',[],...
    'avgAng_ZIA',[],'avgTD',[],'indAng_ZIA',[],'baseAng',[],...
    'intTrajVec',[],'OpposeTrajVec',[],'zeroTrajVec',[]),8,length(valid_files));

%% Iterate through the sessions
for k = 1:size(valid_files,2)

    % Load the IT data
    IT = load(fullfile(int_targ_loc, valid_files{k}));
    IT = IT.IT;
    
    % Collect the session date information
    yr = IT.date(1:4);
    mon = IT.date(5:6);

    % Find and load the center out data.
    pathName = fullfile(dataLoc,IT.subject,yr,mon,IT.date,'translated','trajectoryData');
    dirInfo = dir(pathName);
    fileIdx = find(contains({dirInfo.name},'centerOut'));

    % Load the data
    TD = TrajectoryData().load(fullfile(pathName,dirInfo(fileIdx(end)).name));

    % Normalize the data
    if ~ismember(mean(unique([TD.targPos]')),centerPos)
        centerPos = mean(unique([TD.targPos]','rows'));
    end
    TD  = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

    %% Identify the target directions
    % Calculate the average TD
    aTD = TD.average('avgMode','samp');
    targPos_avg = [aTD.targPos];
    ang_aTD = atan2d(targPos_avg(2,:),targPos_avg(1,:));
    ang_aTD(ang_aTD<0) = ang_aTD(ang_aTD<0) + 360;

    % Determine the angle for each unique trial
    targPos = [TD.targPos];
    angTD = atan2d(targPos(2,:),targPos(1,:));
    angTD(angTD<0) = angTD(angTD<0) + 360;
    [uniAng,~,idxAng] = unique(angTD);
    uniAng =round(uniAng);

    % Set the flowTraj variable as the 0 degree target
    posZero = find(ismembertol(ang_aTD,0,1e-5));
    offsetAng = 45;

    % Iterate through each unique target location and identify the angle
    % difference
    for n = 1:size(uniAng,2)
        posAng = find(ismembertol(ang_aTD,uniAng(n),1e-5));
        tempOff = [aTD(posAng).startPos' aTD(posAng).targPos'];

        % Create the intTraj variable
        intTraj = (tempOff(4:5) - tempOff(1:2))';
        intTraj = intTraj./sqrt(sum(intTraj.^2));

        % Set the flowTraj variable as the target + offset degrees from uniAng(n)
        posPoffset = find(ismembertol(ang_aTD,mod(uniAng(n)+offsetAng,360),1e-5));
        flowTraj = aTD(posPoffset).targPos(1:2) - aTD(posPoffset).startPos(1:2);
        flowTraj = flowTraj./sqrt(sum(flowTraj.^2));

        % Create the antflowTraj variable as the target - offset degrees from uniAng(n)
        posNoffset = find(ismembertol(ang_aTD,mod(uniAng(n)-offsetAng,360),1e-5));
        antFlowTraj = aTD(posNoffset).targPos(1:2) - aTD(posNoffset).startPos(1:2);
        antFlowTraj = antFlowTraj./sqrt(sum(antFlowTraj.^2));

        tempTD = TD(idxAng==n);

        % Calculate the initial angles
        [AD_co(n,k).avgAng_ZIA,~, AD_co(n,k).baseAng] = util.calc_angle_diff(aTD(posAng),...
            aTD(posZero),'intTraj',intTraj,'flowTraj',flowTraj,...
            'antFlowTraj',antFlowTraj,'compInt',1);

        % Save the data
        [AD_co(n,k).indAng_ZIA] = util.calc_angle_diff(tempTD,aTD(posZero),...
            'intTraj',intTraj,'flowTraj',flowTraj,...
            'antFlowTraj',antFlowTraj,'compInt',1);
        AD_co(n,k).indAng_ZIA(end+1,:) = [tempTD.successful];
        AD_co(n,k).avgTD = aTD(posAng);
        AD_co(n,k).TD = tempTD;
        AD_co(n,k).targAng = uniAng(n);
        AD_co(n,k).intTrajVec = intTraj;
        AD_co(n,k).zeroTrajVec = flowTraj;
        AD_co(n,k).OpposeTrajVec = antFlowTraj;
        AD_co(n,k).subject = IT.subject;
        AD_co(n,k).date = IT.date;
    end
end

% Save the data
if save_file
    if isempty(save_path)
        save('centerOut_AngleAnalysisData_Signed.mat','AD_co')
    else
        save(fullfile(save_path,'centerOut_AngleAnalysisData_Signed.mat'),'AD_co')
    end
end

