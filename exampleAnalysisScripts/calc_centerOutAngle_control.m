% This example script will run through how to calculate the initial angle
% values for the center out trials of a single experiment. This is the data
% used in figure 6 and 7 of the paper. The initial angles are calculated
% for the average and individual trials of each unique target angle. For
% ease of comparison with the constrained task data, the with flow and
% against flow angles are defined as the targets directly to the left and
% right of the target. The final result will be a structure called 
% <AD_co> which will contain the initial measure for each unique target 
% direction. This save structure will also store the average trajectory for
% each target for a visual confirmation. The user has the option to save 
% the data as well.
%
% The batch version of this code is <tube.batch_initial_angle_centerOutControl>
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com
clear,close all

% Add the paper code to the path
pathName = pwd;
addpath(genpath(pathName))

% Set the example session and if you want to save the data
exampleSess = {'20190719'};
centerPos = [0 0 0];
plotAng = 1;                % Determine whether to plot the angle comparison
                            % for each target
saveData = 0;               % Determine whether or not to save the data
savePath = [];              % Saving location for the data.

if saveData
    savePath = uigetdir('',"Where do you want to save the d' data?");
end

% Load in the D structure
dataLoc = serverPath;
D = load(fullfile(dataLoc,'exampleDatasetCatalog.mat'));
D = D.D;

% Determine the session with the correct data
D = D(ismember({D.dataset},exampleSess));
dir_list = db.get_task_datasets(D, {'inttarg_rot_constr_slow'});

% Create the AD file based on the number of sessions
AD_co = repmat(struct('subject',[],'date',[],'targAng',[],...
    'avgAng_ZIA',[],'avgTD',[],'indAng_ZIA',[],'baseAng',[],...
    'intTrajVec',[],'OpposeTrajVec',[],'zeroTrajVec',[]),8,size(dir_list,1));

%% Iterate through the session
for k = 1:size(dir_list,1)
    % Find and load the gradual training data.
    pathName = fullfile(dir_list(k).base,'translated','trajectoryData');
    dirInfo = dir(pathName);
    fileIdx = find(contains({dirInfo.name},'centerOut'));
    if isempty(fileIdx)
        warning('Could not find the correct data set for session: %s',...
            dir_list(k).dataset)
        continue
    end
    fileLoc = fullfile(dirInfo(fileIdx).folder,dirInfo(fileIdx).name);

    % Load the data
    TD = TrajectoryData().load(fileLoc);

    % Normalize the data
    if ~ismember(mean(unique([TD.targPos]')),centerPos)
        centerPos = mean(unique([TD.targPos]','rows'));
    end
    TD  = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);    

    % Identify the target directions
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

    % Set the offset to 45 degrees and find the 0 degree target which we
    % will use as the basline comparison
    posZero = find(ismembertol(ang_aTD,0,1e-5));
    offsetAng = 45;
    
    % Iterate through each unique target location and identify the angle
    % difference
    for n = 1:size(uniAng,2)
        % Create the intTraj variable for the main target angle
        posAng = find(ismembertol(ang_aTD,uniAng(n),1e-5));
        tempOff = [aTD(posAng).startPos' aTD(posAng).targPos'];
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
        
        % Calculate the initial angles for the average
        [avgAng_ZIA,~,baseAng] = util.calc_angle_diff(aTD(posAng),...
            aTD(posZero),'intTraj',intTraj,'flowTraj',flowTraj,...
            'antFlowTraj',antFlowTraj,'compInt',1,'plotAng',plotAng);
        if plotAng
            title(num2str(uniAng(n)))
        end
        % Calculate the initial angles for individual trials
        tempTD = TD(idxAng==n);
        [indAng_ZIA] = util.calc_angle_diff(tempTD,aTD(posZero),...
            'intTraj',intTraj,'flowTraj',flowTraj,...
            'antFlowTraj',antFlowTraj,'compInt',1);
        indAng_ZIA = [indAng_ZIA; [tempTD.successful]];

        % Save the data
        AD_co(n,k).indAng_ZIA(end+1,:) = [tempTD.successful];
        AD_co(n,k).avgTD = aTD(posAng);
        AD_co(n,k).TD = tempTD;
        AD_co(n,k).avgAng_ZIA = avgAng_ZIA;
        AD_co(n,k).indAng_ZIA = indAng_ZIA;
        AD_co(n,k).baseAng = baseAng;
        AD_co(n,k).targAng = uniAng(n);
        AD_co(n,k).intTrajVec = intTraj;
        AD_co(n,k).zeroTrajVec = flowTraj;
        AD_co(n,k).OpposeTrajVec = antFlowTraj;
        AD_co(n,k).subject = dir_list(k).subject;
        AD_co(n,k).date = dir_list(k).dataset;
    end
end

% Save the control angle data structure to a single mat file.
if ~isempty(saveDataPath)
    fileName = sprintf('centerOut_angle_example%s.mat',exampleSess{:});
    save(fullfile(saveDataPath,fileName),'AD_co')
end