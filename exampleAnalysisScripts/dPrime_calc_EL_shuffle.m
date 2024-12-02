% This example script will run through how to calculate the d' early, d'
% late, and the shuffled d' for the IM trials and the SM trials. This is 
% the data that is used in figure 4. The final result will be a structure 
% called <dat_dPrime_EL> which will contain all the d' data for the two
% conditions. The user has the option to save the data as well.
%
% Note: The shuffle data in the file <dPrime_EarlyLate+ShuffledData> and 
% used in the main figures was collected without a seeded random value,
% therefore the shuffle data will be different each iteration, but trends
% will be consistent.
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

clear, close all % Clear and set the workspace

% Add the paper code to the path
pathName = pwd;
addpath(genpath(pathName))

% Set the example session and if you want to save the data
exampleSess = {'20190719'};
saveData = 0;               % Determine whether or not to save the data
savePath = [];              % Saving location for the data.

if saveData
    savePath = uigetdir('',"Where do you want to save the d' data?");
end

% Load in the D structure
dataLoc = serverPath;
D = load(fullfile(dataLoc,'exampleDatasetCatalog.mat'));
D = D.D;

% Determine the sessions with the correct data
dir_list = db.get_task_datasets(D, {'tt_int','tt_rot'});

% Recalculate the d'
invalid = [];
Alpha = 0.05;
bootRep = 1000;
% Create the data structure for the file
dat_dPrime_EL = repmat(struct('subject',[],'dataset',[],'info',[],...
    'dP_all',[],'dP_E',[],'dP_L',[],'dP_shuff1',[],'dP_shuff2',[]),...
    size(dir_list,1),2);

% iterate through intuitive and rotated workspaces conditions
for n = 1:size(dir_list,1)

    % Load the decoder
    [TD_int, P_int, result_int] = util.loadSessionData(dir_list(n,1));
    [TD_rot, P_rot, result_rot] = util.loadSessionData(dir_list(n,2));

    % Normalize the data
    centerPos = [0 0 0];
    if ~ismember(mean(unique([TD_rot.startPos]')),centerPos)
        TD_int = TD_int.normalize(mean(unique([TD_rot.startPos]','rows')));
        TD_rot = TD_rot.normalize(mean(unique([TD_rot.startPos]','rows')));
    else
        TD_int = TD_int.normalize(centerPos);
        TD_rot = TD_rot.normalize(centerPos);
    end
    TD_int = util.preprocessGridTaskTrajData(TD_int,'centerPos',centerPos);
    TD_rot = util.preprocessGridTaskTrajData(TD_rot,'centerPos',centerPos);
    TD_int = TD_int(ismember([TD_int.targPos]',unique([TD_int.startPos]','rows'),'rows'));
    TD_rot = TD_rot(ismember([TD_rot.targPos]',unique([TD_rot.startPos]','rows'),'rows'));

    for k = 1:2
        % Determine which gpfa and TD condition we are using
        if k == 1 % Intuitive trials, intuitive projection
            TD = TD_int;
            P = P_int;
            result = result_int;
            info = 'TD_int + Decode_Int';
        elseif k == 2 % Rotated trials, rotated projection
            TD = TD_rot;
            P = P_rot;
            result = result_rot;
            info = 'TD_SM + Decode_SM';
        end

        % Add the GPFA data, including the time step where the start
        % target was acquired.
        [~,~,TT] = orthogonalize(zeros(P.xDim,1),result.estParams.C);
        TD = util.predictDecodeState_GPFA(TD,P,'TT',TT,'predictPos',1,...
            'spikeCountSrc','decodeSpikeCounts','useT',1,'trunGPFA',1,'timeStep',1);

        % Calculate d' prime per session, full session, early, late,
        % shuffled early, and shuffled late.
        [~,dat_dPrime_EL(n,k).dP_all] = SDT.calcWorkspaceDPrime(TD,'pos');
        szTD = length(TD);
        midPt = round(szTD/2);

        % Early vs Late
        [~,dat_dPrime_EL(n,k).dP_E] = SDT.calcWorkspaceDPrime(TD(1:midPt),'pos');
        [~,dat_dPrime_EL(n,k).dP_L] = SDT.calcWorkspaceDPrime(TD(midPt+1:end),'pos');

        % Shuffled twice
        idx = randperm(szTD);
        [~,dat_dPrime_EL(n,k).dP_shuff1] = SDT.calcWorkspaceDPrime(TD(idx(1:midPt)),'pos');
        [~,dat_dPrime_EL(n,k).dP_shuff2] = SDT.calcWorkspaceDPrime(TD(idx(midPt+1:end)),'pos');

        % Add the experimental information
        dat_dPrime_EL(n,k).subject = TD(1).subject;
        dat_dPrime_EL(n,k).dataset = dir_list(n,1).dataset;
        dat_dPrime_EL(n,k).info = info;
    end
end

% Save the d' data structure to a single mat file.
if ~isempty(savePath)
    fileName = sprintf('dPrimeShuffle_example%s.mat',exampleSess{:});
    save(fullfile(savePath,fileName),'dat_dPrime_EL')
end
