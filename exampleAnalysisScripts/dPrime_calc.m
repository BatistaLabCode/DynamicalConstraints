% This example script will run through how to calculate the d' for IM trial
% data projected into both IM and SM projections. This is the data that is
% used in figure 3. The final result will be a structure called <datROC_WS> 
% that will contain the d', the true positive rate (TPR), the true negative
% rate (TNR), the auc, and t-test results for each condition. There is also
% another variable <dP_WS> that will just store d' values. Lastly, the code
% will plot the distribution and ROC curves for each condition for visual 
% in spection. The user has the option to save the data as well.
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

% Create a structure that will store the session and animal information
sessInfo = repmat(struct('subject',[],'dataset',[],'condition',[]),...
    size(dir_list,1),2);

% Iterate through intuitive mapping (IM) and separation max (SM) workspaces conditions
invalid = []; % Matrix to save any sessions that failed calculations
for n = 1:size(dir_list,1)

    % Load the decoder
    [TD_int, P_int, result_int] = util.loadSessionData(dir_list(n,1));
    [~, P_rot, result_rot] = util.loadSessionData(dir_list(n,2));

    % Normalize the data
    centerPos = [0 0 0];
    if ~ismember(mean(unique([TD_int.startPos]')),centerPos)
        TD_int = TD_int.normalize(mean(unique([TD_int.startPos]','rows')));
    else
        TD_int = TD_int.normalize(centerPos);
    end

    % Exclude any trials that were not part of the two target pair.
    TD_int = TD_int(ismember([TD_int.targPos]',unique([TD_int.startPos]','rows'),'rows'));

    % Determine the new animal name
    if contains(dir_list(n,1).subject,'monkeyD')
        subject = 'Monkey D';
    elseif contains(dir_list(n,1).subject,'monkeyE')
        subject = 'Monkey E';
    elseif contains(dir_list(n,1).subject,'monkeyQ')
        subject = 'Monkey Q';
    end

    for k = 1:2
        % Determine which gpfa and TD condition we are using
        if k == 1 % IM trials, IM projection
            P = P_int;
            result = result_int;
            condStr = 'IM trials, IM projection';
        elseif k == 2 % IM trials, SM projection
            P = P_rot;
            result = result_rot;
            condStr = 'IM trials, SM projection';
        end

        % Add the GPFA data by processing the binned spiking data through
        % the decoder used in the behavior task. Include the time step when
        % the target was acquired.
        [~,~,TT] = orthogonalize(zeros(P.xDim,1),result.estParams.C);
        TD = util.predictDecodeState_GPFA(TD_int,P,'TT',TT,'predictPos',1,...
            'spikeCountSrc','decodeSpikeCounts','useT',1,'trunGPFA',1,'timeStep',1);

        % Calculate d' for the session
        try
            [datROC_WS(n,k),dP_WS(n,k),optMet_WS(n,k)] = SDT.calcWorkspaceDPrime(TD,'pos');
            if ismember({dir_list(n,1).dataset},exampleSess)
                SDT.plotROC(datROC_WS(n,k))
                plt.plotTitle(['ROC Workspace: ' dir_list(n,1).dataset ' ' condStr])
            end
        catch % Save the failed trials
            invalid = [invalid; n k];
        end

        % Add the session and animal information to the function.
        sessInfo(n,k).subject = subject;
        sessInfo(n,k).dataset = dir_list(n,1).dataset;
        sessInfo(n,k).condition = condStr;
    end
end

% Save the d' data structure to a single mat file.
if ~isempty(savePath)
    fileName = sprintf('dPrime_example%s.mat',exampleSess{:});
    save(fullfile(savePath,fileName),'datROC_WS','dP_WS',...
        'optMet_WS','sessInfo')
end