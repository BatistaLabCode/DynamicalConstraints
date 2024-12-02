% This example script will run through how to calculate the SM projection
% for a given dataset and then apply it to the neural activity for a visual
% respresentation. This is the same process used to identify the SM
% projection for all experiment. The code uses data from the example
% session 20190719, but process is valid for other neural data (see steps
% starting at line 58).
% 
% Note: For Monkey E, we assumed equal scaling weights for the objective 
% function (i.e. calcWeights = 0), while the other animals we identify the 
% optimal scaling weights. That said, the identified projections for monkey
% E were nearly identical between the two methods.
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

clear, close all

% Add the paper code to the path
pathName = pwd;
addpath(genpath(pathName))

% Set the example session and if you want to save the data
exampleSess = {'20190719'}; % The example session that we plot the neural
                            %    space trajectories. (From Monkey E)
calcWeights = 1;            % Determine the scaling weights for each 
                            %   objective in the optimization function. If 
                            %   set to 0 than the scale will all be 1.
compare2Online = 1;         % Plot the projection found in the experiment.
tasks = {'tt_int','tt_rot'}; % Tasks 

% Determine the sessions with the correct data
[dataLoc]= serverPath
load(fullfile(dataLoc,'exampleDatasetCatalog.mat'));
dir_list = db.get_task_datasets(D,tasks);

% Identify the example session
maskSessEx = find(ismember({dir_list(:,1).dataset},exampleSess));

% Load the data and the SM decoder
[TD, P, resultInt] = util.loadSessionData(dir_list(maskSessEx,1));
TD = TD.normalize([0 0 0]);
TD = TD(ismember([TD.targPos]',unique([TD.startPos]','rows'),'rows'));
[uniTarg,~,condition] = unique([TD.targPos]','rows');
if size(tasks,2)>=2
    [~, Prot,resultRot] = util.loadSessionData(dir_list(maskSessEx,2));
    [~,~,TTrot] = orthogonalize(zeros(Prot.xDim,1),resultRot.estParams.C);
else % If no SM decoder exists, use the neural space of the IM decoder
    Prot = P;
    [~,~,TTrot] = orthogonalize(zeros(P.xDim,1),resultInt.estParams.C);
end

% Project the neural data into the 10D neural space where the separation
% maximizing projection was identified. We will truncate the data from the
% sample right before the trajectory onset index to the thajectory offset
% index. The rationale for this is that this is the timestep where the
% target was aquired in the previous state, and thus is closest to the
% target position.
arotTD = util.predictDecodeState_GPFA(TD,Prot,'trunGPFA',1,'timeStep',1,...
    'spikeCountSrc','decodeSpikeCounts','useT',1,'TT',TTrot,'predictPos',1);
gp2D = [arotTD.GPFA];

%% Determine the SM projection using the gpfa data. Note that gp2D is 
% easily replaceable. All that is necessary is a clear mapping betwee the 
% neural data and the target pair condition identified as either 1 or 2. 

% Iterate through the number of condition pairs
for n = 1:size(uniTarg,1)/2
    idx = find(ismember(uniTarg,-1*uniTarg(n,:),'rows'));
    mask = ismember(condition,[n idx]);
    tmpCond = condition(mask==1);
    tmpCond(tmpCond == n) = 1;
    tmpCond(tmpCond == idx) = 2;
    optMet = opt.get_data({gp2D(mask==1).xorth},tmpCond);

    % Calculate the scaling weights
    if calcWeights % This will calculate the optimal scaling weight for each objective
        [objFnParams] = opt.calcObjFnWeights(optMet);
    else % This will use equal weights for each objective.
        objFnParams.w_mid = 1;
        objFnParams.w_var = 1;
        objFnParams.w_start = 1;
    end

    % Run the optimization function
    M{n} = opt.stiefelOpt(optMet,@opt.standard_obj,@opt.standard_grad,...
        'objFnParams', objFnParams,'verbose',0);
end

%% Plot the trajectories in the identified projection.

% Define the color maps and set up the figure
    C = util.defineTaskColormap('bc_int');

% Create the optional inputs
optArg = {'plotAxLabels',false,...
    'plotPreOnset', false, ...
    'xSpec', 'xorth', ...
    'plotDim', [1 2], ...
    'colorMode','specified',...
    'lineStr', '-',...
    'lineWidth',.25,...
    'markerSize',3};

% Iterate through the number of condition pairs
for n = 1:size(uniTarg,1)/2
    idx = find(ismember(uniTarg,-1*uniTarg(n,:),'rows'));

    % Apply the identify projection to the neural data
    projDat = opt.applyProjection(arotTD,M{n});

    % Determine the indices for different color schemes
    aIdx = find(C.endTC == find(ismember(round(C.targPos),round(uniTarg(n,:)),'rows')));
    bIdx = find(C.endTC == find(ismember(round(C.targPos),round(uniTarg(idx,:)),'rows')));

    % Create the gpfa variables
    gpa = [projDat.GPFA];
    gpfaA = gpa(condition==n);
    gpfaB = gpa(condition==idx);

    % Plot the projection
    F(n) = figure;
    if compare2Online % Plot side-by-side comparison of the identified vs online projections
        subplot(1,2,1)
        gpfaA.trajectory('col',C.col_light(aIdx,:),optArg{:})
        gpfaB.trajectory('col',C.col_light(bIdx,:),optArg{:})
        view(90,90) % Rotates the space to align with the view in the paper.

        % Set the axis and label information
        axis square,grid on
        xlabel('Neural dimension 2'),ylabel('Neural dimension 1')
        title(sprintf('Identified Projection: Target Pair [%0.0f %0.0f]',n,idx))

        subplot(1,2,2)
        onlineProjDat = opt.applyProjection(arotTD,Prot.ExpInfo.OptProj.M);
        gpaOnline = [onlineProjDat.GPFA];
        gpaOnline(condition==1).trajectory('col',C.col_light(aIdx,:),optArg{:})
        gpaOnline(condition==2).trajectory('col',C.col_light(bIdx,:),optArg{:})
        title(sprintf('Online SM Projection: Target Pair [%0.0f %0.0f]',n,idx))
    else
        gpfaA.trajectory('col',C.col_light(aIdx,:),optArg{:})
        gpfaB.trajectory('col',C.col_light(bIdx,:),optArg{:})
        title(sprintf('Identified Projection: Target Pair [%0.0f %0.0f]',n,idx))
    end
    
    % Set the axis view and label information
    view(90,90) % Rotates the space to align with the view in the paper.
    axis square,grid on
    xlabel('Neural dimension 2'),ylabel('Neural dimension 1')
end