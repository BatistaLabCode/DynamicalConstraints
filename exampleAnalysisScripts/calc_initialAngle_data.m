% This example script will run through how to calculate the initial angle
% values for the unconstrained and constrained trials of an experiment. 
% This is the data used in figure 6 and 7 of the paper. The initial angles 
% are calculated for the average and the individual trials for the 
% unconstrained, each constrained tube boundary, and the two target 
% conditions. These values and the constrained tube boundary sizes are
% stored in the structure <AD_compare>. Simplify angle vectors was also
% created that are used in fig6_fig7 code. The user has the option to save 
% the data as well.
%
% The batch version of this code is <tube.batch_initial_angle>
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

% Find and load the valid IT structures
data_save_loc = fullfile(dataLoc,'ConstrainedPath','mat','int_targ_data');
valid_files = util.findDirContents(data_save_loc, '_int_targ.mat');

for n = 1:length(valid_files)
    IT = load(fullfile(data_save_loc, valid_files{n}));
    IT = IT.IT;
    tD = D(ismember({D.dataset},IT.date)); % Hardcoded since the naming
                                           %   structure is consistent

    % Identify the two target intuitive and rotated trials and then load the
    % the trajectory data and the decoders.
    dir_list = db.get_task_datasets(tD, {'tt_int','tt_rot'});

    TD_int = util.loadSessionData(dir_list(1));
    [~, P_rot, result_rot] = util.loadSessionData(dir_list(2));

    % Normalize the data
    if ~ismember(mean(unique([TD_int.startPos]')),centerPos)
        TD_int = TD_int.normalize(mean(unique([TD_int.startPos]','rows')));
    else
        TD_int = TD_int.normalize(centerPos);
    end

    % Calculate the normal initial angles for the rotated spaces. 
    % (SepMax trials, unconstrained trials, and IM trials)

    % Base comparison is to the SepMax two target task
    [AD_compare(n,1)] = tube.calc_inital_angle_diff(IT,'successOnly',0,...
        'decoderVisual','rot','useRot',1,'compInt',1,'plotFig',0);
    AD_compare(n,1).subject = [AD_compare(n,1).subject ' tt_SM proj_SM'];

    % Base comparison to the IM two target task projected into the
    % SepMax mapping
    tIT = changeTTrotProj(IT,TD_int,P_rot,result_rot,centerPos);
    [AD_compare(n,2)] = tube.calc_inital_angle_diff(tIT,'successOnly',0,...
        'decoderVisual','rot','useRot',1,'compInt',1,'plotFig',0);
    AD_compare(n,2).subject = [AD_compare(n,2).subject ' tt_IM proj_SM'];
end

% Create simplified lists of the angle measurements for comparison.
% SM mapping data
rotDat = [AD_compare(:,1)];
uncon_ang_SM = [rotDat.avgUncon];
ang_unCon2intTarg_r = uncon_ang_SM(2,:); % Angle of the average uncon.
                                         %  trajectory to int target.
tt_ang_SM = [rotDat.compAng];
ang_flow2intTarg_r = tt_ang_SM(1,:); % Angle of the average flow two target
                                       %  trajectory to int target.
ang_antFlow2intTarg_r = tt_ang_SM(3,:); % Angle of the average against flow
                                        %  two target trajectory to int target.

% IM mapping data
intDat = [AD_compare(:,2)];
tt_ang_IM = [intDat.compAng];
ang_flow2intTarg_i = tt_ang_IM(1,:); % Angle of the average flow two target
                                     %  trajectory to int target.
ang_antFlow2intTarg_i = tt_ang_IM(3,:); % Angle of the average against flow
                                        %  two target trajectory to int target.

% Save the inital angle data structures to a single mat file.
if ~isempty(saveDataPath)
    fileName = sprintf('initialAngle_example%s.mat',exampleSess{:});
    save(fullfile(saveDataPath,fileName),'AD_compare',...
        'ang_flow2intTarg_i','ang_flow2intTarg_r','ang_unCon2intTarg_r')
end
%% Functions
function [IT] = changeTTrotProj(IT,TD,P,result,centerPos)
% Remove non-paired targets from trials
TD = TD(ismember([TD.targPos]',unique([TD.startPos]','rows'),'rows'));

% Find orthonormal
[~,~,TT] = orthogonalize(zeros(P.xDim,1),result.estParams.C);

% Add the GPFA data for intuitive trials
TD = util.predictDecodeState_GPFA(TD,P,'spikeCountSrc','decodeSpikeCounts',...
    'useT',1,'TT',TT,'predictPos',1);

% Normalized and set the target center
TDnorm = util.preprocessGridTaskTrajData(TD,'centerPos',centerPos);

IT.TD_tt_rot{1} = TDnorm;
end