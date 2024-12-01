% tube.batch_initial_angle  Run the initial angle analysis on all datasets
%
% Usage:
%   tube.batch_initial_angle
%
% This function iterates over all two-target intuitive, rotated, 
% unconstrained, and tube constrained datasets and 
% runs tube.calc_inital_angle_diff function. 
%
% Optional Inputs:
%   data_loc    Location of the data folder.
%   D           Structure with all the valid experimental data. Save as
%                   filename <<exampleDatasetCatalog.mat>>
%   plotFig     Plot initial angle plots of average trajectories
%   save_file   Determine whether or not to save the AD_compare struct
%   save_path   Determine where to save the AD_compare data
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [AD_compare] = batch_initial_angle(varargin)
% Calculated the initial angles for the inituitive and rotated projection
% of the constrained tasks.
dataLoc = [];
D = [];
plotFig = 0;
save_file = 0;
save_path = [];
task = {'inttarg_rot_constr_slow'};

% Parse optional agruments
assignopts(who,varargin);

% Define the data location
if isempty(dataLoc)
    dataLoc = serverPath();
end

% Load the relevant data
if isempty(D)
    D = load(fullfile(dataLoc,'exampleDatasetCatalog.mat'));
    D = D.D;
end

% Determine the sessions with the correct data
dir_listIT = db.get_task_datasets(D,task);
int_targ_loc = fullfile(dataLoc,'ConstrainedPath','mat','int_targ_data');

% Determine the valid datasets
sessInfo = cell(size(dir_listIT));
for n = 1:size(dir_listIT,1)
    sessInfo{n} = [dir_listIT(n).subject dir_listIT(n).dataset '_int_targ.mat'];
end
valid_files = util.findDirContents(int_targ_loc, '_int_targ.mat');
valid_idx = ismember(sessInfo,valid_files);

% Iterate through the sessions
for n = 1:length(sessInfo)
    if valid_idx(n)
        IT = load(fullfile(int_targ_loc, sessInfo{n}));
        IT = IT.IT;
    else
        continue
    end
    tD = D(ismember({D.dataset},dir_listIT(n).dataset));

    % Identify the two target intuitive and rotated trials and then load the
    % the trajectory data and the decoders.
    dir_list = db.get_task_datasets(tD, {'tt_int','tt_rot'});

    TD_int = util.loadSessionData(dir_list(1));
    [~, P_rot, result_rot] = util.loadSessionData(dir_list(2));

    % Normalize the data
    centerPos = [0 0 0];
    if ~ismember(mean(unique([TD_int.startPos]')),centerPos)
        TD_int = TD_int.normalize(mean(unique([TD_rot.startPos]','rows')));
    else
        TD_int = TD_int.normalize(centerPos);
    end

    % Only include the two target data
    TD_int = TD_int(ismember([TD_int.targPos]',...
        unique([TD_int.startPos]','rows'),'rows'));

    % Calculate the normal initial angles for the rotated spaces. (rotated
    % project, unconstrained. and intuitive tt

    % Base comparison is to the rotated two target task
    [AD_compare(n,1)] = tube.calc_inital_angle_diff(IT,'successOnly',0,...
        'decoderVisual','rot','useRot',1,'compInt',1,'plotFig',plotFig);
    AD_compare(n,1).subject = [AD_compare(n,1).subject ' tt_SM proj_SM'];

    % Base comparison to the intuitive two target task projected into the
    % rotated mapping
    tIT = changeTTrotProj(IT,TD_int,P_rot,result_rot,centerPos);
    [AD_compare(n,2)] = tube.calc_inital_angle_diff(tIT,'successOnly',0,...
        'decoderVisual','rot','useRot',1,'compInt',1,'plotFig',plotFig);
    AD_compare(n,2).subject = [AD_compare(n,2).subject ' tt_IM proj_SM'];
end

% Remove the invalid session
AD_compare = AD_compare(valid_idx==1,:);

% Collect the data into simplified versions and then save it.

% Rotated mapping data
rotDat = [AD_compare(:,1)];
uncon_angs_rot = [rotDat.avgUncon];
ang_unCon2intTarg_r = uncon_angs_rot(2,:);
tt_angs_rot = [rotDat.compAng];
ang_flow2intTarg_r = tt_angs_rot(1,:);
ang_antFlow2intTarg_r = tt_angs_rot(3,:);

% Inituitive mapping data
intDat = [AD_compare(:,2)];
uncon_angs_int = [intDat.avgUncon];
ang_unCon2intTarg_i = uncon_angs_int(2,:);
tt_angs_int = [intDat.compAng];
ang_tt2intTarg_i = tt_angs_int(1,:);
ang_flow2intTarg_i = tt_angs_int(1,:);
ang_antFlow2intTarg_i = tt_angs_int(3,:);

% Save the data
if save_file
    if isempty(save_path)
        save('initialAng_tt_vs_uncon+con_Simplified.mat','AD_compare',...
            'ang_antFlow2intTarg_i','ang_antFlow2intTarg_r',...
            'ang_flow2intTarg_i','ang_flow2intTarg_r','ang_tt2intTarg_i',...
            'ang_unCon2intTarg_i','ang_unCon2intTarg_r')
    else
        save(fullfile(save_path,'initialAng_tt_vs_uncon+con_Simplified.mat'),...
            'AD_compare','ang_antFlow2intTarg_i','ang_antFlow2intTarg_r',...
            'ang_flow2intTarg_i','ang_flow2intTarg_r','ang_tt2intTarg_i',...
            'ang_unCon2intTarg_i','ang_unCon2intTarg_r')
    end
end

end

% Support Function
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