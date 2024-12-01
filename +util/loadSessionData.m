% [data,P,result] = loadSessionData(dir_list,varargin)  Loads the condition 
% relevant data including decoder and latent space parameters
%
% Usage:
%   util.loadSessionData(dir_list)
%
% Inputs:
%   dir_list    Structure with the location information of all files.
%
% Optional Inputs:
%   dataType    Type of data structure to load
%   centerPos   Center of the workspace
%   preprocess  Tag to normalize TrajectoryData
%   idx         Location in a subcell of dir_list
%   gpfaName    Basename of the latent space parameter file
%   loadResult  Tag to load the latent space parameter data
% 
% Outputs:
%   data        Data structure
%   P           BCI decoder associated with data structure
%   result      Associated latent space parameter file.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [data,P,result] = loadSessionData(dir_list,varargin)
%% Loads the data, the decoder, and the results
dataType = 'traj'; %traj, trial, data, export
centerPos = [0,0,0];
preprocess = true;
idx = 1;
gpfaName = 'gpfa_xDim10';
loadResult = true;
% Parse optional agruments
assignopts(who,varargin);

% Load the data
switch dataType
    case('traj')
        data = TrajectoryData().load(dir_list.trajectory{idx}); 
        if preprocess
            [data,~] = util.preprocessGridTaskTrajData(data,'centerPos',centerPos);
        end
        decoderName = data(1).decoderName;
    case('trial')
        data = Trial().load(dir_list.trial{idx});
        decoderName = data(1).Data.decoder.parameters.name;
    case('data')
        data = load(dir_list.translated{idx});
        data = data.Data;
        decoderName = data(1).TrialData.Decoder.Parameters.name;
    case('export')
        data = load(dir_list.export{idx});
        data = data.S;
        decoderName = S(1).decoderName;
end
% Load the decoder
P = load(fullfile(dir_list.base,decoderName));
P = P.bci_params;

% Load the results
if loadResult
    if isfield(P,'ExpInfo')
        runIdx = P.ExpInfo.gpfaRunID;
    else
        runIdx = P.decoderID;
    end
    runDir = sprintf('mat_results/run%03d', runIdx);
    result = load(fullfile(dir_list.base,'analysis',runDir,gpfaName));
else
    result = [];
end