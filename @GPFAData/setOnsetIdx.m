% setOnsetIdx       Set trajectory onset index for GPFAData class
%
% Usage:
%   GP = GP.setOnsetIdx(tOnset)
%
% This function sets the onset index for the neural trajectory in G based
% on the provided onset time tOnset.
%
% Inputs:
%   GP             GPFA object
%   tOnset         Time Onset
%
% Optional Arguments:
%   lag            Define the lag between GP onset and tOnset. Helpful if
%                       you want to include the bin before tOnset or 
%                       exclude the 1st bin.
%
% Author:   Alan Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function GP = setOnsetIdx(GP,tOnset,varargin)

% Parse optional arguments
lag = 0;

assignopts(who,varargin);

% Apply specified onset lag
GP.onset = tOnset + lag;

% Calculate corresponding onset index.  GPFA data should be binned starting
% at the start of the trial, with t=1 corresponding to the first time step
% in the trial.  Thus, the time at the end of the bins will correspond to
% multiples of the binWidth.
GP.onsetIdx = ceil(GP.onset/GP.binWidth);