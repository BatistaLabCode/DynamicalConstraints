% setOnsetIdx       Set trajectory onset index for GPFAData class
%
% Usage:
% G = setOnsetIdx(G,tOnset)
%
% This function sets the onset index for the neural trajectory in G based
% on the provided onset time tOnset.
%
% Author:   Alan Degenhart
% Created:  2018.10.25

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