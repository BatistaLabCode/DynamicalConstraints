% Creates a structure of randomized indices that can be used for any
% bootstrapping analysis of the data.
%
% [resampIndices] = getPermutationSamples(nSz,varargin)
%
% This function returns a structure that will be N x G in size, where N is 
% number of iterations and G is number of groups. This structure will have 
% two fields:
%   1) idx: Indices of the found for the given iteration
%   2) marker: The original condition ID of the data being randomized. 
%
% Usage:
%   util.getPermutationSamples(nSz)
% 
% Inputs:
%   nSz         The size of the data structure to find boopstrapped indices             
% 
% Optional Inputs:
%   nSz         The size of the data structure to find boopstrapped indices 
%   numIter     Number of iteration to produce indices for.
%   numTrl      Number of trials of each iteration of random sampling
%   wReplace    Determine whether to randomly sample with replacement
%   groups      Number of groups to collect random samples of
%   fixRand     Fix the randomization    
%   marker      Condition ID of the data we are randomizing
%  
% Outputs:
%   resampIndices     Structure of random sample indices for every
%                       iteration and group.
%
% Author:   Alan Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [resampIndices] = getPermutationSamples(nSz,varargin)
%% Create a bootstrap randomization function of the gpfa data
numIter = 500;
numTrl = 25;
wReplace = false;
groups = 2;
fixRand = 1;
marker = [];

% Parse optional agruments
assignopts(who,varargin);

if isempty(marker)
    marker = ones(nSz,1);
end

% Set the random stat so that we always getting the same initial projections
if fixRand
    rand('state',0)
end

% create a structure for the random samples
resampIndices = repmat(struct('idx',[],'marker',[]),numIter,groups);
for n = 1:numIter
    % Randomly sample the data.
    if wReplace
        temp = randsample(nSz,groups*numTrl,true);
    elseif groups*numTrl <= nSz
        temp = randsample(nSz,groups*numTrl,false);
    else
        warning(['Number of trials per group greater than the total'...
            ' number of trials. Running with replacement.\n'])
        temp = randsample(nSz,groups*numTrl,true);
    end
    
    pos = linspace(0,groups*numTrl,groups+1);
    
    % Iterater through each group and add it to the structure
    for k = 1:groups
        % Add the index data to the structure and its original cond. marker
        resampIndices(n,k).idx = temp([pos(k)+1]:pos(k+1));
        resampIndices(n,k).marker = marker(temp([pos(k)+1]:pos(k+1)));  
    end
end