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
        % Add the index data to the structure and its original marker
        resampIndices(n,k).idx = temp([pos(k)+1]:pos(k+1));
        resampIndices(n,k).marker = marker(temp([pos(k)+1]:pos(k+1)));  
    end
end