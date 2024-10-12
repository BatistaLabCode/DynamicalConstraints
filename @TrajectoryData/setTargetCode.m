function TD = setTargetCode(TD,tC)
% setTargetCode         Determine unique target code from trajectory data
%
% Author:  Alan D. Degenhart
% Date Created: 2015/01/14
% Last Updated: 2015/01/14

% Get number of trials
nTrials = length(TD);

% Check size of target code.  If is a scalar, set each element of TD to the
% value specified by tC.  If tC is the same size as TD, set each element of
% TD to the corresponding value in tC.
if length(tC) == 1
    tC = repmat(tC,1,nTrials);
end

% If tC is not a scalar, make sure the size matches that of the trajectory
% data.
assert(length(TD) == length(tC), ...
    'Size of target codes do not match that of the provided trajectory data.')

for i = 1:nTrials
    TD(i).targetCode = tC(i);
end