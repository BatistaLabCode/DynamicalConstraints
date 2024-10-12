function TD = setKinematicData(TD,KD)
% [KD] = getKinematicData(TD,kinSource)
%
% Set kinematic data method for TrajectoryData class
%
% Author:  Alan D. Degenhart
% Date Created: 2016/06/21
% Last Updated: 2016/06/21
% Last Update:  Initial version of function

% Make sure that the kinematic data is the appropriate size.
assert(length(TD) == length(KD), ...
    'TrajectoryData object and KinematicData object must be the same length.')

% Loop over all TrajectoryData elements and set KinematicData
nTrials = length(TD);
for i = 1:nTrials
    % Set kinematic data.  This can be either the hand position or the
    % brain-controlled cursor position.  The 'source' field in the
    % kinematic data object specifies this.
    
    kinSource = KD(i).source;
    srcStr = [kinSource 'Kin'];
    TD(i).(srcStr) = KD(i);
end