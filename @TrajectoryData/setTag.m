function TD = setTag(TD,tagVal)
% Set tag for TrajectoryData object
%
% This function sets the tag field for each trial in TrajectoryData object
% TD to the value specified by TAGVAL.

% Loop over trials and set tag value
nTrials = length(TD);
for i = 1:nTrials
    TD(i).tag = tagVal;
end