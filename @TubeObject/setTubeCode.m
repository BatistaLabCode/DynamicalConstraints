function TU = setTubeCode(TU)
% Find unique tube name/rotation combinations

% Get tube name and rotation info
name = {TU.name};
rot = [TU.rotation];

% Find unique name/rotation combinations
[~,~,nameCode] = unique(name);
tubeCond = [nameCode rot'];
[~,~,tubeCode] = unique(tubeCond,'rows');

% Set tube code for all elements
for i = 1:length(TU)
    TU(i).tubeCode = tubeCode(i);
end