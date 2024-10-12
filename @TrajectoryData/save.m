function save(TD,pathName,fileBase)
% Save method for TrajectoryData class.
%
% TD.save(pathName,fileBase)
%
% Inputs
%   pathName        Directory to save data in
%   fileBase        Name of dataset.
%
%
% Author:           Alan D. Degenhart
% Date Created:     N/A
% Date Updated:     2016/08/05
% Last Update:      Adding comments.

% Get data to add to dataset structure
subject = TD(1).subject;
date = TD(1).date;
nTrials = length(TD);

block = nan(1,nTrials);
successCode = nan(1,nTrials);
trialNo = nan(1,nTrials);

% Create data directory if necessary
dataDir = [pathName '/' fileBase];
mkdir(dataDir)

cmd_prog = util.commandline_progress([], 0, nTrials, ...
    'Saving TrajectoryData');
% Loop over all trials and save
for i = 1:nTrials
    util.commandline_progress(cmd_prog, i);
    
    % Get data for trial
    block(i) = TD(i).tag;
    successCode(i) = logical(TD(i).successful);
    trialNo(i) = TD(i).trialID;

    % Set filename for trial
    fName = sprintf('Trl_%0.4d.mat',trialNo(i));
    tempTrial = TD(i);
    save([dataDir '/' fName],'tempTrial')
end

% Save trial info structure
trialInfo.subject = subject;
trialInfo.date = date;
trialInfo.filename = fileBase;
trialInfo.nTrials = nTrials;
trialInfo.trialNo = trialNo;
trialInfo.block = block;
trialInfo.successCode = successCode;

save([dataDir '/' fileBase '_trialInfo.mat'],'trialInfo')