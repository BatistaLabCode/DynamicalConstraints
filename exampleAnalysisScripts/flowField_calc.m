% This example script will run through how to calculate the flow fields for
% an example session 20190719. The end result will be the same as the flow
% field saved in the folder <flowAnalysis/mat>. This is the data that is
% used in figure 5. The final result will be the structure <FlowResults>.
% The user also has the option to create flow field summary plots for each
% flow field condition. The user has the option to save thedata and the
% figures as well.
%
% Created by Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

clear, close all % Clear and set the workspace

% Add the paper code to the path
pathName = pwd;
addpath(genpath(pathName))

% Set the example session and if you want to save the data
exampleSess = {'20190719'};
plotFig = 1;
saveData = 0;         % Determine if the data is saved or not.
saveFig = 0;          % Determine if the figures are saved or not.
saveFigPath = [];     % Where to save the figure
saveDataPath = [];    % where to save the data

if saveData
    saveDataPath = uigetdir('',"Where do you want to save the flow field data?");
end

% Define the colormaps
C_int = util.defineTaskColormap('bc_int');
C_rot = util.defineTaskColormap('bc_rot');

% Load in the D structure
dataLoc = serverPath;
D = load(fullfile(dataLoc,'exampleDatasetCatalog.mat'));
D = D.D;

% Determine the session with the correct data
D = D(ismember({D.dataset},exampleSess));
dir_list = db.get_task_datasets(D, {'tt_int','tt_rot'});

FlowResults = flow.main_analysis(dir_list(1,1).base,...
    dir_list(1,1).trajectory{1},dir_list(1,2).trajectory{1},...
    'data_save_loc',saveDataPath,'saveData', saveData);

% Plot the flow field conditions
if plotFig
    F(1) = flow.plot_flow_field(FlowResults.FF_int,'C',C_int,'arrow_size',7);
    F(2) = flow.plot_flow_field(FlowResults.FF_pred,'C',C_int,'arrow_size',7);
    F(3) = flow.plot_flow_field(FlowResults.FF_rot,'C',C_rot,'arrow_size',7);
    F(4) = flow.plot_flow_field(FlowResults.FF_pred_Int,'C',C_rot,'arrow_size',7);
end

% Save the d' data structure to a single mat file.
if ~isempty(saveDataPath)
    fileName = sprintf('FlowField_example%s.mat',exampleSess{:});
    save(fullfile(saveDataPath,fileName),'FlowResults')
end

if saveFig & plotFig
    if ~isempty(saveFigPath)
        saveFigurePDF(F,saveFigPath)
    else
        saveFigurePDF(F)
    end
end

