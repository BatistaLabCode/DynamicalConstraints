% db.get_dataset_dirs  Get set of directories for selected dataset
%
% Usage:
%   db.get_dataset_dirs(D)
%
% Inputs:
%   D           Structure with all the valid experimental data. Save as
%                   filename <<publicationQualitySessions.mat>>
% 
% Outputs:
%   dir_struct  Structure with the location information of all files.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function dir_struct = get_dataset_dirs(D,varargin)
showWarning = 1; % This will determine whether or not to plot the warnings,
                 % default is to include the warnings, however it is useful
                 % to hide warning when loading examples sessions that are
                 % isolated from the full data pack (ie the code pack).
assignopts(who,varargin);

% Get size of the input dataset structure.  Currently this function
% supports up to 2D structure arrays.
sz = size(D);
assert(numel(sz) <= 2, 'This function only supports 1D and 2D structure arrays as input.')

% Define path to data
[~, loc_path] = serverPath;

% Initialize directory structure.  This will contain the paths to all types
% of data (translated, trialData, trajectoryData, exportData)
dir_struct = repmat(struct('subject',[],...
    'dataset',[],...
    'n_dir',[],...
    'base',[],...
    'translated',[],...
    'trajectory',[]),sz);

% Iterate over elements of D and get directories
for i = 1:sz(1)
    for j = 1:sz(2)
        % Get dataset info
        D_temp = D(i, j);
        dir_list = D_temp.directory_name;
        n_dir = length(dir_list);
        dir_struct(i, j).subject = D_temp.subject;
        dir_struct(i, j).dataset = D_temp.dataset;
        dir_struct(i, j).n_dir = n_dir;

        % Define base path to data
        base_path = fullfile( ...
            loc_path, D_temp.subject, D_temp.dataset(1:4), ...
            D_temp.dataset(5:6), D_temp.dataset);
        dir_struct(i, j).base = base_path;

        % Loop over directories and get paths to each type of data
        % (translated and trajectory)
        translated_dir_list = cell(1, n_dir);
        trajectory_dir_list = cell(1, n_dir);
        for k = 1:n_dir
            % Get translated data paths
            dir_base = sprintf('%s%s_%0.2d_', D_temp.subject, ...
                D_temp.dataset, dir_list(k));
            trans_dir = fullfile(base_path, 'translated');
            dir_contents = util.findDirContents(trans_dir, dir_base);
            data_type = 'translated';
            if length(dir_contents) > 1
                % If more than one file is able, make sure to use the
                % dataset with the snippet info and not in us
                maskSI = contains(dir_contents,'SI');
                maskT = ~contains(dir_contents,'us');
                if sum(maskSI & maskT)==1
                    fname = fullfile(trans_dir, dir_contents{(maskSI & maskT)==1});
                elseif sum(maskSI & maskT)>1
                    pos = find(maskSI & maskT,1,'first');
                    if showWarning
                        warning(['More than one valid dataset/task combination ' ...
                            'found for dataset %s%s task %s (%s).  Using the ' ...
                            'first one with snippet info in ms.\n'], ...
                            D_temp.subject, D_temp.dataset, D_temp.task, data_type)
                    end
                    fname = fullfile(trans_dir, dir_contents{pos});
                else
                    if showWarning
                        warning(['More than one valid dataset/task combination ' ...
                            'found for dataset %s%s task %s (%s).  Using the first one.\n'], ...
                            D_temp.subject, D_temp.dataset, D_temp.task, data_type)
                    end
                    fname = fullfile(trans_dir, dir_contents{1});
                end
            elseif isempty(dir_contents)
                if showWarning
                    warning(['No valid dataset/task combinations found ' ...
                        'for dataset %s%s task %s.\n'], ...
                        D_temp.subject, D_temp.dataset, D_temp.task)
                end
                fname = [];
            else
                fname = fullfile(trans_dir, dir_contents{1});
            end
            translated_dir_list{k} = fname;

            % Get trajectory data paths
            traj_dir = fullfile(trans_dir, 'trajectoryData');
            dir_contents = util.findDirContents(traj_dir, dir_base);
            data_type = 'trajectoryData';
            if length(dir_contents) > 1
                % If more than one file is able, make sure to use the
                % dataset with the snippet info and not in us
                maskSI = contains(dir_contents,'SI');
                maskT = ~contains(dir_contents,'us');
                if sum(maskSI & maskT)==1
                    fname = fullfile(traj_dir, dir_contents{(maskSI & maskT)==1});
                elseif sum(maskSI & maskT)>1
                    pos = find(maskSI & maskT,1,'first');
                    if showWarning
                    warning(['More than one valid dataset/task combination ' ...
                        'found for dataset %s%s task %s (%s).  Using the ' ...
                        'first one with snippet info in ms.\n'], ...
                        D_temp.subject, D_temp.dataset, D_temp.task, data_type)
                    end
                    fname = fullfile(traj_dir, dir_contents{pos});
                else
                    if showWarning
                    warning(['More than one valid dataset/task combination ' ...
                        'found for dataset %s%s task %s (%s).  Using the first one.\n'], ...
                        D_temp.subject, D_temp.dataset, D_temp.task, data_type)
                    end
                    fname = fullfile(traj_dir, dir_contents{1});
                end
            elseif isempty(dir_contents)
                if showWarning
                warning(['No valid dataset/task combinations found ' ...
                    'for dataset %s%s task %s.\n'], ...
                    D_temp.subject, D_temp.dataset, D_temp.task)
                end
                fname = [];
            else
                fname = fullfile(traj_dir, dir_contents{1});
            end
            trajectory_dir_list{k} = fname;

        end

        dir_struct(i, j).translated = translated_dir_list;
        dir_struct(i, j).trajectory = trajectory_dir_list;
    end
end