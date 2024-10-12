% el.flow.batch_createIT  Batch analysis function creating the ID data 
% structure


function invalid_list = batch_createIT(varargin)
save_path = [];
D = [];
task = 'rot_constr';

optArg = assignopts(who,varargin);

if isempty(save_path)
    save_path = uigetdir;
end

if isempty(D)
    load('D:\Figures\EL_NatNeuro_2024_data\publicationQualitySessions.mat')
end

% Get all the tube sessions
maskTask = contains([{D.task}],task);
tD = D(maskTask);
[dataset,idx] = unique([{tD.dataset}]);
subject = {tD(idx).subject};
ds_info = [subject' dataset'];
n_ds = size(ds_info, 1);

% Iterate over experiments and run analysis
invalid_mask = false(1, n_ds);
invalid_list = cell(1, n_ds);
for i = 1:n_ds
    try
        % Run analysis
        fprintf('Processing dataset %s %s ... ', ...
            ds_info{i, 1}, ds_info{i, 2})
        
        if ismember(ds_info{i,1},'Quincy')
            if str2num(ds_info{i,2})>20210101
                centerPos = [-65 -330 0];
            else
                centerPos = [60,-270,0];
            end
        else
            centerPos = [0 0 0];
        end
        % Get IT object for dataset
        [IT] = tube.get_dataset_info(ds_info{i, 1}, ds_info{i, 2},'D',D);
        
        % Load data
        dir_list = db.get_dataset_dirs(D(ismember({D.dataset},dataset(i))))
        [IT] = tube.get_IT_data(IT,'centerPos',centerPos,...
            'dir_TD',dir_list(1));
        
        % Process the data
        [IT, T] = tube.constrained_path_analysis(IT,optArg{:});
        if ~exist(fullfile(save_path,'success_rate_data'))
            mkdir(fullfile(save_path,'success_rate_data'))
        end
        T_name = [IT.subject IT.date '_suc'];
        
        
        % Save results
        [dir_info] = tube.get_dir_info(IT.subject, IT.date,...
            'save_path',save_path);
        f_name_int_targ = [ds_info{i, 1}, ds_info{i, 2}, '_int_targ.mat'];
        save(fullfile(dir_info.int_targ_data_path, f_name_int_targ), 'IT')
        save(fullfile(dir_info.suc_data_path,T_name),'T')

        fprintf('done.\n')
    catch ERR
        fprintf('error encountered.  Skipping ... \n')
        invalid_mask(i) = true;
        invalid_list{i} = sprintf('%s %s', ds_info{i, 1}, ds_info{i, 2});
    end
end

% List any datasets that were skipped due to errors
n_invalid = sum(invalid_mask);
if n_invalid > 0
    fprintf('\nThe following datasets were skipped due to errors:\n\n')
    invalid_list = invalid_list(invalid_mask);
    for i = 1:n_invalid
        fprintf('%s\n', invalid_list{i})
    end
end