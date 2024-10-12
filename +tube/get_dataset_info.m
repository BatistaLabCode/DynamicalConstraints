% tube.get_dataset_info  Get dataset paths for constrained path (tube)
% experiment. Creates and propogrates IntTargExp object (IT).
%
% Usage:
%   db.get_dataset_info(subject,dataset)
%
% 
% This will automatically fill in the intermediate target experiment
% information based on the database structure. 
%
% Optional Input:
%   D           Structure with all the valid experimental data. Save as
%                   filename <<publicationQualitySessions.mat>>
%   inclTT      Determines if two-target data should be saved to IT
%
% Optional Input:
%   IT          IntTargExp object
%   dir_tD      Directory of the data
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [IT,dir_tD] = get_dataset_info(subject,dataset,varargin)

D = [];
inclTT = true;
assignopts(who,varargin);

% Initialize intermediate target experiment object
IT = IntTargExp();
IT.subject = subject;
IT.date = dataset;

% Load the D structure if not defined.
if isempty(D)
    D = el.db.load_info;
end

% Mask the data for the relevant sessions.
mask = ismember([{D.subject}],subject) & ismember([{D.dataset}],dataset);
tD = D(mask);
dir_tD = db.get_dataset_dirs(tD);
task = [{tD.task}];

% Iterate through the task and conditions
n_int = 0;
n_rot = 0;
for n = 1:size(task,2)
    
    for m = 1:size(dir_tD(n).trajectory,2)
        % Find the file name
        tempFile = dir_tD(n).trajectory{m};
        if isempty(tempFile)
            continue
        end
        pos = regexp(tempFile,'\');
        tempFile = tempFile(pos(end)+1:end);
        
        % Load the information into the structure
        switch task{n}
            case {'inttarg_rot_unconstr',...
                    'perturb_hc_unconstr',...
                    'perturb_bc_unconstr'}
                IT.unconstrainedBlockDir = tempFile;
            case {'inttarg_rot_constr_sing',...
                    'perturb_hc_constr_sing',...
                    'perturb_bc_constr_sing'}
            case {'inttarg_rot_constr_slow',...
                    'inttarg_rot_constr_fast',...
                    'tt_int_constr_slow',...
                    'perturb_hc_constr_slow',...
                    'perturb_bc_constr_slow'}
                IT.constrainedBlockDir{m} = tempFile;
            case {'inttarg_rot_unconstr_wo',...
                    'perturb_hc_unconstr_wo',...
                    'perturb_bc_unconstr_wo'}
                IT.unconstrainedWashoutBlockDir = tempFile;
            case {'tt_int'}
              if inclTT
                  n_int = n_int + 1;
                  IT.ttIntBlockDir{n_int} = tempFile;
              end
            case {'tt_rot',...
                    'tt_rot_comp',...
                    'tt_rot_inv'}
               if inclTT
                   n_rot = n_rot + 1;
                   IT.ttRotBlockDir{n_rot} = tempFile;
               end
        end
    end
end

if isempty(IT.unconstrainedBlockDir) & ~isempty(IT.constrainedBlockDir)
    IT.unconstrainedBlockDir = IT.constrainedBlockDir{1};
end
end
