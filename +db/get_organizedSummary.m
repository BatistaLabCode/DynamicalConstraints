function [orgD] = get_organizedSummary(D,varargin)
%% Query the information from the dataset based on the provided information
% to reduce the search results, if no information is provided than simply
% provide a breakdown by subjects, then task, and then date. The data can
% also be organized according to 

% Optional inputs
Subject = cell(0);    % List the subjects that you want to include in the 
                      % summary dataset.
Tasks = cell(0);      % List the task that you want to include in the 
                      % summary dataset.
condDate = [];        % Condition input in case we are looking for a 
                      % specific period of information, ie all data from
                      % October 2018 or all data from 2019. Input as a
                      % cell.***Remember to add the AND/OR property to the
                      % previous state.***
orgStyle = 'Subject'; % How to organize the data, default is organized 
                      % first by subject and then task.   

% Parse optional agruments
assignopts(who,varargin);

switch orgStyle
    case{'Subject'}
        % If subject is empty then create a list for every subject.
        if isempty(Subject)
            Subject = unique([{D.subject}]);
        end
        
        nSub = length(Subject);
        
        % Iterate through each subject.
        orgD = [];
        for n = 1:nSub
            % Mask the relevant data
            sMask = ismember([{D.subject}],Subject(n));
            tD = D(sMask); % Create a temporary D structure for the info.
            
            % Iterate through the list of tasks, if the task variable is
            % empty then create a list of all possible tasks.
            if isempty(Tasks)
                tTasks = unique([{tD.task}]);
            else
                tTasks = Tasks;
            end
            nTask = length(tTasks);
            
            tOrgD = [];
            for m = 1:nTask
                % Mask the relevant task data
                taskMask = ismember([{tD.task}],tTasks(m));
                
                % Check if there is any conditional requirement for the date
                if ~isempty(condDate)
                    % Convert the dataset information to numbers.
                    dateset = cell2mat([{tD.dataset}]');
                    dateset = str2num(dateset);
                    
                    % Determine how many date conditional there are and
                    % create a string
                    evalStr = 'maskDate =';
                    for k = 1:length(condDate)
                            evalStr = [evalStr ' dateset' condDate{k}];
                    end
                    
                    % Run the evaluation
                    eval([evalStr ';'])
                    
                    % Update the mask
                    taskMask = taskMask & maskDate';
                end
                
                tOrgD(m).subject = Subject{n};
                tOrgD(m).task = tTasks{m};
                tOrgD(m).numSessions = sum(taskMask);
                tOrgD(m).sessions = tD(taskMask);
            end
            
            % Add to the final D structure
            orgD = [orgD tOrgD];
        end
    case{'Task'}
        % If subject is empty then create a list for every subject.
        if isempty(Tasks)
            Tasks = unique([{D.task}]);
        end
        
        nTasks = length(Tasks);
        
        % Iterate through each task.
        orgD = [];
        for n = 1:nTasks
            % Mask the relevant data
            tMask = ismember([{D.task}],Tasks(n));
            tD = D(tMask); % Create a temporary D structure for the info.
            
            % Iterate through the list of subjects, if the subject variable
            % is empty then create a list of all possible tasks.
            if isempty(Subject)
                tSubject = unique([{tD.subject}]);
            else
                tSubject = Subject;
            end
            
            nSubject = length(tSubject);
            
            tOrgD = [];
            for m = 1:nSubject
                % Mask the relevant task data
                taskMask = ismember([{tD.subject}],tSubject(m));
                
                % Check if there is any conditional requirement for the date
                if ~isempty(condDate)
                    % Convert the dataset information to numbers.
                    dateset = cell2mat([{tD.dataset}]');
                    dateset = str2num(dateset);
                    
                    % Determine how many date conditional there are and
                    % create a string
                    evalStr = 'maskDate =';
                    for k = 1:length(condDate)
                            evalStr = [evalStr ' dateset' condDate{k}];
                    end
                    
                    % Run the evaluation
                    eval([evalStr ';'])
                    
                    % Update the mask
                    taskMask = taskMask & maskDate';
                end
                
                tOrgD(m).task = Tasks{n};
                tOrgD(m).subject = tSubject{m};
                tOrgD(m).numSessions = sum(taskMask);
                tOrgD(m).sessions = tD(taskMask);
            end
            
            % Add to the final D structure
            orgD = [orgD tOrgD];
        end
    otherwise
        warning('Could not find the defined organization method, defaulting to SUBJECT organization')
        orgD = get_organizedSummary(D,'Subject',Subject,'Tasks',Tasks,...
            'condDate',condDate)
end