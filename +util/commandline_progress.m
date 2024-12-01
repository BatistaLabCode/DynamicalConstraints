% util.commandline_progress  Command-line progress bar
%
% Usage:
%   util.command_progress(cmd_prog, iter, num_iter, msg)
%
% Prints the progress of processing the iterations an analysis function.
% 
% Inputs:
%   cmd_prog        Structure for the progress bar.       
%   iter            Current iteration
%   num_iter        Total number of iterations
%   msg             Display message at the start of the progress bar.
% 
% Outputs:
%   cmd_prog        Structure for the progress bar.
%
% Author:   Alan D. Degenhart
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function [cmd_prog] = commandline_progress(cmd_prog, iter, num_iter, msg)

% If the iteration is 0, set up the progress bar
if iter == 0
    all_iter = 1:num_iter;
    all_prog = 100 * all_iter/num_iter;
    prog_int = 10:10:100;

    % Find the indices for each progress update
    n_updates = length(prog_int);
    prog_iter = nan(1, n_updates);
    for i = 1:n_updates
        [~, prog_iter(i)] = min(abs(all_prog - prog_int(i)));
    end
    
    % Display message
    fprintf('%s: ', msg)
    
    % Set up output structure
    cmd_prog.prog_int = prog_int;
    cmd_prog.prog_iter = prog_iter;
else
    % Update progress
    prog_mask = ismember(cmd_prog.prog_iter, iter);
    if sum(prog_mask) == 1
        fprintf(' ... %d%%', cmd_prog.prog_int(prog_mask))
    end
end