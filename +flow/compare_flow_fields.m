% flow.compare_flow_fields       Compare flow fields between conditions
%
% Usage:
%   [S] = flow.calc_session_stats(FF_base, FF_comp)
%
% This function calculates the MSE distribution when comparing the two flow
% fields. It is possible to calculate other distributions by added
% additional calculations to this function.
%
% Inputs:
%   FF_base        Flow field for one condition
%   FF_comp        Flow field for the second condition
%
% Optional Inputs:
%   cond_str       Condition comparison label
%   saveRaw        Save the original flow fields to the output structure.
%   saveSimp       Save a reduced structure that does not include grid info
% 
% Outputs:
%   A              Structure that includes the stat comparison and the
%                    number of overlapping voxels between conditions.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [A] = compare_flow_fields(FF_base, FF_comp, varargin)
cond_str = [];
saveRaw = 1;
saveSimp = 0; % Determine whether or not to save the simple matrix
assignopts(who, varargin);
% Make sure grids are the same between conditions
assert(FF_base.grid.d_grid == FF_comp.grid.d_grid, ...
    'Grids must match across conditions')

% Loop over grid
n_grid = length(FF_base.grid.grid) - 1;
mse_grid = nan(n_grid, n_grid);
for i = 1:n_grid
    for j = 1:n_grid
        % Get velocities for two conditions
        v_base = [FF_base.grid.VX(i,j), FF_base.grid.VY(i,j)];
        v_comp = [FF_comp.grid.VX(i,j), FF_comp.grid.VY(i,j)];
        
        % Compare angle between conditions.  Need to verify that both
        % velocity vectors are valid (i.e., non-nan)
        if sum(isnan([v_base, v_comp])) == 0
            % Calculate mean squared error
            mean_sq_err = mean((v_base - v_comp).^2);
            
            % Add result to grid
            mse_grid(i, j) = mean_sq_err;
        end
    end
end

% Put results in structure
valid_mask = ~isnan(mse_grid);

if saveSimp
    % Add raw data
    A.condition = cond_str;
    A.n_valid = sum(sum(valid_mask));
    
    % Mean squared error
    A.mseValid = mse_grid(valid_mask)';
    A.mseMean = mean(A.mseValid);
    A.mseMedian = median(A.mseValid);
else
    % Add raw data
    A.condition = cond_str;
    if saveRaw
        A.F_base = FF_base;
        A.F_comp = FF_comp;
    end
    A.n_valid = sum(sum(valid_mask));
    
    % Mean squared error
    A.mse.grid = mse_grid;
    A.mse.valid = mse_grid(valid_mask);
    A.mse.mean = mean(A.mse.valid);
    A.mse.median = median(A.mse.valid);
end

end