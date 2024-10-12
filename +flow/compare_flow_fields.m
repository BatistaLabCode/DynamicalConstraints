% compare_flow_fields       Compare flow fields between conditions
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
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
ang_diff_grid = nan(n_grid, n_grid);
mag_diff_grid = nan(n_grid, n_grid);
mse_grid = nan(n_grid, n_grid);
for i = 1:n_grid
    for j = 1:n_grid
        % Get velocities for two conditions
        v_base = [FF_base.grid.VX(i,j), FF_base.grid.VY(i,j)];
        v_comp = [FF_comp.grid.VX(i,j), FF_comp.grid.VY(i,j)];
        
        % Compare angle between conditions.  Need to verify that both
        % velocity vectors are valid (i.e., non-nan)
        if sum(isnan([v_base, v_comp])) == 0
            ang_base = atan2d(v_base(2), v_base(1));
            ang_comp = atan2d(v_comp(2), v_comp(1));
            ang_diff = abs(ang_base - ang_comp);
            
            % Correct for angles > 180
            if ang_diff > 180; ang_diff = ang_diff - 180; end
            
            % Calculate difference in magnitude
            mag_diff = abs(norm(v_base) - norm(v_comp));
            mean_sq_err = mean((v_base - v_comp).^2);
            
            % Add result to grid
            ang_diff_grid(i, j) = ang_diff;
            mag_diff_grid(i, j) = mag_diff;
            mse_grid(i, j) = mean_sq_err;
        end
    end
end

% Put results in structure
valid_mask = ~isnan(ang_diff_grid);

if saveSimp
    % Add raw data
    A.condition = cond_str;
    A.n_valid = sum(sum(valid_mask));
    
    % Angular error
    A.angValid = ang_diff_grid(valid_mask)';
    A.angMean = mean(A.angValid);
    A.angMedian = median(A.angValid);
    
    % Magnitude error
    A.magValid = mag_diff_grid(valid_mask)';
    A.magMean = mean(A.magValid);
    A.magMedian = median(A.magValid);
    
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
    
    % Angular error
    A.ang.grid = ang_diff_grid;
    A.ang.valid = ang_diff_grid(valid_mask);
    A.ang.mean = mean(A.ang.valid);
    A.ang.median = median(A.ang.valid);
    
    % Magnitude error
    A.mag.grid = mag_diff_grid;
    A.mag.valid = mag_diff_grid(valid_mask);
    A.mag.mean = mean(A.mag.valid);
    A.mag.median = median(A.mag.valid);
    
    % Mean squared error
    A.mse.grid = mse_grid;
    A.mse.valid = mse_grid(valid_mask);
    A.mse.mean = mean(A.mse.valid);
    A.mse.median = median(A.mse.valid);
end

end