% spatialAverage        Calculate spatial average of trajectories.
%
% This function find the spatial average of the provided trajectories xMat.
% In contrast to other averaging methods (which either time-warp or assume
% that the trajectories are temporally aligned, this function assumes that
% the velocity is a function of the position.
%
% For each position in the average, the velocity of the state from that 
% position (which defines the next position in the average trajectory) is
% defined as the vector average of the velocity vectors within radiur r of
% the current point.
%
% Inputs
%   xMat    Matrix of trials to average (trials x dimensions x time points)
%   r       Radius used to identify nearby points
%
% Outputs
%   x_avg   Averaged trajectory (dimensions x time points)
%
% Created:  2019.08.01
% Author:   Alan Degenhart

function x_avg = spatialAverage(xMat, varargin)

% Optional arguments
max_iter = 100;
r_cutoff_frac = 1;
r_factor = 0.15;

assignopts(who, varargin);

% Get data size
nTrials = size(xMat,1);
nDim = size(xMat,2);
nSamp = size(xMat,3);

% Plot all trajectories (debugging purposes)
% figure; hold on;
% for i = 1:nTrials
%     x_temp = squeeze(xMat(i,:,:));
%     plot(x_temp(1,:), x_temp(2,:), 'color', 0.75 * ones(1,3))
% end

% Calculate velocity
dx = diff(xMat,1,3);

% Determine start and end points. Here we assume that the trajectories all
% start at the first provided point and end at the last non-nan point.

% Start point
x_start = squeeze(xMat(:,:,1));
x_start = nanmean(x_start,1);

% End point
x_end = nan(nTrials,nDim);
for i = 1:nTrials
    x_temp = squeeze(xMat(i,:,:));
    idx = find(~isnan(x_temp(1,:)), 1, 'last');
    if isempty(idx)
        idx = 1;
    end
    x_end(i,:) = x_temp(:,idx);
end
x_end = nanmean(x_end,1);

% plot(x_start(1),x_start(2),'ko','MarkerFaceColor','g')
% plot(x_end(1),x_end(2),'ko','MarkerFaceColor','r')

% Determine distance as a factor of the target-to-target distance
r = norm(x_start - x_end)*r_factor;

% Start looping over time points and calculating average
x_avg = nan(nDim, 1000);  % Just something large for now
x_avg(:,1) = x_start;
idx = 1;
stop_flag = false;
while ~stop_flag
    % Find all points that are with radius r of the current point
    mask = false(nTrials,nSamp);
    for i = 1:nTrials
        d = abs(squeeze(xMat(i,:,:)) - x_avg(:,idx));
        d = vecnorm(d);
        mask(i,:) = d < r;
        
        % Set the value for the last non-zero element in d equal to zero,
        % as there is no derivative for this point.
        mask(i,find(~isnan(d),1,'last')) = false;
    end
    
    n_valid_samp = sum(sum(mask));
    dx_valid = nan(nDim,n_valid_samp);
    idx2 = 1;
    for i = 1:nTrials
        dx_temp = squeeze(dx(i,:,mask(i,:)));
        
        % Add data to matrix
        if ~isempty(dx_temp)
            if size(dx_temp,1) == 1; dx_temp = dx_temp'; end
            inds = idx2:(idx2 + size(dx_temp,2) - 1);
            dx_valid(:,inds) = dx_temp;
            idx2 = inds(end) + 1;
        end
    end
    
    % Average velocity and integrate to find next position
    dx_valid = nanmean(dx_valid,2);
    x_temp = x_avg(:,idx) + dx_valid;
    x_avg(:,idx+1) = x_temp;
    idx = idx + 1;
    
    % Determine stopping criteria
    if (norm(x_temp - x_end') < r * r_cutoff_frac) || idx > max_iter
        stop_flag = true;
        x_avg(:,idx+1) = x_end';
    end
end

% Remove NaNs
x_avg = x_avg(:,1:idx+1);

% plot(x_avg(1,:),x_avg(2,:),'k.-')