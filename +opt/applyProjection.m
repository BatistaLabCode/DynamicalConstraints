function TD = applyProjection(TD,P,varargin)
% TD = applyProjection(TD,P,varargin)
%
% Apply projection/transform to GPFA data within the TD. Modification of
% the previous function <<applyDataHighProjection>>, so that we don't need 
% to create a large G structure offline.
%
% Usage:
%   [TD] = opt.applyProjection(TD,P)
%
% Inputs:
%   TD      TrajectoryData data structure
%   P       Projection matrix (nDim x 2)
%
% Optional Inputs:
%   o_pre   Offset vector (pre-transformation, applied to latents)
%   o_post  Offset vector (post-transformation. applied to position)
%   xSpec   Latent state for projection to non-orthonormalized (xsm) or 
%           orthonormalized (xorth, default)
%   savePos Saves the orthonormalized projection to the cursor position.
%
% Outputs:
%   TD      TrajectoryData data structure with the projection applied
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Parse optional arguments
o_pre = [];     % Offset vector (pre-transformation, applied to latents
o_post = [];    % Offset vector (post-transformation. applied to position)
xSpec = 'xorth';% Latent State Space
savePos = 0;    % Save the projection to the cursor position.

assignopts(who,varargin);

% If no offset is provided, specify zero-offset vector
M = size(TD(1).GPFA.(xSpec),1);
D = size(P,2);
if isempty(o_pre)
    o_pre = zeros(M,1);
end
if isempty(o_post)
    o_post = zeros(D,1);
end

% Loop over all trajectories
for i = 1:length(TD)
    X = TD(i).GPFA.(xSpec);
    nSamp = size(X,2);
    TD(i).GPFA.(xSpec) = P'*(X + repmat(o_pre,1,nSamp)) ...
        + repmat(o_post,1,nSamp);
    if savePos
        TD(i).brainKin.pos = TD(i).GPFA.(xSpec)';
    end
end

end