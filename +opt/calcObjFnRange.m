function J = calcObjFnRange(P, obj_params, varargin)
% Uses the calculated optimized data, simplified for offline.
% 
% Optional Inputs:
%   N       Number of iterations
%   objFn   Objective function to iterate over
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

N = 10000;
objFn = 'standard';

assignopts(who, varargin);

% Define dimensionality
R = 2;      % Number of orthonormal dimensions
D = P.D;    % Latent state dimensionality

% Loop over iterations
J = nan(1, N);
for i = 1:N
    % Generate random orthonormal basis
    A = randn(D,R);
    M = orth(A);    % Columns of M will be orthonormal
    
    % Find value of objective function
    switch objFn
        case 'standard'
            J(i) = standard_obj(M, P, obj_params);
        case 'distSq'
            J(i) = distSq_obj(M, P, obj_params);
    end
end
