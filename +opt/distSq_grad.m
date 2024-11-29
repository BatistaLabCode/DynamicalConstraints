function dJ = distSq_grad(M, P, obj_params)
% Gradient function
%
% Inputs
%   M       Current value of variable being optimized over
%   P       Structure containing necessary data
%   obj_params  Structure of objective function parameters
%
% Parameters used from P: mu_A, mu_B, mu_AB, sig_AB, mu_BA, sig_BA
%
% Outputs:
%   dJ       Gradient of objective function
%
% Created by Alan Degenhart and Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

if isempty(obj_params)
    obj_params.w_mid = 1;
    obj_params.w_var = 1;
    obj_params.w_start = 1;
end

w_mid = obj_params.w_mid;  % Weight for trajectory midpoint term
w_var = obj_params.w_var;  % Weight for midpoint variance term
w_start = obj_params.w_start;  % Weight for start point term

p1 = M(:,1);
p2 = M(:,2);

% Unpack parameters for clarity
mu_A = P.mu_A;
mu_B = P.mu_B;
mu_AB = P.mu_AB;
mu_BA = P.mu_BA;
sig_AB = P.sig_AB;
sig_BA = P.sig_BA;

% Compute individual parts of objective function
part1 = w_mid * (p1' * (mu_AB - mu_BA) * (mu_AB - mu_BA)')';  % Maximize distance between centers of midpoints
part2 = - w_var * (2*sig_AB + 2*sig_BA) * p1;  % Minimize variance of projection @ midpoint of trajectory
part3 = w_start * (p2' * (mu_A - mu_B) * (mu_A - mu_B)')';  % Maximize distance between starting points of trajectories

d_p1 = -(part1 + part2);
d_p2 = -part3;
dJ = [d_p1 d_p2];