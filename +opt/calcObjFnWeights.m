function [objFnParams] = calcObjFnWeights(P,varargin)
% Uses the calculated optimized data, simplified for offline.
% Calculate appropriate scaling weights.
% 
% Usage:
%   [objFnParams] = opt.calcObjFnWeights(P)
%
% Inputs:
%   P   Objective function to iterate over
%
% Optional Inputs:
%   objFn   Objective function to iterate over
%   q       Qunatile range
%
% Output:
%   objFnParams Structure that contains the scaling weights.
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart 
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

objFn = 'standard';
q = [0, 1];  % Define quantile range

assignopts(who,varargin);

objFnParams.w_mid = 1;
objFnParams.w_var = 0;
objFnParams.w_start = 0;

J_m = opt.calcObjFnRange(P, objFnParams, 'objFn', objFn);
q_m = quantile(J_m, q);
w_m = 1/abs(diff(q_m));

objFnParams.w_mid = 0;
objFnParams.w_var = 1;
objFnParams.w_start = 0;

J_v = opt.calcObjFnRange(P, objFnParams, 'objFn', objFn);
q_v = quantile(J_v, q);
w_v = 1/abs(diff(q_v));

objFnParams.w_mid = 0;
objFnParams.w_var = 0;
objFnParams.w_start = 1;

J_s = opt.calcObjFnRange(P, objFnParams, 'objFn', objFn);
q_s = quantile(J_s, q);
w_s = 1/abs(diff(q_s));

objFnParams.w_mid = w_m;
objFnParams.w_var = w_v;
objFnParams.w_start = w_s;
end