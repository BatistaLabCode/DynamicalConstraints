function [Xorth, Lorth, TT, DD] = orthogonalize(X, L)
%
% [Xorth, Lorth, TT] = orthogonalize(X, L)
%
% Orthonormalize the columns of the loading matrix and 
% apply the corresponding linear transform to the latent variables.
%
%   yDim: data dimensionality
%   xDim: latent dimensionality
%
% INPUTS:
%
% X        - latent variables (xDim x T)
% L        - loading matrix (yDim x xDim)
%
% OUTPUTS:
%
% Xorth    - orthonormalized latent variables (xDim x T)
% Lorth    - orthonormalized loading matrix (yDim x xDim)
% TT       - linear transform applied to latent variables (xDim x xDim)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu
    
xDim = size(L, 2);

if xDim == 1
    TT    = sqrt(L' * L);
    Lorth = L / TT;
else
    [UU, DD, VV] = svd(L, 'econ');
    % TT is transform matrix
    TT = DD * VV';
    Lorth = UU;
end

% Apply transformation to latents
if ~isempty(X)
    Xorth = TT * X;
else
    Xorth = [];
end