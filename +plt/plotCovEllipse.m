function p = plotCovEllipse(M,C,col,varargin)
% Plot covariance ellipse.
%
% [p] = plotCovEllipse(M,C,col)
%
% Inputs:
%   M       Ellipse mean
%   C       Covariance matrix
%   col     Color to use
%
% Adapted from ellipse.m in DataHigh:
% (c) Benjamin Cowley, Matthew Kaufman, Zachary Butler, Byron Yu, 2012-2013
%
%proj_vecs: 2xN, where N=num dims
%  y = u * x, so cov(y) = u*cov(x)*u'
%cov_y = proj_vecs * cov_matrix * proj_vecs';
%
% Authors:      Alan Degenhart and Erinn Grigsby
% Emails:       erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional Inputs
LineStyle = '-';    % Line Style for covariance ellipse
Marker = 'o';       % Marker for the mean
MarkerSize = 10;    % Marker size
chisquare_val = []; % Chi^2 scale factor
dim = [1,2];        % Dimension to plto along
plotMean = true;    % Plot mean
LineWidth = 1;      % Ellipse's line width

assignopts(who,varargin);

% find the eigenvectors (should be only two)
[u,lam] = pcacov(C);

% find ellipse that'd go on that plane (the two principal comps)

if size(dim,2)==3
    [theta,phi] = meshgrid(0:.01:2*pi,0:.01:pi);
    r = sqrt(diag(lam(dim))) * [cos(theta(:)').*sin(phi(:)');...
                                sin(theta(:)').*sin(phi(:)');...
                                cos(phi(:)')];
else
    theta = 0:.01:2*pi;
    r = sqrt(diag(lam(dim))) * [cos(theta); sin(theta)];
end

if ~isempty(chisquare_val)
    r = chisquare_val*r;
end
% project onto the 2d space 
p = u(:,dim) * r;

% add on the cluster's mean (note this isn't the total mean of all
% stim)
p = p + M * ones(1,size(p,2));

% Plot ellipse
hold on
if plotMean
    if size(dim,2)==3
        k = convhull(p(dim(1),:)',p(dim(2),:)');
        plot3(p(dim(1),k),p(dim(2),k),p(dim(3),k),'color',col,'LineStyle',LineStyle, 'LineWidth', LineWidth)
        plot3(M(dim(1)),M(dim(2)),M(dim(3)),'MarkerSize',MarkerSize,'MarkerFaceColor',col, ...
            'MarkerEdgeColor',col,'Marker',Marker)
    else
        plot(p(dim(1),:),p(dim(2),:),'color',col,'LineStyle',LineStyle, 'LineWidth', LineWidth)
        plot(M(dim(1)),M(dim(2)),'MarkerSize',MarkerSize,'MarkerFaceColor',col, ...
            'MarkerEdgeColor',col,'Marker',Marker)
    end
end