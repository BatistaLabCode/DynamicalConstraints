function TU = calcBoundary(TU,varargin)
% Calculate boundary method for Tube class
% p - tube path
% r - tube radius
% Adapted from getTubeBoundary by Patrick Sadtler
%
% Usage:
%   TU = TU.calcBoundary
%
% Optional Inputs:
%   FLATENDS    Determine whether to plot tubes with rounded(default) or 
%                   flat end caps.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Return if tube path or radius is empty
if isempty(TU.path) || sum(sum(~isnan(TU.window))) == 0
    return
end

p = TU.path(:,1:2); % Currently only valid for 2D

% define default values for the optional arguments
FLATENDS = false;               % default for getting flat or round ends on the tubes

% load the optional arguments
assignopts('ignorecase',who,varargin);

% Set radius from window.  Currently assume that the first non-NAN value in
% the 'window' is used for the entire trial.
w = TU.window(:,1);
w = w(~isnan(w));
if length(unique(w)) > 1
    warning('More than one unique tube tolerance window found.')
end
r = w(1);

% Initialize upper and lower boundaries
nPts = size(p,1);
bU = nan(nPts,2);
bD = nan(nPts,2);

% Check radius.  This should either be a scalar or vector the length of the
% path.
if length(r) == 1
    r = r*ones(nPts,1);
elseif length(r) ~= nPts
    error('Invalid tube radius dimensionality.')
end
if size(r,1) == 1; r = r'; end

% Find the boundaries associated with each point in the tube
for i = 1:nPts
    % Find the tangent unit vector at each point
    if i == 1
        tV = p(2,:) - p(1,:);
    elseif i == nPts
        tV = p(end,:) - p(end-1,:);
    else
        tV = p(i+1,:) - p(i-1,:);
    end
    tV = tV/norm(tV);
    
    % Find the normal vector and corresponding boundary points
    nV = [-tV(2) tV(1)];
    bU(i,:) = p(i,:) + r(i)*nV;
    bD(i,:) = p(i,:) - r(i)*nV;
end

% Handle tube end.  Usually, this is curved, but for some reason a 'flat'
% tube was also an option
if ~FLATENDS
    nAng = 59; % Number of points to use for defining the tube end
    bEnd1 = nan(nAng,2);
    bEnd2 = nan(nAng,2);
    
    % Get the normal vectors
    tU = p(2,:) - p(1,:);
    tU = tU/norm(tU);
    nU = [-tU(2) tU(1)];
    tD = p(end,:) - p(end-1,:);
    tD = tD/norm(tD);
    nD = [-tD(2) tD(1)];

    % Define the rotation angles
    angU = linspace(0,180,nAng+2);
    angU = angU(2:end-1);
    angD = -angU;

    % Rotate normal vectors and find the boundary point
    for i = 1:nAng
        rU = [cosd(angU(i)) -sind(angU(i));sind(angU(i)) cosd(angU(i))];
        rU = rU*nU';
        bEnd1(i,:) = p(1,:) + r(1)*rU';

        rD = [cosd(angD(i)) -sind(angD(i));
            sind(angD(i)) cosd(angD(i))];
        rD = rD*nD';
        bEnd2(i,:) = p(end,:) + r(end)*rD';
    end
else
    bEnd1 = [];
    bEnd2 = [];
end

% Concatenate points into a continuous sequence
b = [bU; bEnd2; flipud(bD); flipud(bEnd1)];
tubeRadiiConcat = [r; r(end)*ones(size(bEnd2,1),1); 
    flipud(r); flipud(r(1)*ones(size(bEnd1,1),1))];

% Remove overlapping points
nPts = size(b,1);
minDist = nan(nPts,1);
for i = 1:size(b,1)
    minDist(i) = round(min(sqrt(sum(bsxfun(@minus,b(i,:),p).^2,2)))*10^4)/10^4;
    % Note: not sure why the scaling is going on here -- proabably to deal
    % with rounding errors
end
mask = round(minDist*10^2)/10^2 < round(tubeRadiiConcat*10^2)/10^2;
b = b(~mask,:);
b(end+1,:) = b(1,:); % Add first point to make a closed region
TU.boundary = b;
