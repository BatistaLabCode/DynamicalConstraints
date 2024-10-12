% tube.get_tube_col  Get colors for tube plotting.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [tube_col] = get_tube_col(uni_start_pos, n_tube,varargin)
C = [];         % Use the default ColMat structure to identify the scale
colStart = [];  % Set the desired brightest color

% Parse optional agruments
assignopts(who,varargin);

if isempty(C) & isempty(colStart)
    C = util.defineTaskColormap('bc');
    scale = norm(uni_start_pos(1,1:2))./100;
    C.targPos = C.targPos.*(scale/.9);
end

% Setup colors - blocks (not currently used)
if isempty(colStart)
    colInfo = util.getColorInfo(C,uni_start_pos,[]);
    colStart = colInfo{2};
end
cHSV = rgb2hsv(colStart);
colGray = ones(1,3)*.75;

% Setup colors - different tubes.  Generate a map of nTube + 1 colors, so
% that the "brightest" color is reserved for the unconstrained trials
cTubeHSV = repmat(cHSV(1,:), n_tube + 1, 1);
cTubeHSV(:,3) = linspace(cTubeHSV(1,3), cTubeHSV(1,3) * 0.15, n_tube + 1);
cTube = hsv2rgb(cTubeHSV);

% Add to output structure
tube_col.tube = cTube;
tube_col.gray = colGray;