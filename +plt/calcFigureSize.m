function [fW,fH,Ax] = calcFigureSize(nRow,nCol,axW,axH,axSp,varargin)
% Creates a figure to the exact specification of the number and size of the
% subplots.
%
% Inputs:
%   nRow        Number of rows
%   nCol        Number of columns
%   axW         Axis width (pixels)
%   axH         Axis height (pixels)
%   sp          Axis spacing (x and y) (pixels)
%
% Optional Inputs:
%   margin      Margins (px) 
%
% Outputs:
%   fW          Figure width
%   fH          Figure Height
%   Ax          subplot axes parameters
%
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional arguments
margin = 75; % Margin (px)

assignopts(who, varargin);

% Calculate figure size
fW = axW*nCol + axSp*(nCol-1) + margin*2;
fH = axH*nRow + axSp*(nRow-1) + margin*2;

% Calculate fractional values (for use with the 'subplotSimple' command
Ax.xMarg = margin/fW;
Ax.xSp = axSp/fW;
Ax.yMarg = margin/fH;
Ax.ySp = axSp/fH;
Ax.nRow = nRow;
Ax.nCol = nCol;