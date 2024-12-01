function saveFigurePDF(F,saveDir,varargin)
% saveFigurePDF(F,saveDir)
%
% Save multiple figures to PDF.
%
% This function saves the figures with handles provided by F to the
% directory specified by SAVEDIR. This function uses the 'name' property of
% the figure as the saved PDF file name, so this property must be set
% before this function is called.
%
% Inputs:
%   F               Figure handle
%   saveDir         Save location path name
%
% Optional Inputs:
%   appendFiles     Combine all figures into a single PDF
%   appendName      Combined PDF's file name
%   closeFigs       Close figures after successful saving
%
% Author:       Alan D. Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Optional arguments
appendFiles = false;
appendName = [];
closeFigs = false;

assignopts(who,varargin);

% If save directory is not specified, prompt user
if (nargin == 1) || isempty(saveDir)
    saveDir = uigetdir;
end

nFig = length(F);
saveFilenames = cell(nFig,1);
for i = 1:nFig
    fName = get(F(i),'Name');
    tempFilename = fullfile(saveDir,[fName '.pdf']);
    printPDF(F(i),tempFilename)
    saveFilenames{i} = tempFilename;
    if closeFigs
        close(F(i))
    end
end

% Save appended figure if desired (this may take some time)
if appendFiles
    outputFilename = fullfile(saveDir,[appendName '.pdf']);
    append_pdfs(outputFilename,saveFilenames{:})
end