function printPDF(F,filename)
% printPDF          Save figure to PDF as displayed.
%
% This function saves the figure specified by figure handle F as a PDF
% file.  As opposed to calling the 'saveas' function, which uses a
% letter-size output, this function re-sizes the paper size to match the
% way the figure is currently displayed.
%
% Inputs:
%   F           Figure handle
%   filename    Name (including path if desired)
%
% Author:       Alan D. Degenhart
% Date created: N/A
% Last updated: N/A
% Last update:  N/A

set(gcf,'renderer','painters') % Needed to save editable figure

pos = get(F,'Position');
dpi = get(0,'ScreenPixelsPerInch');

paperSize = pos(3:4)/dpi; % Convert to paper size (in inches)

set(F,'PaperSize',paperSize)
set(F,'PaperPosition',[0 0 paperSize])

dpiStr = sprintf('-r%d',dpi);
print(F,'-dpdf',filename,dpiStr)