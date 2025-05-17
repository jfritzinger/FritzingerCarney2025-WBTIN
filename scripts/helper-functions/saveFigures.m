function saveFigures(filename)
%SAVEFIGURES Save current figure to standardized location and format
%   SAVEFIGURES(FILENAME) saves the current active figure to the WBTIN
%   project's designated save path as a high-resolution TIFF file.
%
%   Input:
%       filename - Base filename (without extension) as string/char
%
%   Outputs:
%       None (saves file to disk)

[~, ~, savepath, ~] = getPathsWBTIN();
% print(fullfile(savepath, [filename '.png']),...
% 	'-dpng', '-r600')
print('-vector', fullfile(savepath, [filename '.tif']),...
	'-dtiff', '-r600')
% print('-vector', fullfile(savepath, [filename '.eps']),...
% 	'-depsc', '-r600')


end