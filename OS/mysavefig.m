function mysavefig(Path, plotName)
% MYSAVEFIG.m saves the selected figure using a path and a plotName and
% save the .png, .pdf, and the .fig files of the plot.
% input - Path, a string of the path to save the plot
%         plotName, a string without extension name
% output - 'Path/plotName.png' file
%          'Path/plotName.pdf' file
%          'Path/plotName.fig' file
% Hewenxuan Li, Dec, 2023 @Cornell
pngName = strjoin({plotName,'.png'},'');
pdfName = strjoin({plotName,'.pdf'},'');
figName = strjoin({plotName,'.fig'},'');
print(strjoin({Path, pngName},filesep),'-dpng', '-r600')
print(strjoin({Path, pdfName},filesep),'-dpdf')
savefig(strjoin({Path, figName},filesep))
end