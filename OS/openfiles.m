function [file, path, selectedfile] = openfiles(filetype, description)
% OPENFILE() opens a file through graphical use interface
% syntax: [file, path, selectedfile] = openfile(filetype)
% 
% input: filetype (string type) - extension name of the file you would like
% to load, e.g., '*.mat' selects the matrix file, '*.txt' selects the text
% file.
% output: file (characters type) - file name that you selected using the
%         GUI 
%         path (characters type) - file path that you selected using the
%         GUI
%         selectedfile (character type) - full file path to which the file
%         you selected from the GUI
%
% 
% Hewenxuan Li, Feb 2021
% Last modified: 
% October 1, 2021, assertion is added to replace the if condition

if nargin<2
    description = [];
end

extnname = split(filetype,'.'); % get the extension name
extnname = extnname{end};       % ----------------------

[file, path] = uigetfile(filetype,['Select the file', description],[],'MultiSelect','on');

% if isequal(file, 0)
%     disp('No file has selected! User canceled the selection!')
% else
%     selectedfile = fullfile(path,file);
%     disp(['---- loading:',file,' ----'])
%     disp(['Selected File:', selectedfile])
% end

% IF NO FILE SELECTED, FLAG == 1
NoFileFlag = isequal(file, 0);
assert(NoFileFlag == 0,'No File Selected! openfile.m aborted!')

selectedfile = fullfile(path,file);
if ~ischar(file)
    disp('===== Mulitple Files Selected =====')
    % Check if all the selected files are data frame files
    for i = 1:size(selectedfile,2)
        disp([file{i}, ' will be loaded!'])
    end
else
    disp('===== Single Files Selected =====')
    disp([file, ' will be loaded!'])
    tempfile = file;
    file = cell(1,1);
    file{1,1} = tempfile;
end

file = file(~cellfun('isempty',file));

% BEGIN FOR selected data files
selectedfile = fullfile(path,file);

end