function checkdir(dirName)
% Check if a directory exists base on the given directory name 'dirName'
% check if directory exists
if ~exist(dirName, 'dir')
    % if directory does not exist, create it
    mkdir(dirName);
    disp(['Directory ' dirName ' created.']);
else
    disp(['Directory ' dirName ' already exists.']);
end