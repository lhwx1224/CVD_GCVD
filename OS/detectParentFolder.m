function parentFolder = detectParentFolder()
    % Get the full path of the current MATLAB file
    currentFilePath = mfilename('fullpath');

    % Extract the parent folder from the current file path
    [FilePath, ~, ~] = fileparts(currentFilePath);
    
    % Finding parent folder 
    folders = split(FilePath,filesep);
    parentFolder = folders{end-1};
    % Display the parent folder
    disp(['Parent folder: ' parentFolder]);
end