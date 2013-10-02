% find potential conflict between Matlab functions
%
% A. Delorme, May 25th, 2012

matlabRoot = fileparts(fileparts(fileparts(fileparts(which('ispc')))));
folders = path;
delim = find(folders == ':');
delim = [ 0 delim length(folders)+1 ];

for iFolder = 1:length(delim)-1
    currentFolder = folders(delim(iFolder)+1:delim(iFolder+1)-1);
    
    folderContent = dir(currentFolder);
    
    if isempty(findstr(matlabRoot, currentFolder))
        
        rmpath(currentFolder);
        folderContent = { folderContent.name };
        
        for iFile = 1:length(folderContent)
            
            currentFile = folderContent{iFile};
            if length(currentFile) > 2 && strcmpi(currentFile(end-1:end), '.m')
                
                if ~strcmpi(currentFile, 'Contents.m')
                    if ~isempty(which(currentFile))
                        fullFileName1 = fullfile(currentFolder, currentFile);
                        fullFileName2 = which(currentFile);
                        fprintf('Potential conflict between %s and %s\n', fullFileName1, fullFileName2);
                    end;
                end;
            end;
        end;
        addpath(currentFolder);
        
    end;
end;