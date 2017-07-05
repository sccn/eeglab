% recursively get the list of Matlab file in 
% a given directory
%
% A. Delorme, 2017

function fileList = getmfilelist(folder)

fileList = {};
if nargin < 1
    folder = pwd;
end

allFiles = dir(folder);

for iFile = 1:length(allFiles)
    
    if ~isequal(allFiles(iFile).name, '..') && ~isequal(allFiles(iFile).name, '.')
        
        if exist(fullfile(folder, allFiles(iFile).name), 'dir')
            fileList2 = getmfilelist(fullfile(folder, allFiles(iFile).name));
            fileList = [ fileList fileList2 ];
        elseif isequal(allFiles(iFile).name(end-1:end), '.m')
            fileList{end+1} = fullfile(folder, allFiles(iFile).name);
        end
        
    end
end
        
        