% script to find potential conflict between Matlab functions in the current
% path
%
% A. Delorme, May 25th, 2012

% Copyright (C) A. Delorme, May 25th, 2012
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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
                    end
                end
            end
        end
        addpath(currentFolder);
        
    end
end
