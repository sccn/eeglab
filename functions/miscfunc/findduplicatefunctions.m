% script to find potential conflict between Matlab functions in the current
% path
%
% A. Delorme, May 25th, 2012

% Copyright (C) A. Delorme, May 25th, 2012
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

folders = path;
delim = find(folders == ':');
delim = [ 0 delim length(folders)+1 ];

for iFolder = 1:length(delim)-1
    currentFolder = folders(delim(iFolder)+1:delim(iFolder+1)-1);
    
    folderContent = dir(currentFolder);
    
    if ~contains( currentFolder, matlabroot)
        
        rmpath(currentFolder);
        folderContent = { folderContent.name };
        
        for iFile = 1:length(folderContent)
            
            currentFile = folderContent{iFile};
            if length(currentFile) > 2 && strcmpi(currentFile(end-1:end), '.m')
                
                if ~strcmpi(currentFile, 'Contents.m')
                    tmp = which(currentFile);

                    if ~isempty(tmp)
                        fullFileName1 = fullfile(currentFolder, currentFile);
                        fullFileName2 = which(currentFile);
                        %                         cond1 = contains(fullFileName1, 'toolbox/matlab') || contains(fullFileName1, 'toolbox\matlab');
                        cond1 = contains(fullFileName1, matlabroot);
                        if cond1
                            sadffdas
                        end
                        %                         if ~cond1
                        fprintf('Potential conflict between %s and %s\n', fullFileName1, fullFileName2);
                    end
                end
            end
        end
        addpath(currentFolder);
        
    end
end
