% recursively get the list of Matlab file in 
% a given directory
%
% A. Delorme, 2017

% Copyright (C) A. Delorme, 2017
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
        
        
