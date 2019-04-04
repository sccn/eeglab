% compare EEGLAB distributions

% Copyright (C) 2006 Arnaud Delorme
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

function compareeeglabdistrib(folder1, folder2)

fileList1 = getmfilelist(folder1);
fileList2 = getmfilelist(folder2);

% strip list from folders
for iFile = 1:length(fileList1), fileList1{iFile} = fileList1{iFile}(length(folder1)+1:end); end
for iFile = 1:length(fileList2), fileList2{iFile} = fileList2{iFile}(length(folder2)+1:end); end

for iFile1 = 1:length(fileList1)
    iFile2 = strmatch(fileList1{iFile1}, fileList2, 'exact');
    if isempty(iFile2)
        fprintf('Missing file %s in folder %s\n', fileList1{iFile1}, folder2);
    else
        dir1 = dir(fullfile(folder1, fileList1{iFile1}));
        dir2 = dir(fullfile(folder2, fileList2{iFile2}));
        if length(dir1) ~= 1 || length(dir2) ~= 1 
            error('Size error');
        end
        if dir1.bytes ~= dir2.bytes
            fprintf('Files %s differ\n', fileList2{iFile2});
            %fprintf('---------------------\n');
            %system([ '/usr/bin/diff ' fullfile(folder1, fileList1{iFile1}) ' ' fullfile(folder2, fileList2{iFile2}) ]);
        end
        fileList2(iFile2) = [];
    end
end

for iFile2 = 1:length(fileList2)
    fprintf('Missing file %s in folder %s\n', fileList2{iFile2}, folder1);
end
    


            
            
        
    
