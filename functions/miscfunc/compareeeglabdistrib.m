% compare EEGLAB distributions

% Copyright (C) 2006 Arnaud Delorme
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
    


            
            
        
    
