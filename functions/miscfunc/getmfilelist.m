% recursively get the list of Matlab file in 
% a given directory
%
% A. Delorme, 2017

% Copyright (C) A. Delorme, 2017
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
        
        
