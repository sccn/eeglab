% function checking if a specific file is a script or a function

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

function bool = isscript(fileName)

fid = fopen(fileName, 'r');
cont = true;
while cont && ~feof(fid)
    l = strtok(fgetl(fid));
    
    if ~isempty(l) && l(1) ~= '%'
        if strcmpi(l, 'function')
            bool = false; return;
        else
            bool = true; return;
        end
    end
end
bool = true; return;
