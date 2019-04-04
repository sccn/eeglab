% remove all path with a given parent path
% varargin contains a list of path to exclude

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

function removepath(parentpath, varargin)

if isempty(parentpath), return; end
folder = path;
if ispc, sep = ';'; else sep = ':'; end
indSep = find(folder == sep);

indSep = [ 0 indSep length(folder)+1 ];
for iSep = 1:length(indSep)-1
    curPath = folder(indSep(iSep)+1:indSep(iSep+1)-1);
    if ~isempty(strfind(curPath, parentpath)) && ~any(strmatch(curPath, varargin, 'exact'))
        rmpath(curPath);
    end
end
