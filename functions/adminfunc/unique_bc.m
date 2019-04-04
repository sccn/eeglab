% unique_bc - unique backward compatible with Matlab versions prior to 2013a

% Copyright (C) 2013 Arnaud Delorme
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

function [C,IA,IB] = unique_bc(A,varargin)

errorFlag = error_bc;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));
if v > 7.19, v = floor(v) + rem(v,1)/10; end

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end
end

if v >= 7.14
    [C,IA,IB] = unique(A,varargin{:},'legacy');
    if errorFlag
        [C2,IA2] = unique(A,varargin{:});
        if ~isequal(C, C2) || ~isequal(IA, IA2) || ~isequal(IB, IB2)
            warning('backward compatibility issue with call to unique function');
        end
    end
else
    [C,IA,IB] = unique(A,varargin{:});
end