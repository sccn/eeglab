% size() - size of memory mapped underlying array
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Nov. 2008

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD
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

function [s varargout] = size(obj,dim)
    
    if obj.transposed
        if length(obj.dimensions) ~= 2 && length(obj.dimensions) ~= 3
            error('Transposed array must be 2 or 3 dims');
        end;
        if length(obj.dimensions) == 2 tmpdimensions = [obj.dimensions(2) obj.dimensions(1)];
        else                           tmpdimensions = [obj.dimensions(3) obj.dimensions(1:end-1)];
        end;
    else
        tmpdimensions = obj.dimensions;
    end;
    
    s = tmpdimensions;
    
    if nargin > 1
        if dim >length(s)
            s = 1;
        else
            s = s(dim);
        end;
    else
        if nargout > 1
            s = [s 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
            alls = s;
            s = s(1);
            for index = 1:max(nargout,1)-1
                varargout{index} = alls(index+1);
            end;
        end;
    end;
    
    