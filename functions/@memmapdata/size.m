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

function [s s2 s3] = size(a,dim)
    
    if isnumeric(a.data), 
        s = size(a.data)
    else
        s = a.data.format{2};
        if strcmpi(a.fileformat, 'transposed')
            if length(s) == 2, s = s([2 1]);
            elseif length(s) == 3
                s = [s(3) s(1) s(2)];
            end;
        end;
    end;
    
    if nargin > 1
        s = [s 1];
        s = s(dim);
    end;
    
    if nargout > 2
        s3 = s(3);
    end;
    if nargout > 1
        s2 = s(2);
        s  = s(1);
    end;
