% sum() - sum of memory mapped underlying array
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

function sumval = sum(obj,dim)
    
    if nargin < 2
        dim = 1;
    end;
    
    s1  = size(obj);
    ss.type = '()';
    ss.subs(1:length(s1)) = { ':' };
    for index = 1:s1(dim)
        ss.subs{dim} = index;
        if index == 1
             sumval = subsref(obj, ss);
        else sumval = sumval + subsref(obj, ss);
        end;
    end;
