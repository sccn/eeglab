% reshape() - reshape of memory mapped underlying array
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

function a = reshape(a,d1,d2,d3)
    
    % decode length
    % -------------
    if nargin > 3
        d1 = [ d1 d2 d3 ];
    elseif nargin > 2
        d1 = [ d1 d2 ];
    end;
    
    if prod(size(a)) ~= prod(d1)
        error('Wrong dimensions for reshaping');
    end;
    
    if ~strcmpi(a.fileformat, 'transposed')
        a.data.format{2} = d1;
    else
        if length(d1) == 1
            a.data.format{2} = d1;
        elseif length(d1) == 2
            a.data.format{2} = [d1(2) d1(1)];
        else
            a.data.format{2} = d1([2 3 1]);
        end;
    end;
