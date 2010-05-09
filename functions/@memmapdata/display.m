% display() - display an EEG data class underlying structure
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

function b = display(a)
 
    i.type = '()';
    i.subs = { ':' ':' ':' };
    b = subsref(a, i); % note that subsref cannot be called directly
    return;
    
  
    %struct(a)
    %return;

    if ~strcmpi(a.fileformat, 'transposed')
        a.data.data.x;
    else
        permute(a, [3 1 2]);
    end;
