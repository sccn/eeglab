% subsasgn() - define index assignment for eegdata objects
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

function b = subsasgn(a,index,val)
    
    i.type = '()';
    i.subs = { ':' ':' ':' };
    b = subsref(a, i); % note that subsref cannot be called directly
    c = subsref(val, i);
    b = builtin('subsasgn', b, index, c);
    return;
    
    
switch index.type
 case '()'
  switch length(index.subs)
   case 1, a.data(index.subs{1}) = val;
   case 2, a.data(index.subs{1}, index.subs{2}) = val;
   case 3, a.data(index.subs{1}, index.subs{2}, index.subs{3}) = val;
   case 4, a.data(index.subs{1}, index.subs{2}, index.subs{3}, index.subs{4}) = val;
  end;
 case '.'
  switch index.subs
   case 'srate'
    a.srate = val;
   case 'nbchan'
    a.nbchan = val;
   otherwise
    error('Invalid field name')
   end
end
