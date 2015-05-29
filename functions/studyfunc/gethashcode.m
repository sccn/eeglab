% gethashcode() - get hash code for input array
%
% Usage:
%         >> md5code = gethashcode(array);
%
% Inputs:
%       array - array of number. gets converted to double.
%
% Output:
%  md5code   - MD5 output code
%
% Author: Arnaud Delorme, 2015-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function res= gethashcode(str)

import java.security.*;
import java.math.*;
md = MessageDigest.getInstance('MD5');
hash = md.digest(double(str));
bi = BigInteger(1, hash);
res = char(bi.toString(16));
