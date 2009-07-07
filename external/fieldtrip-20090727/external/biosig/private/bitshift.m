function [B] = bitshift(A,k,N)
% BITSHIFT performs a shift of the binary bitpattern 
%	of the integer number represented by A.  
%	k > 0 performs a left shift, k<0 performs a right shift. 
%	N ensures that B < 2^N
%    
%   B = bitshift(A, k [,N])
%	returns B = mod(A*2^k,2^N) 
%

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%       $Revision: 1.1 $
%       $Id: bitshift.m,v 1.1 2009-07-07 02:23:48 arno Exp $
%	Copyright (C) 2004 by Alois Schloegl <a.schloegl@ieee.org>	
%       This function is part of the NaN-toolbox
%       http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/NaN/


A = fix(A);
k = fix(k);

if nargin < 3,
	% allocate output memory and check size of argument
	B = fix(A.*(2.^k));	% size of input arguments do not fit - if an error occurs in this line  

	%B = fix(2.^(log2(A)+k));	% size of input arguments do not fit - if an error occurs in this line  
else
	% allocate output memory and check size of argument
	B = rem(fix(A.*(2.^k)), 2^N);	% size of input arguments do not fit - if an error occurs in this line  

	%B = fix(2.^(mod(log2(A)+k,N-)));	% size of input arguments do not fit - if an error occurs in this line  
end;
