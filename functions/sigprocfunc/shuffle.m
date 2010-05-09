% shuffle() - shuffle a given dimension in an array
%
% Usage: >> Y = shuffle(X)
%        >> [Y = shuffle(X, DIM)
% 
% Inputs: 
%   X   - input array
%   DIM - dimension index (default is firt non-singleton dimention)
%
% Outputs: 
%    Y - shuffled array
%    I - forward indices (Y = X(I) if 1D)
%    J - reverse indices (X(J) = Y if 1D)
%
% Author: Arnaud Delorme, SCCN/INC/UCSD USA, Dec 2000

% Copyright (C) Arnaud Delorme, SCCN/INC/UCSD USA, Dec 2000
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

function [x, i, j]=shuffle( y, dim)

if nargin < 1
	help shuffle;
	return;
end;
if nargin < 2
	if size(y,1) ~= 1
		dim = 1;
	else
		if size(y,2) ~= 1
			dim = 2;
		else
			dim = 3;
		end;
	end;
end;
	
r =size(y, dim);
i=1:r;
for j=1:r,
  p=fix(1+r*rand);
  a=i(p);
  i(p)=i(j);
  i(j)=a;
end;
switch dim
	case 1
		x = y(i,:,:,:,:);
	case 2
		x = y(:,i,:,:,:);
	case 3
		x = y(:,:,i,:,:);
	case 4
		x = y(:,:,:,i,:);
	case 5
		x = y(:,:,:,:,i);
end;		
[tmp j] = sort(i); % unshuffle

return;
