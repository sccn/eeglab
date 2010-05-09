% convolve() - convolve two matrices (normalize by the sum of convolved 
%              elements to compensate for border effects).
%
% Usage:
%   >> r = convolve( a, b );
%
% Inputs:
%   a          - first input vector
%   b          - second input vector
%
% Outputs:
%   r          - result of the convolution
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: conv2(), conv()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function R = convolve( A, B )

% size of B must be odd

% convolve A and B (normalize by the sum of convolved elements)
% -------------------------------------------------------------
R = zeros( size(A) );
for index = 1:length(A)
	sumconvo = 0;
	for indexB = 1:length(B)
        indexconvo = indexB-1-round(length(B)/2);
		indexA = index + indexconvo+1;
		if indexA > 0 & indexA <= length(A)
            R(index) = R(index) + B(indexB)*A(indexA);
            sumconvo = sumconvo + B(indexB);
		end;
	end;
	R(index) = R(index) / sumconvo;
end;

return;
