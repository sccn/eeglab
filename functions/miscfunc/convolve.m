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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% 01-25-02 reformated help & license -ad 

function R = convolve( A, B )

% size of B must be odd

% convolve A and B (normalize by the sum of convolved elements)
% -------------------------------------------------------------
SIZECONVO = floor(size(B(:),1)/2);
R = zeros( size(A) );
for index = 1:size(A(:),1)
	sumconvo = 0;
	for indexconvo = -SIZECONVO:SIZECONVO
		indexA = index + indexconvo;
		indexB = indexconvo+SIZECONVO+1;
		if indexA > 0 & indexA < size(A(:),1)
			if indexB > 0 & indexB < size(B(:),1)
				R(index) = R(index) + B(indexB)*A(indexA);
				sumconvo = sumconvo + B(indexconvo+SIZECONVO+1);
			end;	
		end;
	end;
	R(index) = R(index) / sumconvo;
end;

return;
