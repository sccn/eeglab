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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
		if indexA > 0 && indexA <= length(A)
            R(index) = R(index) + B(indexB)*A(indexA);
            sumconvo = sumconvo + B(indexB);
		end
	end
	R(index) = R(index) / sumconvo;
end

return;
