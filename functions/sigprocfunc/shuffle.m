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

function [x, i, j]=shuffle( y, dim)

if nargin < 1
	help shuffle;
	return;
end
if nargin < 2
	if size(y,1) ~= 1
		dim = 1;
	else
		if size(y,2) ~= 1
			dim = 2;
		else
			dim = 3;
		end
	end
end
	
r =size(y, dim);
a = rand(1,r);
[tmp i] = sort(a);
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
