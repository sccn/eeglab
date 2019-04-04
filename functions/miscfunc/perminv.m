% perminv() - returns the inverse permutation vector
%
% Usage: >> [invvec] = perminverse(vector);
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 11-30-96 

% Copyright (C) 11-30-96 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 4-4-97 shortened name to perminv() -sm
% 4-7-97 allowed row vector, added tests -sm
% 01-25-02 reformated help & license -ad 

function [invvec]=perminv(vector)

[n,c] = size(vector);
if n>1 && c>1,
    fprintf('perminv(): input must be a vector.\n');
	return
end
transpose=0;
if c>1
	vector = vector';
	transpose =1;
end

invvec = zeros(size(vector));
for i=1:length(vector)
  invvec(vector(i)) = i;
end

if transpose==1,
	invvec = invvec';
end
