% perminv() - returns the inverse permutation vector
%
% Usage: >> [invvec] = perminverse(vector);
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 11-30-96 

% Copyright (C) 11-30-96 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 4-4-97 shortened name to perminv() -sm
% 4-7-97 allowed row vector, added tests -sm
% 01-25-02 reformated help & license -ad 

function [invvec]=perminv(vector)

[n,c] = size(vector);
if n>1 & c>1,
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
end;

if transpose==1,
	invvec = invvec';
end
