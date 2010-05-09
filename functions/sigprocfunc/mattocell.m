% mattocell() - convert matrix to cell array
%
% Usage: >> C = mattocell( M );
%
% Author: Arnaud Delorme, CNL / Salk Institute, Jan 25 2002
%
% Note: this function overload the nnet toolbox function mattocell, 
% but does not have all its capacities. You can delete the current  
% function if you have the toolbox.

% Copyright (C) Jan 25 2002 Arnaud Delorme, CNL / Salk Institute  
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

function C = mattocell( M, varargin );

if nargin < 1
	help mattocell;
	return;
end;
if isempty(M)
	C = [];
	return;
end;
C = cell(size(M));
for i=1:size(M,1)
    for j=1:size(M,2)
        C{i,j} = M(i,j);
    end;
end;   
