% cell2mat() - convert cell array to matrix
%
% Usage: >> M = cell2mat( C );
%
% Author: Arnaud Delorme, CNL / Salk Institute, Jan 25 2002
%
% Note: this function overload the nnet toolbox function cell2mat, 
% but does not have all its capacities. You can delete the current  
% function if you have the toolbox.

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

function M = cell2mat( C, varargin );

if nargin < 1
	help cell2mat;
	return;
end;

for i=1:size(C,1)
    for j=1:size(C,2)
        M(i,j) = C{i,j};
    end;
end;   