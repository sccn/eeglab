% pophelp() - Same as hthelp() but does not crash under windows.
%
% Usage: >> pophelp( function );
%
% Inputs:
%   function  - string for a Matlab function name 
%               (with or without the '.m' extension).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab() 

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
% 01-25-02 reformatted help & license -ad 
% 3/19/02 Edited help message (below) -sm 

function pophelp( funct, header );

if nargin <2
	header = 1;
end;
	
if header == 1
	doc   = { }; 
else
	doc = {};
end;			 

if findstr( funct, '.m')
	fid = fopen( funct, 'r');
else
	fid = fopen( [funct '.m'], 'r');
end;
if fid == -1
	error('File not found');
end;

str = fgets( fid );
while (str(1) == '%')
	str = deblank(str(1:end-1));

	doc = { doc{:} str(2:end) };
	str = fgets( fid );
end;

textgui(doc);
return;
