% pop_delset() - Delete a dataset from the variable containing
%                all datasets.
%
% Usage: >> ALLEEG = pop_delset(ALLEEG, indices);
%
% Inputs:
%   ALLEEG   - array of EEG datasets
%   indices  - indices of datasets to delete. None -> a pop_up window asks 
%              the user to choose. Index < 0 -> it's positive is given as 
%              the default in the pop-up window (ex: -3 -> default 3).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_copyset(), eeglab()

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

% load a set and store it in the current set
% ------------------------------------------
function [ALLSET, command] = pop_delset(ALLSET, set_in);

command = '';
if nargin < 1
	help pop_delset;
	return;
end;
if isempty( ALLSET )
	error('Cannot delete dataset. Restart eeglab to clear all dataset information');
    return;
end;    

if nargin < 2 | set_in < 0
	% which set to delete
	% -----------------
	promptstr    = { 'Dataset(s) to delete:' };
	if nargin == 2
		inistr       = { int2str(-set_in) };
	else
		inistr       = { '1' };
	end;
	result       = inputdlg2( promptstr, 'Delete dataset -- pop_delset()', 1,  inistr, 'pop_delset');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	set_in   	 = eval( [ '[' result{1} ']' ] );
end;

if isempty(set_in)
	return;
end;	

A = fieldnames( ALLSET );
A(:,2) = cell(size(A));
A = A';
for i = set_in
    try
   		ALLSET(i) = struct(A{:});
		%ALLSET = setfield(ALLSET, {set_in}, A{:}, cell(size(A)));
	catch
		error('Error: no such dataset');
		return;
	end;
end;
    
command = sprintf('%s = pop_delset( %s, [%s] );', inputname(1), inputname(1), int2str(set_in));
return;
