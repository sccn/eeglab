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

% load a set and store it in the current set
% ------------------------------------------
function [ALLSET, command] = pop_delset(ALLSET, set_in);

command = '';
if nargin < 1
	help pop_delset;
	return;
end
if isempty( ALLSET )
	error('Cannot delete dataset. Restart eeglab to clear all dataset information');
    return;
end;    

if nargin < 2 || set_in < 0
	% which set to delete
	% -----------------
	promptstr    = { 'Dataset(s) to delete:' };
	if nargin == 2
		inistr       = { int2str(-set_in) };
	else
		inistr       = { '1' };
	end
	result       = inputdlg2( promptstr, 'Delete dataset -- pop_delset()', 1,  inistr, 'pop_delset');
	size_result  = size( result );
	if size_result(1) == 0 return; end
	set_in   	 = eval( [ '[' result{1} ']' ] );
end

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
	end
end
    
% command = sprintf('%s = pop_delset( %s, [%s] );', inputname(1), inputname(1), int2str(set_in));
command = sprintf('EEG = pop_delset( EEG, [%s] );', int2str(set_in));
return;
