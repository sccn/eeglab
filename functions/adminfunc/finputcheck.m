% finputcheck() - check matlab function { 'key', 'val' } inputs
%
% Usage: >> struct = finputcheck( varargin, fieldlist, callingfunc );
%
% Input:
%   varargin  - varargin arguement from a function call with 'key', 'val'
%               arguements.
%   fieldlist - 3 or 4 columns cell array with one row per variable. The first
%               column contain the variable name, the second one the type, 
%                the third the accepted value range and the fourth one the 
%               defaultvalue. Ex:
%                  { 'varanme1' { 'value1' 'value2' } 'defaultvaluevar1' }
%                  { 'varanme2' { int1 int2 } 'defaultvaluevar2' }
%                  etc...
%               allowed types are 'boolean', 'integer', 'real', 'string', 
%               'cell' or 'struct'
%               the fifth column may contain the size (can be a cell array 
%               of size), but is optional.
%  callingfunc - calling function name for error messages. Default is none.
%
% Outputs:
%   struct    - checked structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu
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

function g = finputcheck( vararg, fieldlist, callfunc )

	if nargin < 2
		help finputcheck;
		return;
	end;
	if nargin < 3
		callfunc = '';
	end;

	NAME = 1;
	TYPE = 2;
	VALS = 3
	DEF  = 4;
	SIZE = 5;
	
	% create structure
	% ----------------
	if ~isempty(vararg)
		for index=1:length(vararg)
			if iscell(vararg{index})
				vararg{index} = {vararg{index}};
			end;
		end;
		try
			g = struct(vararg{:});
		catch
			error([ callfunc 'error: bad ''key'', ''val'' sequence' ]); 
		end;
	end;
	
	for index = 1:size(fieldlist,NAME)
		% check if present
		% ----------------
		if ~isfield(g, fieldlist{index, NAME})
			g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
		end;
		tmpval = getfield( g, {1}, fieldlist{index, NAME});
		
		% check type
		% ----------
		switch fieldlist{index, TYPE}
		 case { 'integer' 'real' 'boolean' }, 
		  if ~isnumeric(tmpval)
			error([ callfunc 'error: argument ''' fieldlist{index, NAME} ''' must be numeric' ]); 
		  end;
		  if strcmp(fieldlist{index, TYPE}, 'boolean')
			  if tmpval ~=0 & tmpval ~= 1
				  error([ callfunc 'error: argument ''' fieldlist{index, NAME} ''' must be 0 or 1' ]); 
			  end;  
		  else 
			  if strcmp(fieldlist{index, TYPE}, 'integer')
				  if ~isempty(fieldlist{index, VALS})
					  if ~ismember(tmpval, fieldlist{index, VALS})
						  error([ callfunc 'error: wrong value for argument ''' fieldlist{index, NAME} '''' ]); 
					  end;
				  end;
			  else % real
				  if ~isempty(fieldlist{index, VALS})
					  if tmpval < fieldlist{index, VALS}(1) | tmpval > fieldlist{index, VALS}(2)
						  error([ callfunc 'error: value out of range for argument ''' fieldlist{index, NAME} '''' ]); 
					  end;
				  end;
			  end;
		  end;  
			
		  
		 case 'string'
		  if ~isstr(tmpval)
			error([ callfunc 'error: argument ''' fieldlist{index, NAME} ''' must be a string' ]); 
		  end;
		  if ~isempty(fieldlist{index, VALS})
			  if isempty(strmatch(lower(tmpval)), fieldlist{index, VALS})
				  error([ callfunc 'error: wrong value for argument''' fieldlist{index, NAME} '''' ]); 
			  end;
		  end;

		  
		 case 'cell'
		  if ~isstr(tmpval)
			error([ callfunc 'error: argument ''' fieldlist{index, NAME} ''' must be a cell array' ]); 
		  end;
		  
		  
		 case 'struct'
		  if ~isstr(tmpval)
			error([ callfunc 'error: argument ''' fieldlist{index, NAME} ''' must be a structure' ]); 
		  end;
		  
		  
		 case '';
		 otherwise, error([ 'finputcheck error: unrecognized type ''' fieldlist{index, NAME} '''' ]);
		end;
	end;
	
	% check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	allowedfields =  fieldlist(:, 1)';
	for index=1:length(allfields)
		switch allfields{index}
		 case allowedfields;
		 otherwise, error([ callfunc 'error: undefined argument ''' fieldlist{index, NAME} '''']);
		end;
	end;
