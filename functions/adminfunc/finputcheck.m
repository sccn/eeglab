% finputcheck() - check matlab function { 'key', 'val' } inputs
%
% Usage: >> struct = finputcheck( varargin, fieldlist );
%        >> [struct varargin] = finputcheck( varargin, fieldlist, ... 
%                                                      callingfunc, mode );
%
% Input:
%   varargin  - varargin arguement from a function call with 'key', 'val'
%               arguements.
%   fieldlist - 3 or 4 columns cell array with one row per variable. The first
%               column contain the variable name, the second one the type, 
%                the third the accepted value range and the fourth one the 
%               defaultvalue. For instance
%                  { 'varanme1' { 'value1' 'value2' } 'defaultvaluevar1' }
%                  { 'varanme2' { int1 int2 } 'defaultvaluevar2' }
%                  etc...
%               allowed types are 'boolean', 'integer', 'real', 'string', 
%               'cell' or 'struct'
%               the fifth column may contain the size (can be a cell array 
%               of size), but is optional.
%  callingfunc - calling function name for error messages. Default is none.
%  mode        - ['ignore'|'error'] ignore keywords that are not specified in
%                the fieldlist cell array or generate an error. Default is
%                'error'.
% Outputs:
%   struct     - checked structure
%   varargin   - residual varagin with non-recognized arguments
%
% Note: in case of error, a string is returned with the error message
%       instead of a structure.
%
% Example:
%	finputcheck(varargin, ...
%               { 'title'         'string'   []                       ''; ...
%                 'percent'       'real'     [0 1]                     1 ; ...
%                 'elecamp'       'integer'  [1:10]                   [] });
%   'title' is a string with no default value
%   'percent' is a real number in between 0 and 1 and default value 1
%   'elecamp' is an integer that can tak value between 1 and 10
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
% Revision 1.11  2002/09/30 15:29:23  arno
% autorizing cell arrays for types
%
% Revision 1.10  2002/09/30 00:42:08  arno
% debug input arguments
%
% Revision 1.9  2002/07/29 18:00:53  arno
% debugging for NaN
%
% Revision 1.8  2002/07/29 17:24:22  arno
% header
%
% Revision 1.7  2002/07/20 19:10:41  arno
% debugging output
%
% Revision 1.6  2002/07/19 17:58:11  arno
% returning non-matched 'key' 'val' arguments
%
% Revision 1.5  2002/07/19 17:46:53  arno
% g empty if no varargin
%
% Revision 1.4  2002/07/19 16:27:14  arno
% adding ignore mode
%
% Revision 1.3  2002/07/10 02:18:32  arno
% header info
%
% Revision 1.2  2002/07/10 02:17:27  arno
% debugging error message passing
%
% Revision 1.1  2002/07/10 01:03:19  arno
% Initial revision
%

function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode )

	if nargin < 2
		help finputcheck;
		return;
	end;
	if nargin < 3
		callfunc = '';
	else 
		callfunc = [callfunc ' ' ];
	end;
    if nargin < 4
        mode = 'do not ignore';
    end;
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
	
	varargnew = {};
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
			g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
		end;
	else 
		g = [];
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
        if ~iscell( fieldlist{index, TYPE} )
            res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
            if isstr(res), g = res; return; end;
        else 
            testres = 0;
            tmplist = fieldlist;
            for it = 1:length( fieldlist{index, TYPE} )
                res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
                if ~isstr(res{it}), testres = 1; end;
            end;
            if testres == 0, g = strvcat(res{:}); return; end;
        end;
	end;
    
    % check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	for index=1:length(allfields)
		if isempty(strmatch(allfields{index}, fieldlist(:, 1)'))
			if ~strcmpi(mode, 'ignore')
				g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
			end;
			varargnew{end+1} = allfields{index};
			varargnew{end+1} = getfield(g, {1}, allfields{index});
		end;
	end;


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
    g = [];
    
    switch fieldtype
     case { 'integer' 'real' 'boolean' }, 
      if ~isnumeric(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
      end;
      if strcmp(fieldtype, 'boolean')
          if tmpval ~=0 & tmpval ~= 1
              g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
          end;  
      else 
          if strcmp(fieldtype, 'integer')
              if ~isempty(fieldval)
                  if (isnan(tmpval) & ~any(isnan(fieldval))) ...
                          & (~ismember(tmpval, fieldval))
                      g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                  end;
              end;
          else % real
              if ~isempty(fieldval)
                  if tmpval < fieldval(1) | tmpval > fieldval(2)
                      g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                  end;
              end;
          end;
      end;  
      
      
     case 'string'
      if ~isstr(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
      end;
      if ~isempty(fieldval)
          if isempty(strmatch(lower(tmpval), fieldval))
              g = [ callfunc 'error: wrong value for argument''' fieldname '''' ]; return;
          end;
      end;

      
     case 'cell'
      if ~iscell(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
      end;
      
      
     case 'struct'
      if ~isstruct(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
      end;
      
      
     case '';
     otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
    end;
