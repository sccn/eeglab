% h() - history function.           
%
% Usage:
%   >> h( arg );
%   >> h( arg1, arg2 );
%
% Inputs:
%   - With no argument, it return the command history.
%   - arg is a string:   with a string argument it pulls the command 
%                        onto the stack.
%   - arg is a number>0: execute the element in the stack at the 
%                        required position.
%   - arg is a number<0: unstack the required number of elements
%   - arg is 0         : clear stack
%   - arg1 is 'find' and arg2 is a string, try to find the closest command
%     in the stack containing the string
%
% Global variables used:
%   LASTCOM   - last command
%   ALLCOM    - all the commands   
%
% See also:
%  eeglab() (a graphical interface for eeg plotting, space frequency
%  decomposition, ICA, ... under Matlab for which this command
%  was designed).

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% To increase/decrease the maximum depth of the stack, edit the eeg_consts file
 
% $Log: not supported by cvs2svn $
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function str = h( command, str );

mode = 1; % mode = 1, full print, mode = 0, truncated print

try
	eeg_global;
	eeg_consts;
catch
	global LASTCOM;
	global ALLCOM;
end;

if nargin < 1
	if isempty(ALLCOM)
		fprintf('No history\n');
	else	
      for index = 1:length(ALLCOM)
         if mode == 0, txt = ALLCOM{ index }; fprintf('%d: ', index);
         else          txt = ALLCOM{ length(ALLCOM)-index+1 };
         end;   
         if (length(txt) > 72) & (mode == 0)
				fprintf('%s...\n', txt(1:70) );
			else
				fprintf('%s\n', txt );
			end;				
		end;
	end;	
elseif nargin == 1
	if isempty( command )
		return;
	end;
	if isstr( command )
		if isempty(ALLCOM)
			ALLCOM = { command };
		else	
			ALLCOM          = { command ALLCOM{:}};
		end;	
		LASTCOM         = command;
	else	
		if command == 0
			ALLCOM = [];
		else if command < 0
				ALLCOM = ALLCOM( -command+1:end ); % unstack elements
				h;
			else
				txt = ALLCOM{command};
				if length(txt) > 72
					fprintf('%s...\n', txt(1:70) );
				else
					fprintf('%s\n', txt );
				end;				
				evalin( 'base', ALLCOM{command} ); % execute element
				h( ALLCOM{command} );    % add to history
			end;
		end;	
	end;		
else % nargin == 2
	if strcmp(command, 'find')
      for index = 1:length(ALLCOM)
		  if ~isempty(findstr(ALLCOM{index}, str))
			  str = ALLCOM{index};  
			  return;
		  end;
	  end;
	  str = [];
	end;
end;
