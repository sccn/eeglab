% eegh() - history function.           
%
% Usage:
%   >> eegh( arg );
%   >> eegh( arg1, arg2 );
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
%   - arg1 is a string and arg2 is a structure, also add the history to
%     the structure in filed 'history'.
%
% Global variables used:
%   LASTCOM   - last command
%   ALLCOM    - all the commands   
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2001
%
% See also:
%  eeglab() (a graphical interface for eeg plotting, space frequency
%  decomposition, ICA, ... under Matlab for which this command
%  was designed).

% Copyright (C) 2001 Arnaud Delorme, SCCN/INC/UCSD, arno@salk.edu
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
 
function str = eegh( command, str );

mode = 1; % mode = 1, full print, mode = 0, truncated print

global ALLCOM;

%if nargin == 2
%    fprintf('2: %s\n', command);
%elseif nargin == 1
%    fprintf('1: %s\n', command);
%end;

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
    if nargout > 0
        str = strvcat(ALLCOM);
    end;
elseif nargin == 1
	if isempty( command )
		return;
	end;
	if isstr( command )
		if isempty(ALLCOM)
			ALLCOM = { command };
		else	
			ALLCOM = { command ALLCOM{:}};
		end;	
        global LASTCOM;
		LASTCOM  = command;
	else	
		if command == 0
			ALLCOM = [];
		else if command < 0
				ALLCOM = ALLCOM( -command+1:end ); % unstack elements
			else
				txt = ALLCOM{command};
				if length(txt) > 72
					fprintf('%s...\n', txt(1:70) );
				else
					fprintf('%s\n', txt );
				end;				
				evalin( 'base', ALLCOM{command} ); % execute element
				eegh( ALLCOM{command} );    % add to history
			end;
		end;	
	end;		
else % nargin == 2
    if ~isstruct(str)
        if strcmp(command, 'find')
            for index = 1:length(ALLCOM)
                if ~isempty(findstr(ALLCOM{index}, str))
                    str = ALLCOM{index};  
                    return;
                end;
            end;
            str = [];
        end;
    else
        % warning also some code present in eeg_store and pop_newset
        eegh(command); % add to history
        if ~isempty(command)
            if length(str) == 1
                str = eeg_hist(str, command);
            else
                for i = 1:length(str)
                    str(i) = eeg_hist(str(i), [ '% multiple datasets command: ' command ]);
                end;
            end;
        end;
    end;
end;
