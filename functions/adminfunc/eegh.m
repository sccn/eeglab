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

% To increase/decrease the maximum depth of the stack, edit the eeg_consts file
 
function str = eegh( command, str );

mode = 1; % mode = 1, full print, mode = 0, truncated print

global ALLCOM;

%if nargin == 2
%    fprintf('2: %s\n', command);
%elseif nargin == 1
%    fprintf('1: %s\n', command);
%end

if nargin < 1
	if isempty(ALLCOM)
		fprintf('No history\n');
	else	
      for index = 1:length(ALLCOM)
         if mode == 0, txt = ALLCOM{ index }; fprintf('%d: ', index);
         else          txt = ALLCOM{ length(ALLCOM)-index+1 };
         end;   
         if (length(txt) > 72) && (mode == 0)
				fprintf('%s...\n', txt(1:70) );
			else
				fprintf('%s\n', txt );
			end;				
		end
	end;	
    if nargout > 0
        str = strvcat(ALLCOM);
    end
elseif nargin == 1
	if isempty( command )
		return;
	end
	if ischar( command )
        if ~isempty(ALLCOM) && isequal(ALLCOM{1}, command), return; end
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
			end
		end;	
	end;		
else % nargin == 2
    if ~isstruct(str)
        if strcmp(command, 'find')
            for index = 1:length(ALLCOM)
                if ~isempty(findstr(ALLCOM{index}, str))
                    str = ALLCOM{index};  
                    return;
                end
            end
            str = [];
        end
    else
        % warning also some code present in eeg_store and pop_newset
        if ~isempty(ALLCOM) && isequal(ALLCOM{1}, command), return; end
        eegh(command); % add to history
        if ~isempty(command)
            if length(str) == 1
                str = eeg_hist(str, command);
            else
                for i = 1:length(str)
                    str(i) = eeg_hist(str(i), [ '% multiple datasets command: ' command ]);
                end
            end
        end
    end
end


