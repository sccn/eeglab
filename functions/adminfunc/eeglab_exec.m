% eeglab_exec() - apply eeglab function to a dataset and perform appropriate
%                 checks. LEGACY FUNCTION. NOT USED ANYWHERE IN EEGLAB.
%
% Usage:
%   >> OUTEEG = eeglab_exec(funcname, INEEG, 'key1', value1, 'key2', value2 ...);
%
% Inputs:
%   funcname      - [string] name of the function
%   INEEG         - EEGLAB input dataset(s)
%
% Optional inputs
%   'params'      - [cell array] funcname parameters.
%   'warning'     - ['on'|'off'] warning pop-up window if several dataset
%                   stored on disk that will be automatically overwritten.
%                   Default is 'on'.
%
% Outputs:
%   OUTEEG        - output dataset(s)
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-
%
% see also: eeglab()

% Copyright (C) 2006 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function eeglab_exec( func, varargin);

    eeg_global;
    % global variables are EEG ALLEEG CURRENTSET ALLCOM LASTCOM STUDY CURRENTSTUDY

    com = '';
    if nargin < 2
        help eeglab_exec;
        return;
    end;
    
    %try
        % check input parameters
        % ----------------------
        g = finputcheck( varargin, { ...
                            'input'    { 'string' 'cell' }  { {} {} }        {};
                            'output'   { 'string' 'cell' }  { {} {} }        {};
                            'check'    { 'string' 'cell' }  { {} {} }        {};
                            'backup'   'string'  { 'on' 'off' }              'off';
                            'load'     'string'  { 'eeg' 'study' 'none' }    'none';
                            'datasethist' 'string'  { 'on' 'off' }           'off';
                            'store'    'string'  { 'neweeg' 'sameeeg' 'loadeeg' 'samestudy' 'loadstudy' 'none' }     'none' }, ...
                            'eeglab_exec');
        if isstr(g), error(g); end;
        if ~iscell(g.input),  g.input  = { g.input  }; end;
        if ~iscell(g.output), g.output = { g.output }; end;
        if ~iscell(g.check),  g.check  = { g.check  }; end;
        
        % check dataset
        % -------------
        LASTCOM = '';
        if ~isempty(g.check)
            [EEG LASTCOM] = eeg_checkset(EEG, g.check{:});
        end;
        
        % save history from check
        % -----------------------
        if ~isempty(LASTCOM),
            if LASTCOM(1) == -1, LASTCOM = ''; return; end;
        end;
        eegh(LASTCOM);
        
        % basic command to run
        % --------------------
        if isstr(func)
            eval( [ 'func = @' func ';' ] );
        end;

        % backup data on disk if necessary
        % --------------------------------
        if strcmpi(g.backup, 'on')
            if CURRENTSET ~= 0,
               [ ALLEEG EEG ] = eeg_store(ALLEEG, EEG, CURRENTSET, 'savegui');
               eegh('[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET, ''savedata'');');
            end;
        end;

        % run command
        % -----------
        switch length(g.output)
         case 0, feval(func, g.input{:}); % Matlab 6 and higher
         case 1, out{1} = feval(func, g.input{:}); % Matlab 6 and higher
         case 2, [out{1} out{2}] = feval(func, g.input{:}); % Matlab 6 and higher
         case 3, [out{1} out{2} out{3}]  = feval(func, g.input{:}); % Matlab 6 and higher
         case 4, [out{1} out{2} out{3} out{4}] = feval(func, g.input{:}); % Matlab 6 and higher
        end;
        
        % assign outputs
        % --------------
        for index = length(g.output):-1:1 % so LASTCOM gets assigned first
            
            % test if empty output (for instance when loading data)
            % --------------------
            if strcmpi(g.output{index}, 'EEG') | strcmpi(g.output{index}, 'ALLEEG') | strcmpi(g.output{index}, 'STUDY')
                if isempty(out{index}), return;
                end;
            end;
            
            eval( [ g.output{index} '= out{index};' ]); % this is where global variables are modified
        end;
        
        % store dataset history (LASTCOM can be empty)
        % ---------------------
        if strcmpi(g.datasethist, 'on')
            EEG = eegh(LASTCOM, EEG);
        else
            eegh(LASTCOM);
        end;
        
        % test for LASTCOM
        % ----------------
        if ~isempty(LASTCOM)
            
            % loading EEG data
            % ----------------
            switch g.store
             case 'loadeeg', [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, 0); 
                       eegh('[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);');
             case 'neweeg',  [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); eegh(LASTCOM);
             case 'sameeeg', [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 
                       eegh('[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);');
             case 'samestudy',;
             case 'loadstudy', CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)]; 
                       eegh([ 'CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:' int2str(length(EEG)) '];' ]); 
             case 'none',;
            end;
            
            evalin('base', 'eeglab redraw');
        end;
        
        disp('Done.');
        
