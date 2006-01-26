% eeglab_exec() - apply eeglab function to a dataset and perform appropriate
%                 checks.
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

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
        
% $$$     catch, 
% $$$         
% $$$         % handling errords
% $$$         % ----------------
% $$$         tmp = lasterror;
% $$$         tmperr = [ 'EEGLAB error in function ' tmp.stack(1).name '() at line ' int2str(tmp.stack(1).line) ':' 10 10 lasterr ];
% $$$         tmperr = [ tmperr 10 10 'If you think it is a bug, send a detailed' 10 ...
% $$$                                 'description of how to reproduce the problem' 10 ...
% $$$                                 'and a (small) test dataset to eeglab@sccn.ucsd.edu' ];
% $$$         errordlg2(tmperr, 'EEGLAB error');
% $$$         
% $$$     end;
    
        
% checking strings
% ----------------
e_try = 'try,';
e_catch = 'catch, errordlg2(lasterr, ''EEGLAB error''); LASTCOM= ''''; clear EEGTMP ALLEEGTMP STUDYTMP; end;';
nocheck           = e_try;
ret = 'if ~isempty(LASTCOM), if LASTCOM(1) == -1, LASTCOM = ''''; return; end; end;';
check             = ['[EEG LASTCOM] = eeg_checkset(EEG, ''data'');' ret ' eegh(LASTCOM);' e_try];
checkcont         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''contdata'');' ret ' eegh(LASTCOM);' e_try];
checkica          = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'');' ret ' eegh(LASTCOM);' e_try];
checkepoch        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'');' ret ' eegh(LASTCOM);' e_try];
checkevent        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''event'');' ret ' eegh(LASTCOM);' e_try];
checkbesa         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''besa'');' ret ' eegh(''% no history yet for BESA dipole localization'');' e_try];
checkepochica     = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'');' ret ' eegh(LASTCOM);' e_try];
checkplot         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''chanloc'');' ret ' eegh(LASTCOM);' e_try];
checkicaplot      = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'', ''chanloc'');' ret ' eegh(LASTCOM);' e_try];
checkepochplot    = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''chanloc'');' ret ' eegh(LASTCOM);' e_try];
checkepochicaplot = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'', ''chanloc'');' ret ' eegh(LASTCOM);' e_try];

% check string and backup old dataset
% -----------------------------------
backup =     [ 'if CURRENTSET ~= 0,' ...
               '    [ ALLEEG EEG ] = eeg_store(ALLEEG, EEG, CURRENTSET, ''savegui'');' ...
               '    eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET, ''''savedata'''');'');' ...
               'end;' ];

nocheck_back           = [ backup e_try ];
check_back             = ['[EEG LASTCOM] = eeg_checkset(EEG, ''data'');' ret 'eegh(LASTCOM);' backup e_try];
checkcont_back         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''contdata'');' ret ' eegh(LASTCOM);' backup e_try];
checkica_back          = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'');' ret ' eegh(LASTCOM);' backup e_try];
checkepoch_back        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'');' ret ' eegh(LASTCOM);' backup e_try];
checkevent_back        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''event'');' ret ' eegh(LASTCOM);' backup e_try];
checkbesa_back         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''besa'');' ret ' eegh(''% no history yet for BESA dipole localization'');' backup e_try];
checkepochica_back     = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'');' ret ' eegh(LASTCOM);' backup e_try];
checkplot_back         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''chanloc'');' ret ' eegh(LASTCOM);' backup e_try];
checkicaplot_back      = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'', ''chanloc'');' ret ' eegh(LASTCOM);' backup e_try];
checkepochplot_back    = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''chanloc'');' ret ' eegh(LASTCOM);' backup e_try];
checkepochicaplot_back = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'', ''chanloc'');' ret ' eegh(LASTCOM);' backup e_try];

storecall    = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);'');';
storeload    = '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, 0); eegh(''[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);'');';
storenewcall = '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); eegh(LASTCOM);';
storeallcall = [ 'if ~isempty(ALLEEG) & ~isempty(ALLEEG(1).data), ALLEEG = eeg_checkset(ALLEEG);' ...
                 'EEG = eeg_checkset(EEG); eegh(''ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_checkset(EEG);''); end;' ];

ifeegtmp     =  'if ~isempty(LASTCOM) & ~isempty(EEGTMP),';
ifeeg        =  'if ~isempty(LASTCOM) & ~isempty(EEG),';
ifeegnh      =  'if ~isempty(LASTCOM) & ~isempty(EEG) & ~isempty(findstr(''='',LASTCOM)),';

e_newnonempty_nh     = [e_catch 'eegh(LASTCOM);' ifeegtmp 'EEG = EEGTMP;' storenewcall 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_load_nh            = [e_catch 'eegh(LASTCOM);' ifeegtmp 'EEG = EEGTMP;' storeload    'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_newset_nh          = [e_catch 'eegh(LASTCOM);' ifeeg                    storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store_nh           = [e_catch 'eegh(LASTCOM);' ifeegnh                  storecall    'disp(''Done.''); end; eeglab(''redraw'');'];
e_storeall_nh        = [e_catch 'eegh(LASTCOM);' ifeeg                    storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist_nh            = [e_catch 'eegh(LASTCOM);'];
e_histdone_nh        = [e_catch 'eegh(LASTCOM); if ~isempty(LASTCOM), disp(''Done.''); end;' ];

% same as above but also save history in dataset
% ----------------------------------------------
e_newnonempty   = [e_catch 'EEG = eegh(LASTCOM, EEGTMP);' ifeegtmp 'EEG = EEGTMP;' storenewcall 'disp(''Done.''); end; clear EEGTMP; eeglab(''redraw'');'];
e_load          = [e_catch 'EEG = eegh(LASTCOM, EEGTMP);' ifeegtmp 'EEG = EEGTMP;' storeload    'disp(''Done.''); end; clear EEGTMP; eeglab(''redraw'');'];
e_newset        = [e_catch 'EEG = eegh(LASTCOM, EEG);' ifeeg                    storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store         = [e_catch 'EEG = eegh(LASTCOM, EEG);' ifeegnh                  storecall    'disp(''Done.''); end; eeglab(''redraw'');'];
e_storeall      = [e_catch 'EEG = eegh(LASTCOM, EEG);' ifeeg                    storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist          = [e_catch 'EEG = eegh(LASTCOM, EEG);'];
e_histdone      = [e_catch 'EEG = eegh(LASTCOM, EEG); if ~isempty(LASTCOM), disp(''Done.''); end;' ];

% study checking
% --------------
e_load_study    = [e_catch 'eegh(LASTCOM); if ~isempty(LASTCOM), ALLEEG = ALLEEGTMP; EEG = ALLEEG; CURRENTSET = [1:length(EEG)]; eegh(''CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:' int2str(length(EEG)) '];''); STUDY = STUDYTMP; CURRENTSTUDY = 1; disp(''Done.''); end; clear ALLEEGTMP STUDYTMP; eeglab(''redraw'');'];

% build structures for plugins
% ----------------------------
trystrs.no_check                 = e_try;
trystrs.no_check_back            = [ backup e_try ];
trystrs.check_data               = check;
trystrs.check_ica                = checkica;
trystrs.check_cont               = checkcont;
trystrs.check_epoch              = checkepoch;
trystrs.check_event              = checkevent;
trystrs.check_epoch_ica          = checkepochica;
trystrs.check_chanlocs           = checkplot;
trystrs.check_epoch_chanlocs     = checkepochplot;
trystrs.check_epoch_ica_chanlocs = checkepochicaplot;
catchstrs.add_to_hist            = e_hist;
catchstrs.store_and_hist         = e_store;
catchstrs.new_and_hist           = e_newset;
catchstrs.new_non_empty          = e_newnonempty;
