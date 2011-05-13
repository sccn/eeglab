% eeg_eval() - apply eeglab function to a collection of input datasets
%
% Usage:
%   >> OUTEEG = eeg_eval(funcname, INEEG, 'key1', value1, 'key2', value2 ...);
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
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005
% 
% see also: eeglab()

% Copyright (C) 2005 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = eeg_eval( funcname, EEG, varargin);
    
    com = '';
    if nargin < 2
        help eeg_eval;
        return;
    end;
    
    % check input parameters
    % ----------------------
    g = finputcheck( varargin, { 'params'   'cell'    {}               {};
                                 'warning'  'string'  { 'on','off' }   'on' }, 'eeg_eval');
    if isstr(g), error(g); end;
    
    % warning pop up
    % --------------
    eeglab_options;
    if strcmpi(g.warning, 'on')
        if ~option_storedisk
            res = questdlg2(strvcat( 'When processing multiple datasets, it is not', ...
                                 'possible to enter new names for the newly created', ...
                                 'datasets and old datasets are overwritten.', ...
                                 'You may still cancel this operation though.'), ...
                 'Multiple dataset warning', 'Cancel', 'Proceed', 'Proceed');
        else
            res = questdlg2( [ 'Data files on disk will be automatically overwritten.' 10 ...
                                'Are you sure you want to proceed with this operation?' ], ...
                            'Confirmation', 'Cancel', 'Proceed', 'Proceed');
        end;
        switch lower(res),
         case 'cancel', return;
         case 'proceed',;
        end;
    end;
 
    % execute function
    % ----------------
    v = version;
    if str2num(v(1)) == '5' % Matlab 5
        command = [ 'TMPEEG = ' funcname '( TMPEEG, ' vararg2str(g.params) ');' ];
    else
        eval( [ 'func = @' funcname ';' ] );
    end;
        
    NEWEEG = [];
    for i = 1:length(EEG)
        fprintf('Processing group dataset %d of %d named: %s ****************\n', i, length(EEG), EEG(i).setname);
        TMPEEG    = eeg_retrieve(EEG, i);
        if v(1) == '5', eval(command);                      % Matlab 5
        else            TMPEEG = feval(func, TMPEEG, g.params{:}); % Matlab 6 and higher
        end;
        TMPEEG = eeg_checkset(TMPEEG);
        TMPEEG.saved = 'no';
        if option_storedisk
            TMPEEG = pop_saveset(TMPEEG, 'savemode', 'resave');
            TMPEEG = update_datafield(TMPEEG);
        end;
        NEWEEG = eeg_store(NEWEEG, TMPEEG, i);
        if option_storedisk
            NEWEEG(i).saved = 'yes'; % eeg_store by default set it to no
        end;
    end;
    EEG = NEWEEG;
    
    % history
    % -------
    if nargout > 1
        com = sprintf('%s = %s( %s,%s);', funcname, inputname(2), funcname, inputname(2), vararg2str(g.params));
    end;

function EEG = update_datafield(EEG);
    if ~isfield(EEG, 'datfile'), EEG.datfile = ''; end;
    if ~isempty(EEG.datfile)
        EEG.data = EEG.datfile;
    else 
        EEG.data = 'in set file';
    end;
    EEG.icaact = [];
