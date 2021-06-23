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

function [EEG, com] = eeg_eval( funcname, EEG, varargin)
    
    com = '';
    if nargin < 2
        help eeg_eval;
        return;
    end
    
    % check input parameters
    % ----------------------
    g = finputcheck( varargin, { 'params'   'cell'    {}               {};
                                 'warning'  'string'  { 'on','off' }   'off' }, 'eeg_eval');
    if ischar(g), error(g); end
    
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
        end
        switch lower(res)
            case 'cancel', return;
            case 'proceed'
        end
        vPar = ver('parallel');
        if option_parallel && ~isempty(vPar)
            res = questdlg2(strvcat( 'You have selected to use the parallel toolbox,', ...
                                 'to process multiple datasets. If you saturate the ', ...
                                 'memory, this could cause your computer to become.', ...
                                 'unresponsive or even crash.'), ...
                 'Multiple dataset warning', 'Cancel', 'Proceed', 'Proceed');
        end
        switch lower(res)
            case 'cancel', return;
            case 'proceed'
        end
    end
 
    % execute function
    % ----------------
    v = version;
    if str2num(v(1)) == '5' % Matlab 5
        command = [ 'TMPEEG = ' funcname '( TMPEEG, ' vararg2str(g.params) ');' ];
    elseif isstr(funcname)
        eval( [ 'func = @' funcname ';' ] );
    else
        func = funcname;
    end
    
%     notCompatibleFunc = { 
%         @clean_artifacts
%     };
    
    NEWEEG = EEG;
    if option_parallel
        disp('Using the parallel toolbox to process multiple datasets (change in File > Preferences)');
        myPool = gcp;
        tmpoptions_store = option_storedisk;
        parfor i = 1:length(EEG)
            fprintf('Processing group dataset %d of %d named: %s ****************\n', i, length(EEG), EEG(i).setname);
            TMPEEG    = eeg_retrieve(EEG, i);
            TMPEEG = feval(func, TMPEEG, g.params{:});
            TMPEEG = eeg_checkset(TMPEEG);
            TMPEEG.saved = 'no';
            if tmpoptions_store
                TMPEEG = pop_saveset(TMPEEG, 'savemode', 'resave');
                TMPEEG = update_datafield(TMPEEG);
                NEWEEG(i) = TMPEEG;
                NEWEEG(i).saved = 'yes'; % eeg_store by default set it to no
            else
                NEWEEG(i) = TMPEEG;
            end
        end
    else
        for i = 1:length(EEG)
            fprintf('Processing group dataset %d of %d named: %s ****************\n', i, length(EEG), EEG(i).setname);
            TMPEEG    = eeg_retrieve(EEG, i);
            TMPEEG = feval(func, TMPEEG, g.params{:});
            TMPEEG = eeg_checkset(TMPEEG);
            TMPEEG.saved = 'no';
            if option_storedisk
                TMPEEG = pop_saveset(TMPEEG, 'savemode', 'resave');
                TMPEEG = update_datafield(TMPEEG);
            end
            NEWEEG = eeg_store(NEWEEG, TMPEEG, i);
            if option_storedisk
                NEWEEG(i).saved = 'yes'; % eeg_store by default set it to no
            end
        end
    end
    EEG = NEWEEG;
    
    % history
    % -------
    if nargout > 1
        com = sprintf('%s = %s( %s,%s);', inputname(2), char(funcname), inputname(2), vararg2str(g.params));
    end

function EEG = update_datafield(EEG)
    if ~isfield(EEG, 'datfile'), EEG.datfile = ''; end
    if ~isempty(EEG.datfile)
        EEG.data = EEG.datfile;
    else 
        EEG.data = 'in set file';
    end
    EEG.icaact = [];
