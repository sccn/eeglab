% eeglab() - Matlab graphic user interface environment for 
%   electrophysiological data analysis incorporating the ICA/EEG toolbox 
%   (Makeig et al.) developed at CNL / The Salk Institute, 1997-2001. 
%   Released 11/2002- as EEGLAB (Delorme, Makeig, et al.) at the Swartz Center 
%   for Computational Neuroscience, Institute for Neural Computation, 
%   University of California San Diego (http://sccn.ucsd.edu/). 
%   User feedback welcome: email eeglab@sccn.ucsd.edu
%
% Authors: Arnaud Delorme and Scott Makeig, with substantial contributions
%   from Colin Humphries, Sigurd Enghoff, Tzyy-Ping Jung, plus
%   contributions 
%   from Tony Bell, Te-Won Lee, Luca Finelli and many other contributors. 
%
% Description:
%   EEGLAB is Matlab-based software for processing continuous or event-related 
%   EEG or other physiological data. It is designed for use by both novice and 
%   expert Matlab users. In normal use, the EEGLAB graphic interface calls 
%   graphic functions via pop-up function windows. The EEGLAB history mechanism 
%   can save the resulting Matlab calls to disk for later incorporation into 
%   Matlab scripts.  A single data structure ('EEG') containing all dataset 
%   parameters may be accessed and modified directly from the Matlab commandline. 
%   EEGLAB now recognizes "plugins," sets of EEGLAB functions linked to the EEGLAB
%   main menu through an "eegplugin_[name].m" function (Ex. >> help eeplugin_besa.m). 
%
% Usage: 1) To (re)start EEGLAB, type
%            >> eeglab           % Ignores any loaded datasets
%        2) To redaw and update the EEGLAB interface, type
%            >> eeglab redraw    % Scans for non-empty datasets
%            >> eeglab rebuild   % Closes and rebuilds the EEGLAB window
%            >> eeglab versions  % State EEGLAB version number
%
%   >> type "license.txt" % the GNU public license
%   >> web http://sccn.ucsd.edu/eeglab/tutorial/ % the EEGLAB tutorial
%   >> help eeg_checkset  % the EEG dataset structure
%
% GUI Functions calling eponymous processing and plotting functions:
% ------------------------------------------------------------------
% <a href="matlab:helpwin pop_eegfilt">pop_eegfilt</a>   - bandpass filter data (eegfilt())
% <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>   - scrolling multichannel data viewer (eegplot())
% <a href="matlab:helpwin pop_eegthresh">pop_eegthresh</a> - simple thresholding method (eegthresh())
% <a href="matlab:helpwin pop_envtopo">pop_envtopo</a>   - plot ERP data and component contributions (envtopo())
% <a href="matlab:helpwin pop_epoch">pop_epoch</a>     - extract epochs from a continuous dataset (epoch())
% <a href="matlab:helpwin pop_erpimage">pop_erpimage</a>  - plot single epochs as an image (erpimage())
% <a href="matlab:helpwin pop_jointprob">pop_jointprob</a> - reject epochs using joint probability (jointprob())
% <a href="matlab:helpwin pop_loaddat">pop_loaddat</a>   - load Neuroscan .DAT info file (loaddat())
% <a href="matlab:helpwin pop_loadcnt">pop_loadcnt</a>   - load Neuroscan .CNT data (lndcnt())
% <a href="matlab:helpwin pop_loadeeg">pop_loadeeg</a>   - load Neuroscan .EEG data (loadeeg())
% <a href="matlab:helpwin pop_loadbva">pop_loadbva</a>   - load Brain Vision Analyser matlab files
% <a href="matlab:helpwin pop_plotdata">pop_plotdata</a>  - plot data epochs in rectangular array (plotdata())
% <a href="matlab:helpwin pop_readegi">pop_readegi</a>   - load binary EGI data file (readegi())
% <a href="matlab:helpwin pop_rejkurt">pop_rejkurt</a>   - compute data kurtosis (rejkurt())
% <a href="matlab:helpwin pop_rejtrend">pop_rejtrend</a>  - reject EEG epochs showing linear trends  (rejtrend())
% <a href="matlab:helpwin pop_resample">pop_resample</a>  - change data sampling rate (resample())
% <a href="matlab:helpwin pop_rmbase">pop_rmbase</a>    - remove epoch baseline (rmbase())
% <a href="matlab:helpwin pop_runica">pop_runica</a>    - run infomax ICA decomposition (runica())
% <a href="matlab:helpwin pop_newtimef">pop_newtimef</a>  - event-related time-frequency (newtimef())
% <a href="matlab:helpwin pop_timtopo">pop_timtopo</a>   - plot ERP and scalp maps  (timtopo())
% <a href="matlab:helpwin pop_topoplot">pop_topoplot</a>  - plot scalp maps (topoplot())
% <a href="matlab:helpwin pop_snapread">pop_snapread</a>  - read Snapmaster .SMA files (snapread())
% <a href="matlab:helpwin pop_newcrossf">pop_newcrossf</a> - event-related cross-coherence (newcrossf())
% <a href="matlab:helpwin pop_spectopo">pop_spectopo</a>  - plot all channel spectra and scalp maps (spectopo())
% <a href="matlab:helpwin pop_plottopo">pop_plottopo</a>  - plot a data epoch in a topographic array (plottopo())
% <a href="matlab:helpwin pop_readedf">pop_readedf</a>   - read .EDF EEG data format (readedf())
% <a href="matlab:helpwin pop_headplot">pop_headplot</a>  - plot a 3-D data scalp map (headplot())
% <a href="matlab:helpwin pop_reref">pop_reref</a>     - re-reference data (reref())
% <a href="matlab:helpwin pop_signalstat">pop_signalstat</a> - plot signal or component statistic (signalstat())
%
% Other GUI functions:
% -------------------
% <a href="matlab:helpwin pop_chanevent">pop_chanevent</a>      - import events stored in data channel(s)
% <a href="matlab:helpwin pop_comments">pop_comments</a>       - edit dataset comment ('about') text
% <a href="matlab:helpwin pop_compareerps">pop_compareerps</a>    - compare two dataset ERPs using plottopo()
% <a href="matlab:helpwin pop_prop">pop_prop</a>           - plot channel or component properties (erpimage, spectra, map)
% <a href="matlab:helpwin pop_copyset">pop_copyset</a>        - copy dataset
% <a href="matlab:helpwin pop_dispcomp">pop_dispcomp</a>       - display component scalp maps with reject buttons
% <a href="matlab:helpwin pop_editeventfield">pop_editeventfield</a> - edit event fields
% <a href="matlab:helpwin pop_editeventvals">pop_editeventvals</a>  - edit event values
% <a href="matlab:helpwin pop_editset">pop_editset</a>        - edit dataset information
% <a href="matlab:helpwin pop_export">pop_export</a>         - export data or ica activity to ASCII file
% <a href="matlab:helpwin pop_expica">pop_expica</a>         - export ica weights or inverse matrix to ASCII file
% <a href="matlab:helpwin pop_icathresh">pop_icathresh</a>      - choose rejection thresholds (in development)
% <a href="matlab:helpwin pop_importepoch">pop_importepoch</a>    - import epoch info ASCII file
% <a href="matlab:helpwin pop_importevent">pop_importevent</a>    - import event info ASCII file
% <a href="matlab:helpwin pop_importpres">pop_importpres</a>     - import Presentation info file
% <a href="matlab:helpwin pop_importev2">pop_importev2</a>      - import Neuroscan ev2 file
% <a href="matlab:helpwin pop_loadset">pop_loadset</a>        - load dataset
% <a href="matlab:helpwin pop_mergeset">pop_mergeset</a>       - merge two datasets
% <a href="matlab:helpwin pop_rejepoch">pop_rejepoch</a>       - reject pre-identified epochs in a EEG dataset
% <a href="matlab:helpwin pop_rejspec">pop_rejspec</a>        - reject based on spectrum (computes spectrum -% eegthresh)
% <a href="matlab:helpwin pop_saveh">pop_saveh</a>          - save EEGLAB command history
% <a href="matlab:helpwin pop_saveset">pop_saveset</a>        - save dataset
% <a href="matlab:helpwin pop_select">pop_select</a>         - select data (epochs, time points, channels ...)
% <a href="matlab:helpwin pop_selectevent">pop_selectevent</a>    - select events
% <a href="matlab:helpwin pop_subcomp">pop_subcomp</a>        - subtract components from data
%
% Non-GUI functions use for handling the EEG structure:
% ----------------------------------------------------
% <a href="matlab:helpwin eeg_checkset">eeg_checkset</a>       - check dataset parameter consistency
% <a href="matlab:helpwin eeg_context">eeg_context</a>        - return info about events surrounding given events
% <a href="matlab:helpwin pop_delset">pop_delset</a>         - delete dataset
% <a href="matlab:helpwin pop_editoptions">pop_editoptions</a>    - edit the option file
% <a href="matlab:helpwin eeg_emptyset">eeg_emptyset</a>       - empty dataset
% <a href="matlab:helpwin eeg_epochformat">eeg_epochformat</a>    - convert epoch array to structure
% <a href="matlab:helpwin eeg_eventformat">eeg_eventformat</a>    - convert event array to structure
% <a href="matlab:helpwin eeg_getepochevent">eeg_getepochevent</a>  - return event values for a subset of event types
% <a href="matlab:helpwin eeg_global">eeg_global</a>         - global variables
% <a href="matlab:helpwin eeg_multieegplot">eeg_multieegplot</a>   - plot several rejections (using different colors)
% <a href="matlab:helpwin eeg_options">eeg_options</a>        - option file
% <a href="matlab:helpwin eeg_rejsuperpose">eeg_rejsuperpose</a>   - use by rejmenu to superpose all rejections
% <a href="matlab:helpwin eeg_rejmacro">eeg_rejmacro</a>       - used by all rejection functions
% <a href="matlab:helpwin pop_rejmenu">pop_rejmenu</a>        - rejection menu (with all rejection methods visible)
% <a href="matlab:helpwin eeg_retrieve">eeg_retrieve</a>       - retrieve dataset from ALLEEG
% <a href="matlab:helpwin eeg_store">eeg_store</a>          - store dataset into ALLEEG
%
% Help functions:
% --------------
% <a href="matlab:helpwin eeg_helpadmin">eeg_helpadmin</a>      - help on admin function
% <a href="matlab:helpwin eeg_helphelp">eeg_helphelp</a>       - help on help
% <a href="matlab:helpwin eeg_helpmenu">eeg_helpmenu</a>       - EEG help menus
% <a href="matlab:helpwin eeg_helppop">eeg_helppop</a>        - help on pop_ and eeg_ functions
% <a href="matlab:helpwin eeg_helpsigproc">eeg_helpsigproc</a>    - help on

% Copyright (C) 2001 Arnaud Delorme and Scott Makeig, Salk Institute, 
% arno@salk.edu, smakeig@ucsd.edu.
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

function varargout = eeglab( onearg )

if nargout > 0
    varargout = { [] [] 0 {} };
    %[ALLEEG, EEG, CURRENTSET, ALLCOM]
end;

% add the paths
% -------------
eeglabpath = which('eeglab.m');
eeglabpath = eeglabpath(1:end-length('eeglab.m'));
if strcmpi(eeglabpath, './') || strcmpi(eeglabpath, '.\'), eeglabpath = [ pwd filesep ]; end;

% solve BIOSIG problem
% --------------------
pathtmp = which('wilcoxon_test');
if ~isempty(pathtmp)
    try,
        rmpath(pathtmp(1:end-15));
    catch, end;
end;

% test for local SCCN copy
% ------------------------
if ~iseeglabdeployed2
    addpath(eeglabpath);
    if exist( fullfile( eeglabpath, 'functions', 'adminfunc') ) ~= 7
        warning('EEGLAB subfolders not found');
    end;
end;

% determine file format
% ---------------------
fileformat = 'maclinux';
comp = computer;
try
    if strcmpi(comp(1:3), 'GLN') | strcmpi(comp(1:3), 'MAC') | strcmpi(comp(1:3), 'SOL')
        fileformat = 'maclinux';
    elseif strcmpi(comp(1:5), 'pcwin')
        fileformat = 'pcwin';
    end;
end;

% add paths
% ---------
if ~iseeglabdeployed2
    myaddpath( eeglabpath, 'eeg_checkset.m',   [ 'functions' filesep 'adminfunc'        ]);
    myaddpath( eeglabpath, 'memmapdata.m', 'functions');
    myaddpath( eeglabpath, 'readeetraklocs.m', [ 'functions' filesep 'sigprocfunc'      ]);
    myaddpath( eeglabpath, 'supergui.m',       [ 'functions' filesep 'guifunc'          ]);
    myaddpath( eeglabpath, 'pop_study.m',      [ 'functions' filesep 'studyfunc'        ]);
    myaddpath( eeglabpath, 'pop_loadbci.m',    [ 'functions' filesep 'popfunc'          ]);
    myaddpath( eeglabpath, 'statcond.m',       [ 'functions' filesep 'statistics'       ]);
    myaddpath( eeglabpath, 'timefreq.m',       [ 'functions' filesep 'timefreqfunc'     ]);
    myaddpath( eeglabpath, 'icademo.m',        [ 'functions' filesep 'miscfunc'         ]);
    myaddpath( eeglabpath, 'eeglab1020.ced',   [ 'functions' filesep 'resources'        ]);
    myaddpath( eeglabpath, 'startpane.m',      [ 'functions' filesep 'javachatfunc' ]);
    eeglab_options;
    
    % remove path to to fmrlab if neceecessary
    path_runica = fileparts(which('runica'));
    if length(path_runica) > 6 && strcmpi(path_runica(end-5:end), 'fmrlab')
        rmpath(path_runica);
    end;

    % add path if toolboxes are missing
    % ---------------------------------
    tmppath = path;
    signalpath = fullfile(eeglabpath, 'functions', 'octavefunc', 'signal');
    optimpath  = fullfile(eeglabpath, 'functions', 'octavefunc', 'optim');
    if ~license('test','signal_toolbox') || option_donotusetoolboxes
        warning('off', 'MATLAB:dispatcher:nameConflict');
        addpath( signalpath );
    else
        if ~isempty(findstr( signalpath, tmppath)) rmpath( signalpath ); end;
    end;
    if ~license('test','optim_toolbox') && ~ismatlab
        addpath( optimpath );
    else
        if ~isempty(findstr( optimpath, tmppath)) rmpath( optimpath ); end;
    end;
    
    myaddpath( eeglabpath, 'eegplugin_dipfit', 'plugins');
else
    eeglab_options;
end;

if nargin == 1 && strcmp(onearg, 'redraw')
    if evalin('base', 'exist(''EEG'')', '0') == 1
        evalin('base', 'eeg_global;');
    end;
end;
eeg_global;

% for the history function
% ------------------------
comtmp = 'warning off MATLAB:mir_warning_variable_used_as_function';
evalin('base'  , comtmp, '');
evalin('caller', comtmp, '');
    
evalin('base', 'eeg_global;');
if nargin < 1 | exist('EEG') ~= 1
	clear global EEG ALLEEG CURRENTSET ALLCOM LASTCOM STUDY;
    CURRENTSTUDY = 0;
	eeg_global;
	EEG = eeg_emptyset;
	eegh('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;');
    if ismatlab && get(0, 'screendepth') <= 8
        disp('Warning: screen color depth too low, some colors will be inaccurate in time-frequency plots');
    end;
end;

if nargin == 1
	if strcmp(onearg, 'versions')
        disp( [ 'EEGLAB v' eeg_getversion ] );
	elseif strcmp(onearg, 'nogui')
        if nargout < 1, clear ALLEEG; end; % do not return output var
        return;
	elseif strcmp(onearg, 'redraw')
        if ~ismatlab,return; end;
		W_MAIN = findobj('tag', 'EEGLAB');
		if ~isempty(W_MAIN)
			updatemenu;
            if nargout < 1, clear ALLEEG; end; % do not return output var
			return;
		else
			eegh('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab(''rebuild'');');
		end;
	elseif strcmp(onearg, 'besa');
		disp('Besa option deprecated. Download the BESA plugin to add the BESA menu.');
        eegh('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;');
	else
        eegh('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab(''rebuild'');');
	end;
else 
    onearg = 'rebuild';
end;
ALLCOM = ALLCOM;
try, eval('colordef white;'); catch end;

% default option folder
% ---------------------
if ~iseeglabdeployed2
    eeglab_options;
    fprintf('eeglab: options file is %s%seeg_options.m\n', homefolder, filesep);
end;

% checking strings
% ----------------
e_try             = 'try,';
e_catch           = 'catch, eeglab_error; LASTCOM= ''''; clear EEGTMP ALLEEGTMP STUDYTMP; end;';
nocheck           = e_try;
ret               = 'if ~isempty(LASTCOM), if LASTCOM(1) == -1, LASTCOM = ''''; return; end; end;';
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

storecall    = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);'');';
storenewcall = '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''study'', ~isempty(STUDY)+0); eegh(LASTCOM);';
storeallcall = [ 'if ~isempty(ALLEEG) & ~isempty(ALLEEG(1).data), ALLEEG = eeg_checkset(ALLEEG);' ...
                 'EEG = eeg_retrieve(ALLEEG, CURRENTSET); eegh(''ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_retrieve(ALLEEG, CURRENTSET);''); end;' ];

testeegtmp   =  'if exist(''EEGTMP'') == 1, EEG = EEGTMP; clear EEGTMP; end;'; % for backward compatibility
ifeeg        =  'if ~isempty(LASTCOM) & ~isempty(EEG),';
ifeegnh      =  'if ~isempty(LASTCOM) & ~isempty(EEG) & ~isempty(findstr(''='',LASTCOM)),';

% nh = no dataset history
% -----------------------
e_storeall_nh   = [e_catch 'eegh(LASTCOM);' ifeeg storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist_nh       = [e_catch 'eegh(LASTCOM);'];

% same as above but also save history in dataset
% ----------------------------------------------
e_newset        = [e_catch 'EEG = eegh(LASTCOM, EEG);' testeegtmp ifeeg   storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store         = [e_catch 'EEG = eegh(LASTCOM, EEG);' ifeegnh storecall    'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist          = [e_catch 'EEG = eegh(LASTCOM, EEG);'];
e_histdone      = [e_catch 'EEG = eegh(LASTCOM, EEG); if ~isempty(LASTCOM), disp(''Done.''); end;' ];

% study checking
% --------------
e_load_study    = [e_catch 'eegh(LASTCOM); if ~isempty(LASTCOM), ALLEEG = ALLEEGTMP; EEG = ALLEEG; CURRENTSET = [1:length(EEG)]; eegh(''CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];''); STUDY = STUDYTMP; CURRENTSTUDY = 1; disp(''Done.''); end; clear ALLEEGTMP STUDYTMP; eeglab(''redraw'');'];

% build structures for plugins
% ----------------------------
trystrs.no_check                 = e_try;
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
catchstrs.new_non_empty          = e_newset;
catchstrs.update_study           = e_load_study;

% create eeglab figure
% --------------------
javaobj = eeg_mainfig(onearg);
ptopoplot = fileparts(which('cbar'));

% detecting icalab
% ----------------
if exist('icalab')
    disp('ICALAB toolbox detected (algo. added to "run ICA" interface)');
end;

biosigflag = 0;
if iseeglabdeployed2
    biosigflag = 1;
else
    % adding all folders in external
    % ------------------------------
    disp(['Adding path to all EEGLAB functions']);
    p = which('eeglab.m');
    p = p(1:findstr(p,'eeglab.m')-1);
    if strcmpi(p, './') || strcmpi(p, '.\'), p = [ pwd filesep ]; end;
    dircontent  = dir([ p 'external' ]);
    dircontent  = { dircontent.name };
    for index = 1:length(dircontent)
        if dircontent{index}(1) ~= '.'
            if exist([p 'external' filesep dircontent{index}]) == 7
                addpath([p 'external' filesep dircontent{index}]);
                disp(['Adding path to eeglab' filesep 'external' filesep dircontent{index}]);
            end;
            if ~isempty(findstr('fieldtrip', lower(dircontent{index})))
                addpathexist([p 'external' filesep dircontent{index} filesep 'forward' ]);
                addpathexist([p 'external' filesep dircontent{index} filesep 'inverse' ]);
                addpathexist([p 'external' filesep dircontent{index} filesep 'utilities' ]);
                addpathexist([p 'external' filesep dircontent{index} ]);
                disp(['Adding path to eeglab' filesep 'external' filesep dircontent{index} ' subfolders' ]);
            end;
            if ~isempty(findstr('biosig', lower(dircontent{index})))
                addpathexist([p 'external' filesep dircontent{index} filesep 't200_FileAccess' ]);
                addpathexist([p 'external' filesep dircontent{index} filesep 't250_ArtifactPreProcessingQualityControl' ]);
                addpathexist([p 'external' filesep dircontent{index} filesep 'doc' ]);

                biosigflag = 1;
            end;
        end;
    end;

    % check for older version of Fieldtrip and presence of topoplot
    % -------------------------------------------------------------
    if ismatlab
        ptopoplot2 = fileparts(which('topoplot'));
        if ~strcmpi(ptopoplot, ptopoplot2),
            %disp('  Warning: duplicate function topoplot.m in Fieldtrip and EEGLAB');
            %disp('  EEGLAB function will prevail and call the Fieldtrip one when appropriate');
            addpath(ptopoplot);
        end;
    end;
end;

cb_importdata  = [ nocheck '[EEG LASTCOM] = pop_importdata;'   e_newset ];
cb_readegi     = [ nocheck '[EEG LASTCOM] = pop_readegi;'      e_newset ];
cb_readsegegi  = [ nocheck '[EEG LASTCOM] = pop_readsegegi;'   e_newset ];
cb_readegiepo  = [ nocheck '[EEG LASTCOM] = pop_importegimat;' e_newset ];
cb_loadbci     = [ nocheck '[EEG LASTCOM] = pop_loadbci;'      e_newset ];
cb_snapread    = [ nocheck '[EEG LASTCOM] = pop_snapread;'     e_newset ]; 
cb_loadcnt     = [ nocheck '[EEG LASTCOM] = pop_loadcnt;'      e_newset ]; 
cb_loadeeg     = [ nocheck '[EEG LASTCOM] = pop_loadeeg;'      e_newset ]; 
cb_biosig      = [ nocheck '[EEG LASTCOM] = pop_biosig; '      e_newset ]; 
cb_fileio      = [ nocheck '[EEG LASTCOM] = pop_fileio; '      e_newset ]; 
cb_fileio2     = [ nocheck '[EEG LASTCOM] = pop_fileiodir;'   e_newset ]; 

cb_importepoch = [ checkepoch   '[EEG LASTCOM] = pop_importepoch(EEG);'   e_store ];
cb_loaddat     = [ checkepoch   '[EEG LASTCOM]= pop_loaddat(EEG);'        e_store ]; 
cb_importevent = [ check        '[EEG LASTCOM] = pop_importevent(EEG);'   e_store ];
cb_chanevent   = [ check        '[EEG LASTCOM]= pop_chanevent(EEG);'      e_store ]; 
cb_importpres  = [ check        '[EEG LASTCOM]= pop_importpres(EEG);'     e_store ]; 
cb_importev2   = [ check        '[EEG LASTCOM]= pop_importev2(EEG);'      e_store ]; 
cb_export      = [ check        'LASTCOM = pop_export(EEG);'              e_histdone ];
cb_expica1     = [ check        'LASTCOM = pop_expica(EEG, ''weights'');' e_histdone ]; 
cb_expica2     = [ check        'LASTCOM = pop_expica(EEG, ''inv'');'     e_histdone ];
cb_expevents   = [ check        'LASTCOM = pop_expevents(EEG);'           e_histdone ];
cb_expdata     = [ check        'LASTCOM = pop_writeeeg(EEG);'            e_histdone ]; 

cb_loadset     = [ nocheck '[EEG LASTCOM] = pop_loadset;'                                e_newset];
cb_saveset     = [ check   '[EEG LASTCOM] = pop_saveset(EEG, ''savemode'', ''resave'');' e_store ];
cb_savesetas   = [ check   '[EEG LASTCOM] = pop_saveset(EEG);'                           e_store ];
cb_delset      = [ nocheck '[ALLEEG LASTCOM] = pop_delset(ALLEEG, -CURRENTSET);'         e_hist_nh 'eeglab redraw;' ];
cb_study1      = [ nocheck 'pop_stdwarn; [STUDYTMP ALLEEGTMP LASTCOM] = pop_study([], ALLEEG         , ''gui'', ''on'');' e_load_study]; 
cb_study2      = [ nocheck 'pop_stdwarn; [STUDYTMP ALLEEGTMP LASTCOM] = pop_study([], isempty(ALLEEG), ''gui'', ''on'');' e_load_study]; 
cb_studyerp    = [ nocheck 'pop_stdwarn; [STUDYTMP ALLEEGTMP LASTCOM] = pop_studyerp;' e_load_study]; 
cb_loadstudy   = [ nocheck 'pop_stdwarn; [STUDYTMP ALLEEGTMP LASTCOM] = pop_loadstudy;'                                   e_load_study]; 
cb_savestudy1  = [ check   '[STUDYTMP ALLEEGTMP LASTCOM] = pop_savestudy(STUDY, EEG, ''savemode'', ''resave'');'      e_load_study];
cb_savestudy2  = [ check   '[STUDYTMP ALLEEGTMP LASTCOM] = pop_savestudy(STUDY, EEG);'                                e_load_study];
cb_clearstudy  =           'LASTCOM = ''STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];''; eval(LASTCOM); eegh( LASTCOM ); eeglab redraw;';
cb_editoptions = [ nocheck 'if isfield(ALLEEG, ''nbchan''), LASTCOM = pop_editoptions(length([ ALLEEG.nbchan ]) >1);' ...
                                'else                            LASTCOM = pop_editoptions(0); end;'                  e_storeall_nh];

cb_saveh1      = [ nocheck 'LASTCOM = pop_saveh(EEG.history);' e_hist_nh];
cb_saveh2      = [ nocheck 'LASTCOM = pop_saveh(ALLCOM);'      e_hist_nh];
cb_runsc       = [ nocheck 'LASTCOM = pop_runscript;'          e_hist   ];
cb_quit        = [ 'close(gcf); disp(''To save the EEGLAB command history  >> pop_saveh(ALLCOM);'');' ...
                   'clear global EEG ALLEEG LASTCOM CURRENTSET;'];

cb_editset     = [ check      '[EEG LASTCOM] = pop_editset(EEG);'        e_store];
cb_editeventf  = [ checkevent '[EEG LASTCOM] = pop_editeventfield(EEG);' e_store];
cb_editeventv  = [ checkevent '[EEG LASTCOM] = pop_editeventvals(EEG);'  e_store];
cb_comments    = [ check      '[EEG.comments LASTCOM] =pop_comments(EEG.comments, ''About this dataset'');' e_store];
cb_chanedit    = [ 'disp(''IMPORTANT: After importing/modifying data channels, you must close'');' ...
                   'disp(''the channel editing window for the changes to take effect in EEGLAB.'');' ...
                   'disp(''TIP: Call this function directy from the prompt, ">> pop_chanedit([]);"'');' ...
                   'disp(''     to convert between channel location file formats'');' ...
                   '[EEG TMPINFO TMP LASTCOM] = pop_chanedit(EEG); if ~isempty(LASTCOM), EEG = eeg_checkset(EEG, ''chanlocsize'');' ...
                   'clear TMPINFO TMP; EEG = eegh(LASTCOM, EEG);' storecall 'end; eeglab(''redraw'');'];
cb_select      = [ check      '[EEG LASTCOM] = pop_select(EEG);'                     e_newset];
cb_rmdat       = [ checkevent '[EEG LASTCOM] = pop_rmdat(EEG);'                      e_newset];
cb_selectevent = [ checkevent '[EEG TMP LASTCOM] = pop_selectevent(EEG); clear TMP;' e_newset ];
cb_copyset     = [ check      '[ALLEEG EEG CURRENTSET LASTCOM] = pop_copyset(ALLEEG, CURRENTSET); eeglab(''redraw'');' e_hist_nh];
cb_mergeset    = [ check      '[EEG LASTCOM] = pop_mergeset(ALLEEG);' e_newset];

cb_resample    = [ check      '[EEG LASTCOM] = pop_resample(EEG);' e_newset];
cb_eegfilt     = [ check      '[EEG LASTCOM] = pop_eegfilt(EEG);'  e_newset];
cb_interp      = [ check      '[EEG LASTCOM] = pop_interp(EEG); '  e_newset];
cb_reref       = [ check      '[EEG LASTCOM] = pop_reref(EEG);'    e_newset];
cb_eegplot     = [ checkcont  '[LASTCOM] = pop_eegplot(EEG, 1);'   e_hist];
cb_epoch       = [ check      '[EEG tmp LASTCOM] = pop_epoch(EEG); clear tmp;' e_newset check '[EEG LASTCOM] = pop_rmbase(EEG);' e_store];
cb_rmbase      = [ check      '[EEG LASTCOM] = pop_rmbase(EEG);'   e_store];
cb_runica      = [ check      '[EEG LASTCOM] = pop_runica(EEG);'   e_store];
cb_subcomp     = [ checkica   '[EEG LASTCOM] = pop_subcomp(EEG);'  e_newset];
%cb_chanrej     = [ check      'pop_rejchan(EEG); LASTCOM = '''';'  e_hist];
cb_chanrej     = [ check      '[EEG tmp1 tmp2 LASTCOM] = pop_rejchan(EEG); clear tmp1 tmp2;'  e_hist];
cb_autorej     = [ check      '[EEG tmpp LASTCOM] = pop_autorej(EEG); clear tmpp;'  e_hist];

cb_rejmenu1    = [ check      'pop_rejmenu(EEG, 1); LASTCOM = '''';'    e_hist];
cb_eegplotrej1 = [ check      '[LASTCOM] = pop_eegplot(EEG, 1);'        e_hist];
cb_eegthresh1  = [ checkepoch '[TMP LASTCOM] = pop_eegthresh(EEG, 1); clear TMP;' e_hist];
cb_rejtrend1   = [ checkepoch '[EEG LASTCOM] = pop_rejtrend(EEG, 1);'   e_store];
cb_jointprob1  = [ checkepoch '[EEG LASTCOM] = pop_jointprob(EEG, 1);'  e_store];
cb_rejkurt1    = [ checkepoch '[EEG LASTCOM] = pop_rejkurt(EEG, 1);'    e_store];
cb_rejspec1    = [ checkepoch '[EEG Itmp LASTCOM] = pop_rejspec(EEG, 1); clear Itmp;' e_store];
cb_rejsup1     = [ checkepochica '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 1,1,1,1,1,1,1,1); eegh(LASTCOM);' ...
                     'LASTCOM = ''EEG.reject.icarejmanual = EEG.reject.rejglobal;''; eval(LASTCOM);' e_store ];
cb_rejsup2     = [ checkepoch '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 1,1,1,1,1,1,1,1); EEG = eegh(LASTCOM, EEG);' ...
                   '[EEG LASTCOM] = pop_rejepoch(EEG);' e_newset];

cb_selectcomps = [ checkicaplot  '[EEG LASTCOM] = pop_selectcomps(EEG);'  e_store];
cb_rejmenu2    = [ checkepochica 'pop_rejmenu(EEG, 0); LASTCOM ='''';'    e_hist];
cb_eegplotrej2 = [ checkica      '[LASTCOM] = pop_eegplot(EEG, 0);'       e_hist];
cb_eegthresh2  = [ checkepochica '[TMP LASTCOM] = pop_eegthresh(EEG, 0); clear TMP;' e_hist];
cb_rejtrend2   = [ checkepochica '[EEG LASTCOM] = pop_rejtrend(EEG, 0);'  e_store];
cb_jointprob2  = [ checkepochica '[EEG LASTCOM] = pop_jointprob(EEG, 0);' e_store];
cb_rejkurt2    = [ checkepochica '[EEG LASTCOM] = pop_rejkurt(EEG, 0);'   e_store];
cb_rejspec2    = [ checkepochica '[EEG Itmp LASTCOM] = pop_rejspec(EEG, 1); clear Itmp;'   e_store];
cb_rejsup3     = [ checkepochica '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); eegh(LASTCOM);' ...
                   'LASTCOM = ''EEG.reject.rejmanual = EEG.reject.rejglobal;''; eval(LASTCOM);' e_store ];
cb_rejsup4     = [ checkepochica '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); EEG = eegh(LASTCOM, EEG);' ...
                   '[EEG LASTCOM] = pop_rejepoch(EEG);'             e_newset ];

cb_topoblank1  = [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''',  ' ...
                   '''''electrodes'''', ''''labelpoint'''', ''''chaninfo'''', EEG.chaninfo);'']; eval(LASTCOM);' e_hist];
cb_topoblank2  = [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''',  ' ...
                   '''''electrodes'''', ''''numpoint'''', ''''chaninfo'''', EEG.chaninfo);'']; eval(LASTCOM);' e_hist];
cb_eegplot1    = [ check         'LASTCOM = pop_eegplot(EEG, 1, 1, 1);'   e_hist];
cb_spectopo1   = [ check         'LASTCOM = pop_spectopo(EEG, 1);'        e_hist];
cb_prop1       = [ checkplot     'LASTCOM = pop_prop(EEG,1);'             e_hist];
cb_erpimage1   = [ checkepoch    'LASTCOM = pop_erpimage(EEG, 1, eegh(''find'',''pop_erpimage(EEG,1''));' e_hist];
cb_timtopo     = [ checkplot     'LASTCOM = pop_timtopo(EEG);'            e_hist];
cb_plottopo    = [ check         'LASTCOM = pop_plottopo(EEG);'           e_hist];
cb_topoplot1   = [ checkplot     'LASTCOM = pop_topoplot(EEG, 1);'        e_hist];
cb_headplot1   = [ checkplot     '[EEG LASTCOM] = pop_headplot(EEG, 1);'  e_store];
cb_comperp1    = [ checkepoch    'LASTCOM = pop_comperp(ALLEEG);'         e_hist];

cb_eegplot2    = [ checkica      '[LASTCOM] = pop_eegplot(EEG, 0, 1, 1);' e_hist];
cb_spectopo2   = [ checkicaplot  'LASTCOM = pop_spectopo(EEG, 0);'        e_hist];
cb_topoplot2   = [ checkicaplot  'LASTCOM = pop_topoplot(EEG, 0);'        e_hist];
cb_headplot2   = [ checkicaplot  '[EEG LASTCOM] = pop_headplot(EEG, 0);'  e_store];
cb_prop2       = [ checkicaplot  'LASTCOM = pop_prop(EEG,0);'             e_hist];
cb_erpimage2   = [ checkepochica 'LASTCOM = pop_erpimage(EEG, 0, eegh(''find'',''pop_erpimage(EEG,0''));' e_hist];
cb_envtopo1    = [ checkica      'LASTCOM = pop_envtopo(EEG);'            e_hist];
cb_envtopo2    = [ checkica      'if length(ALLEEG) == 1, error(''Need at least 2 datasets''); end; LASTCOM = pop_envtopo(ALLEEG);' e_hist];
cb_plotdata2   = [ checkepochica '[tmpeeg LASTCOM] = pop_plotdata(EEG, 0); clear tmpeeg;' e_hist];
cb_comperp2    = [ checkepochica 'LASTCOM = pop_comperp(ALLEEG, 0);'      e_hist];

cb_signalstat1 = [ check         'LASTCOM = pop_signalstat(EEG, 1);' e_hist];
cb_signalstat2 = [ checkica      'LASTCOM = pop_signalstat(EEG, 0);' e_hist];
cb_eventstat   = [ checkevent    'LASTCOM = pop_eventstat(EEG);'     e_hist];
cb_timef1      = [ check         'LASTCOM = pop_newtimef(EEG, 1, eegh(''find'',''pop_newtimef(EEG,1''));'  e_hist];
cb_crossf1     = [ check         'LASTCOM = pop_newcrossf(EEG, 1,eegh(''find'',''pop_newcrossf(EEG,1''));' e_hist];
cb_timef2      = [ checkica      'LASTCOM = pop_newtimef(EEG, 0, eegh(''find'',''pop_newtimef(EEG,0''));'  e_hist];
cb_crossf2     = [ checkica      'LASTCOM = pop_newcrossf(EEG, 0,eegh(''find'',''pop_newcrossf(EEG,0''));' e_hist];

cb_study3      = [ nocheck '[STUDYTMP ALLEEGTMP LASTCOM] = pop_study(STUDY, ALLEEG, ''gui'', ''on'');'  e_load_study];
cb_studydesign = [ nocheck '[STUDYTMP LASTCOM] = pop_studydesign(STUDY, ALLEEG); ALLEEGTMP = ALLEEG;'   e_load_study];     
cb_precomp     = [ nocheck '[STUDYTMP ALLEEGTMP LASTCOM] = pop_precomp(STUDY, ALLEEG);'                 e_load_study];
cb_chanplot    = [ nocheck '[STUDYTMP LASTCOM] = pop_chanplot(STUDY, ALLEEG); ALLEEGTMP=ALLEEG;'        e_load_study];
cb_precomp2    = [ nocheck '[STUDYTMP ALLEEGTMP LASTCOM] = pop_precomp(STUDY, ALLEEG, ''components'');' e_load_study];
cb_preclust    = [ nocheck '[STUDYTMP ALLEEGTMP LASTCOM] = pop_preclust(STUDY, ALLEEG);'                e_load_study];
cb_clust       = [ nocheck '[STUDYTMP ALLEEGTMP LASTCOM] = pop_clust(STUDY, ALLEEG);'                   e_load_study];
cb_clustedit   = [ nocheck 'ALLEEGTMP = ALLEEG; [STUDYTMP LASTCOM] = pop_clustedit(STUDY, ALLEEG);'     e_load_study];
% 
% % add STUDY plugin menus
% if exist('eegplugin_stderpimage')
%     structure.uilist = { { } ...
%         {'style' 'pushbutton' 'string' 'Plot ERPimage'    'Callback' 'stderpimageplugin_plot(''onecomp'', gcf);' } { }  ...
%         {'style' 'pushbutton' 'string' 'Plot ERPimage(s)' 'Callback' 'stderpimageplugin_plot(''oneclust'', gcf);' } };
%     structure.geometry = { [1] [1 0.3 1] };
%     arg = vararg2str( { structure } );
%     cb_clustedit   = [ nocheck 'ALLEEGTMP = ALLEEG; [STUDYTMP LASTCOM] = pop_clustedit(STUDY, ALLEEG, [], ' arg ');'     e_load_study];
% end;

% menu definition
% --------------- 
if ismatlab
    % defaults
    % --------
    % startup:on
    % study:off
    % chanloc:off
    % epoch:on
    % continuous:on
    
    on          = 'study:on';
    onnostudy   = '';
    ondata      = 'startup:off';
    onepoch     = 'startup:off;continuous:off';
    ondatastudy = 'startup:off;study:on';
    onchannel   = 'startup:off;chanloc:on';
    onepochchan = 'startup:off;continuous:off;chanloc:on';
    onstudy     = 'startup:off;epoch:off;continuous:off;study:on';
    
    W_MAIN = findobj('tag', 'EEGLAB');
    EEGUSERDAT = get(W_MAIN, 'userdata');
    set(W_MAIN, 'MenuBar', 'none');
    file_m   = uimenu( W_MAIN,   'Label', 'File'                                    , 'userdata', on);
    import_m = uimenu( file_m,   'Label', 'Import data'                             , 'userdata', onnostudy); 
    neuro_m  = uimenu( import_m, 'Label', 'Using EEGLAB functions and plugins'      , 'tag', 'import data' , 'userdata', onnostudy); 
    epoch_m  = uimenu( file_m,   'Label', 'Import epoch info', 'tag', 'import epoch', 'userdata', onepoch); 
    event_m  = uimenu( file_m,   'Label', 'Import event info', 'tag', 'import event', 'userdata', ondata); 
    exportm  = uimenu( file_m,   'Label', 'Export'           , 'tag', 'export'      , 'userdata', ondata); 
    edit_m   = uimenu( W_MAIN,   'Label', 'Edit'                                    , 'userdata', ondata);
    tools_m  = uimenu( W_MAIN,   'Label', 'Tools',             'tag', 'tools'       , 'userdata', ondatastudy);
    plot_m   = uimenu( W_MAIN,   'Label', 'Plot',              'tag', 'plot'        , 'userdata', ondata);
    loc_m    = uimenu( plot_m,   'Label', 'Channel locations'                       , 'userdata', onchannel);
    std_m    = uimenu( W_MAIN,   'Label', 'Study', 'tag', 'study'                   , 'userdata', onstudy);
    set_m    = uimenu( W_MAIN,   'Label', 'Datasets'                                , 'userdata', ondatastudy);
    help_m   = uimenu( W_MAIN,   'Label', 'Help'                                    , 'userdata', on);

    uimenu( neuro_m, 'Label', 'From ASCII/float file or Matlab array' , 'CallBack', cb_importdata);
    uimenu( neuro_m, 'Label', 'From Netstation .mff (FILE-IO toolbox)', 'CallBack', cb_fileio2,    'Separator', 'on'); 
    uimenu( neuro_m, 'Label', 'From Netstation binary simple file'    , 'CallBack', cb_readegi); 
    uimenu( neuro_m, 'Label', 'From Multiple seg. Netstation files'   , 'CallBack', cb_readsegegi); 
    uimenu( neuro_m, 'Label', 'From Netstation Matlab files'          , 'CallBack', cb_readegiepo); 
    uimenu( neuro_m, 'Label', 'From BCI2000 ASCII file'               , 'CallBack', cb_loadbci,    'Separator', 'on'); 
    uimenu( neuro_m, 'Label', 'From Snapmaster .SMA file'             , 'CallBack', cb_snapread,   'Separator', 'on'); 
    uimenu( neuro_m, 'Label', 'From Neuroscan .CNT file'              , 'CallBack', cb_loadcnt,    'Separator', 'on'); 
    uimenu( neuro_m, 'Label', 'From Neuroscan .EEG file'              , 'CallBack', cb_loadeeg); 

    % BIOSIG MENUS
    % ------------
    if biosigflag
        uimenu( neuro_m, 'Label', 'From Biosemi BDF file (BIOSIG toolbox)', 'CallBack' , cb_biosig, 'Separator', 'on'); 
        uimenu( neuro_m, 'Label', 'From EDF/EDF+/GDF files (BIOSIG toolbox)', 'CallBack', cb_biosig); 
    end;

    uimenu( epoch_m, 'Label', 'From Matlab array or ASCII file'       , 'CallBack', cb_importepoch);
    uimenu( epoch_m, 'Label', 'From Neuroscan .DAT file'              , 'CallBack', cb_loaddat); 
    uimenu( event_m, 'Label', 'From Matlab array or ASCII file'       , 'CallBack', cb_importevent);
    uimenu( event_m, 'Label', 'From data channel'                     , 'CallBack', cb_chanevent); 
    uimenu( event_m, 'Label', 'From Presentation .LOG file'           , 'CallBack', cb_importpres); 
    uimenu( event_m, 'Label', 'From E-Prime ASCII (text) file'        , 'CallBack', cb_importevent);
    uimenu( event_m, 'Label', 'From Neuroscan .ev2 file'              , 'CallBack', cb_importev2); 
    uimenu( exportm, 'Label', 'Data and ICA activity to text file'    , 'CallBack', cb_export);
    uimenu( exportm, 'Label', 'Weight matrix to text file'            , 'CallBack', cb_expica1); 
    uimenu( exportm, 'Label', 'Inverse weight matrix to text file'    , 'CallBack', cb_expica2);
    uimenu( exportm, 'Label', 'Events to text file'                   , 'CallBack', cb_expevents);
    uimenu( exportm, 'Label', 'Data to EDF/BDF/GDF file'              , 'CallBack', cb_expdata, 'separator', 'on'); 

    uimenu( file_m, 'Label', 'Load existing dataset'                  , 'userdata', onnostudy,   'CallBack', cb_loadset, 'Separator', 'on'); 
    uimenu( file_m, 'Label', 'Save current dataset(s)'                , 'userdata', ondatastudy, 'CallBack', cb_saveset);
    uimenu( file_m, 'Label', 'Save current dataset as'                , 'userdata', ondata,      'CallBack', cb_savesetas);
    uimenu( file_m, 'Label', 'Clear dataset(s)'                       , 'userdata', ondata,      'CallBack', cb_delset);

    std2_m = uimenu( file_m, 'Label', 'Create study'                  , 'userdata', on     , 'Separator', 'on'); 
    uimenu( std2_m,  'Label', 'Using all loaded datasets'             , 'userdata', ondata , 'Callback', cb_study1); 
    uimenu( std2_m,  'Label', 'Browse for datasets'                   , 'userdata', on     , 'Callback', cb_study2); 
    uimenu( std2_m,  'Label', 'Simple ERP STUDY'                      , 'userdata', on     , 'Callback', cb_studyerp); 

    uimenu( file_m, 'Label', 'Load existing study'                    , 'userdata', on     , 'CallBack', cb_loadstudy,'Separator', 'on' ); 
    uimenu( file_m, 'Label', 'Save current study'                     , 'userdata', onstudy, 'CallBack', cb_savestudy1);
    uimenu( file_m, 'Label', 'Save current study as'                  , 'userdata', onstudy, 'CallBack', cb_savestudy2);
    uimenu( file_m, 'Label', 'Clear study'                            , 'userdata', ondatastudy, 'CallBack', cb_clearstudy);
    uimenu( file_m, 'Label', 'Memory and other options'               , 'userdata', on     , 'CallBack', cb_editoptions, 'Separator', 'on');

    hist_m = uimenu( file_m, 'Label', 'History scripts'               , 'userdata', on     , 'Separator', 'on');
    uimenu( hist_m, 'Label', 'Save dataset history script'            , 'userdata', ondata     , 'CallBack', cb_saveh1);
    uimenu( hist_m, 'Label', 'Save session history script'            , 'userdata', ondatastudy, 'CallBack', cb_saveh2);    
    uimenu( hist_m, 'Label', 'Run script'                             , 'userdata', on         , 'CallBack', cb_runsc);    

    uimenu( file_m, 'Label', 'Quit'                                   , 'userdata', on     , 'CallBack', cb_quit, 'Separator', 'on');

    uimenu( edit_m, 'Label', 'Dataset info'                           , 'userdata', ondata, 'CallBack', cb_editset);
    uimenu( edit_m, 'Label', 'Event fields'                           , 'userdata', ondata, 'CallBack', cb_editeventf);
    uimenu( edit_m, 'Label', 'Event values'                           , 'userdata', ondata, 'CallBack', cb_editeventv);
    uimenu( edit_m, 'Label', 'About this dataset'                     , 'userdata', ondata, 'CallBack', cb_comments);
    uimenu( edit_m, 'Label', 'Channel locations'                      , 'userdata', ondata, 'CallBack', cb_chanedit);
    uimenu( edit_m, 'Label', 'Select data'                            , 'userdata', ondata, 'CallBack', cb_select, 'Separator', 'on');
    uimenu( edit_m, 'Label', 'Select data using events'               , 'userdata', ondata, 'CallBack', cb_rmdat);
    uimenu( edit_m, 'Label', 'Select epochs or events'                , 'userdata', ondata, 'CallBack', cb_selectevent);
    uimenu( edit_m, 'Label', 'Copy current dataset'                   , 'userdata', ondata, 'CallBack', cb_copyset, 'Separator', 'on');
    uimenu( edit_m, 'Label', 'Append datasets'                        , 'userdata', ondata, 'CallBack', cb_mergeset);
    uimenu( edit_m, 'Label', 'Delete dataset(s)'                      , 'userdata', ondata, 'CallBack', cb_delset);

    uimenu( tools_m, 'Label', 'Change sampling rate'                  , 'userdata', ondatastudy, 'CallBack', cb_resample);

    filter_m = uimenu( tools_m, 'Label', 'Filter the data'            , 'userdata', ondatastudy, 'tag', 'filter');
    uimenu( filter_m, 'Label', 'Basic FIR filter'                     , 'userdata', ondatastudy, 'CallBack', cb_eegfilt);

    uimenu( tools_m, 'Label', 'Re-reference'                          , 'userdata', ondata, 'CallBack', cb_reref);
    uimenu( tools_m, 'Label', 'Interpolate electrodes'                , 'userdata', ondata, 'CallBack', cb_interp);
    uimenu( tools_m, 'Label', 'Reject continuous data by eye'         , 'userdata', ondata, 'CallBack', cb_eegplot);
    uimenu( tools_m, 'Label', 'Extract epochs'                        , 'userdata', ondata, 'CallBack', cb_epoch, 'Separator', 'on');
    uimenu( tools_m, 'Label', 'Remove baseline'                       , 'userdata', ondatastudy, 'CallBack', cb_rmbase);
    uimenu( tools_m, 'Label', 'Run ICA'                               , 'userdata', ondatastudy, 'CallBack', cb_runica, 'foregroundcolor', 'b', 'Separator', 'on');
    uimenu( tools_m, 'Label', 'Remove components'                     , 'userdata', ondata, 'CallBack', cb_subcomp);
    uimenu( tools_m, 'Label', 'Automatic channel rejection'           , 'userdata', ondata, 'CallBack', cb_chanrej, 'Separator', 'on');
    uimenu( tools_m, 'Label', 'Automatic epoch rejection'             , 'userdata', onepoch, 'CallBack', cb_autorej);
    rej_m1 = uimenu( tools_m, 'Label', 'Reject data epochs'           , 'userdata', onepoch);
    rej_m2 = uimenu( tools_m, 'Label', 'Reject data using ICA'        , 'userdata', ondata );

    uimenu( rej_m1, 'Label', 'Reject data (all methods)'              , 'userdata', onepoch, 'CallBack', cb_rejmenu1);
    uimenu( rej_m1, 'Label', 'Reject by inspection'                   , 'userdata', onepoch, 'CallBack', cb_eegplotrej1);
    uimenu( rej_m1, 'Label', 'Reject extreme values'                  , 'userdata', onepoch, 'CallBack', cb_eegthresh1);
    uimenu( rej_m1, 'Label', 'Reject by linear trend/variance'        , 'userdata', onepoch, 'CallBack', cb_rejtrend1);
    uimenu( rej_m1, 'Label', 'Reject by probability'                  , 'userdata', onepoch, 'CallBack', cb_jointprob1);
    uimenu( rej_m1, 'Label', 'Reject by kurtosis'                     , 'userdata', onepoch, 'CallBack', cb_rejkurt1);
    uimenu( rej_m1, 'Label', 'Reject by spectra'                      , 'userdata', onepoch, 'CallBack', cb_rejspec1);
    uimenu( rej_m1, 'Label', 'Export marks to ICA reject'             , 'userdata', onepoch, 'CallBack', cb_rejsup1, 'separator', 'on');
    uimenu( rej_m1, 'Label', 'Reject marked epochs'                   , 'userdata', onepoch, 'CallBack', cb_rejsup2, 'separator', 'on', 'foregroundcolor', 'b');
    uimenu( rej_m2, 'Label', 'Reject components by map'               , 'userdata', ondata , 'CallBack', cb_selectcomps);
    uimenu( rej_m2, 'Label', 'Reject data (all methods)'              , 'userdata', onepoch, 'CallBack', cb_rejmenu2, 'Separator', 'on');
    uimenu( rej_m2, 'Label', 'Reject by inspection'                   , 'userdata', onepoch, 'CallBack', cb_eegplotrej2);
    uimenu( rej_m2, 'Label', 'Reject extreme values'                  , 'userdata', onepoch, 'CallBack', cb_eegthresh2);
    uimenu( rej_m2, 'Label', 'Reject by linear trend/variance'        , 'userdata', onepoch, 'CallBack', cb_rejtrend2);
    uimenu( rej_m2, 'Label', 'Reject by probability'                  , 'userdata', onepoch, 'CallBack', cb_jointprob2);
    uimenu( rej_m2, 'Label', 'Reject by kurtosis'                     , 'userdata', onepoch, 'CallBack', cb_rejkurt2);
    uimenu( rej_m2, 'Label', 'Reject by spectra'                      , 'userdata', onepoch, 'CallBack', cb_rejspec2);
    uimenu( rej_m2, 'Label', 'Export marks to data reject'            , 'userdata', onepoch, 'CallBack', cb_rejsup3, 'separator', 'on');
    uimenu( rej_m2, 'Label', 'Reject marked epochs'                   , 'userdata', onepoch, 'CallBack', cb_rejsup4, 'separator', 'on', 'foregroundcolor', 'b');

    uimenu( loc_m,  'Label', 'By name'                                , 'userdata', onchannel, 'CallBack', cb_topoblank1);
    uimenu( loc_m,  'Label', 'By number'                              , 'userdata', onchannel, 'CallBack', cb_topoblank2);
    uimenu( plot_m, 'Label', 'Channel data (scroll)'                  , 'userdata', ondata , 'CallBack', cb_eegplot1, 'Separator', 'on');
    uimenu( plot_m, 'Label', 'Channel spectra and maps'               , 'userdata', ondata , 'CallBack', cb_spectopo1);
    uimenu( plot_m, 'Label', 'Channel properties'                     , 'userdata', ondata , 'CallBack', cb_prop1);
    uimenu( plot_m, 'Label', 'Channel ERP image'                      , 'userdata', onepoch, 'CallBack', cb_erpimage1);

    ERP_m = uimenu( plot_m, 'Label', 'Channel ERPs'                   , 'userdata', onepoch);
    uimenu( ERP_m,  'Label', 'With scalp maps'                        , 'CallBack', cb_timtopo);
    uimenu( ERP_m,  'Label', 'In scalp/rect. array'                   , 'CallBack', cb_plottopo);

    topo_m = uimenu( plot_m, 'Label', 'ERP map series'                , 'userdata', onepochchan);
    uimenu( topo_m, 'Label', 'In 2-D'                                 , 'CallBack', cb_topoplot1);
    uimenu( topo_m, 'Label', 'In 3-D'                                 , 'CallBack', cb_headplot1);
    uimenu( plot_m, 'Label', 'Sum/Compare ERPs'                       , 'userdata', onepoch, 'CallBack', cb_comperp1);

    uimenu( plot_m, 'Label', 'Component activations (scroll)'         , 'userdata', ondata , 'CallBack', cb_eegplot2,'Separator', 'on');
    uimenu( plot_m, 'Label', 'Component spectra and maps'             , 'userdata', ondata , 'CallBack', cb_spectopo2);

    tica_m = uimenu( plot_m, 'Label', 'Component maps'                , 'userdata', onchannel);
    uimenu( tica_m, 'Label', 'In 2-D'                                 , 'CallBack', cb_topoplot2);
    uimenu( tica_m, 'Label', 'In 3-D'                                 , 'CallBack', cb_headplot2);
    uimenu( plot_m, 'Label', 'Component properties'                   , 'userdata', ondata , 'CallBack', cb_prop2);
    uimenu( plot_m, 'Label', 'Component ERP image'                    , 'userdata', onepoch, 'CallBack', cb_erpimage2);

    ERPC_m = uimenu( plot_m, 'Label', 'Component ERPs'                , 'userdata', onepoch);
    uimenu( ERPC_m, 'Label', 'With component maps'                    , 'CallBack', cb_envtopo1);
    uimenu( ERPC_m, 'Label', 'With comp. maps (compare)'              , 'CallBack', cb_envtopo2);
    uimenu( ERPC_m, 'Label', 'In rectangular array'                   , 'CallBack', cb_plotdata2);
    uimenu( plot_m, 'Label', 'Sum/Compare comp. ERPs'                 , 'userdata', onepoch, 'CallBack', cb_comperp2);

    stat_m = uimenu( plot_m, 'Label', 'Data statistics', 'Separator', 'on', 'userdata', ondata );
    uimenu( stat_m, 'Label', 'Channel statistics'                     , 'CallBack', cb_signalstat1);
    uimenu( stat_m, 'Label', 'Component statistics'                   , 'CallBack', cb_signalstat2);
    uimenu( stat_m, 'Label', 'Event statistics'                       , 'CallBack', cb_eventstat);

    spec_m = uimenu( plot_m, 'Label', 'Time-frequency transforms', 'Separator', 'on', 'userdata', ondata);
    uimenu( spec_m, 'Label', 'Channel time-frequency'                 , 'CallBack', cb_timef1);
    uimenu( spec_m, 'Label', 'Channel cross-coherence'                , 'CallBack', cb_crossf1);
    uimenu( spec_m, 'Label', 'Component time-frequency'               , 'CallBack', cb_timef2,'Separator', 'on');     
    uimenu( spec_m, 'Label', 'Component cross-coherence'              , 'CallBack', cb_crossf2);

    uimenu( std_m,  'Label', 'Edit study info'                        , 'userdata', onstudy, 'CallBack', cb_study3);
    uimenu( std_m,  'Label', 'Select/Edit study design(s)'            , 'userdata', onstudy, 'CallBack', cb_studydesign);
    uimenu( std_m,  'Label', 'Precompute channel measures'            , 'userdata', onstudy, 'CallBack', cb_precomp, 'separator', 'on');
    uimenu( std_m,  'Label', 'Plot channel measures'                  , 'userdata', onstudy, 'CallBack', cb_chanplot);
    uimenu( std_m,  'Label', 'Precompute component measures'          , 'userdata', onstudy, 'CallBack', cb_precomp2, 'separator', 'on');
    clust_m = uimenu( std_m, 'Label', 'PCA clustering (original)'     , 'userdata', onstudy);
    uimenu( clust_m,  'Label', 'Build preclustering array'            , 'userdata', onstudy, 'CallBack', cb_preclust);
    uimenu( clust_m,  'Label', 'Cluster components'                   , 'userdata', onstudy, 'CallBack', cb_clust);
    uimenu( std_m,  'Label', 'Edit/plot clusters'                     , 'userdata', onstudy, 'CallBack', cb_clustedit);

    if ~iseeglabdeployed2
        newerVersionMenu = uimenu( help_m, 'Label', 'Upgrade to the Latest Version'          , 'userdata', on, 'ForegroundColor', [0.6 0 0]);
        uimenu( help_m, 'Label', 'About EEGLAB'                           , 'userdata', on, 'CallBack', 'pophelp(''eeglab'');');
        uimenu( help_m, 'Label', 'About EEGLAB help'                      , 'userdata', on, 'CallBack', 'pophelp(''eeg_helphelp'');');
        uimenu( help_m, 'Label', 'EEGLAB menus'                           , 'userdata', on, 'CallBack', 'eeg_helpmenu;','separator','on');

        help_1 = uimenu( help_m, 'Label', 'EEGLAB functions');
        uimenu( help_1, 'Label', 'Admin functions'                        , 'userdata', on, 'Callback', 'eeg_helpadmin;');	
        uimenu( help_1, 'Label', 'Interactive pop_ functions'             , 'userdata', on, 'Callback', 'eeg_helppop;');	
        uimenu( help_1, 'Label', 'Signal processing functions'            , 'userdata', on, 'Callback', 'eeg_helpsigproc;');	

        uimenu( help_m, 'Label', 'EEGLAB license'                         , 'userdata', on, 'CallBack', 'pophelp(''eeglablicense.txt'', 1);');
    else
        uimenu( help_m, 'Label', 'About EEGLAB'                           , 'userdata', on, 'CallBack', 'abouteeglab;');
        uimenu( help_m, 'Label', 'EEGLAB license'                         , 'userdata', on, 'CallBack', 'pophelp(''eeglablicense.txt'', 1);');
    end;

    uimenu( help_m, 'Label', 'Web tutorial'                               , 'userdata', on, 'CallBack', 'tutorial;', 'Separator', 'on');
    uimenu( help_m, 'Label', 'Email the EEGLAB team'                      , 'userdata', on, 'CallBack', 'web(''mailto:eeglab@sccn.ucsd.edu'');');
end;

if iseeglabdeployed2
    disp('Adding FIELDTRIP toolbox functions');
    disp('Adding BIOSIG toolbox functions');
    disp('Adding FILE-IO toolbox functions');
    funcname = {  'eegplugin_VisEd' ...
                  'eegplugin_eepimport' ...
                  'eegplugin_bdfimport' ...
                  'eegplugin_brainmovie' ...
                  'eegplugin_bva_io' ...
                  'eegplugin_ctfimport' ...
                  'eegplugin_dipfit' ...
                  'eegplugin_erpssimport' ...
                  'eegplugin_fmrib' ...
                  'eegplugin_iirfilt' ...
                  'eegplugin_ascinstep' ...
                  'eegplugin_loreta' ...
                  'eegplugin_miclust' ...
                  'eegplugin_4dneuroimaging' };
    for indf = 1:length(funcname)
        try 
            vers = feval(funcname{indf}, gcf, trystrs, catchstrs);
            disp(['EEGLAB: adding "' vers '" plugin' ]);  
        catch
            feval(funcname{indf}, gcf, trystrs, catchstrs);
            disp(['EEGLAB: adding plugin function "' funcname{indf} '"' ]);   
        end;
    end;
else
    % looking for eeglab plugins
    % --------------------------
    p = which('eeglab.m');
    p = p(1:findstr(p,'eeglab.m')-1);
    if strcmpi(p, './') || strcmpi(p, '.\'), p = [ pwd filesep ]; end;
    dircontent  = dir(fullfile(p, 'plugins'));
    dircontent  = { dircontent.name };

    % scan plugin folder
    % ------------------
    for index = 1:length(dircontent)

        % find function
        % -------------
        funcname = '';
        if exist([p 'plugins' filesep dircontent{index}]) == 7
            if ~strcmpi(dircontent{index}, '.') & ~strcmpi(dircontent{index}, '..')
                addpath([ p 'plugins' filesep dircontent{index} ]);
                tmpdir = dir([ p 'plugins' filesep dircontent{index}]);
                for tmpind = 1:length(tmpdir)
                    % find plugin function in subfolder
                    % ---------------------------------
                    if ~isempty(findstr(tmpdir(tmpind).name, 'eegplugin')) & tmpdir(tmpind).name(end) == 'm'
                        funcname = tmpdir(tmpind).name(1:end-2);
                        tmpind = length(tmpdir);
                    end;

                    % special case of eeglab subfolder (for BIOSIG)
                    % --------------------------------
                    if strcmpi(tmpdir(tmpind).name, 'eeglab')
                        addpath([ p 'plugins' filesep dircontent{index} filesep 'eeglab' ],'-end');
                        tmpdir2 = dir([ p 'plugins' filesep dircontent{index} filesep 'eeglab' ]);
                        for tmpind2 = 1:length(tmpdir2)
                            if ~isempty(findstr(tmpdir2(tmpind2).name, 'eegplugin')) ...
                                    & tmpdir2(tmpind2).name(end) == 'm'
                                funcname = tmpdir2(tmpind2).name(1:end-2);
                                tmpind2  = length(tmpdir2);
                                tmpind   = length(tmpdir);
                            end;
                        end;
                    end;
                end;
            end;
        else 
            if ~isempty(findstr(dircontent{index}, 'eegplugin')) & dircontent{index}(end) == 'm'
                funcname = dircontent{index}(1:end-2); % remove .m
            end;
        end;

        % execute function
        % ----------------
        if ~isempty(funcname)
            vers = ''; % version
            try,
                eval( [ 'vers =' funcname '(gcf, trystrs, catchstrs);' ]);
                disp(['EEGLAB: adding "' vers '" plugin (see >> help ' funcname ')' ]);    
           catch
                try,
                    eval( [ funcname '(gcf, trystrs, catchstrs)' ]);
                    disp(['EEGLAB: adding plugin function "' funcname '"' ]);    
                catch
                    disp([ 'EEGLAB: error while adding plugin "' funcname '"' ] ); 
                    disp([ '   ' lasterr] );
                end;
            end;
        end;
    end;
end; % iseeglabdeployed2

if ~ismatlab, return; end;
% add other import ...
% --------------------
cb_others = [ 'warndlg2(strvcat(''Several EEGLAB plugins (not included by default) are available to import specific'',' ...
                               '''data formats. To download plugins go to http://sccn.ucsd.edu/wiki/EEGLAB_Plugins'',' ...
                               '''  '',' ...
                               '''The FILEIO and BIOSIG toolboxes interface also allow to import in EEGLAB a wide variety'',' ...
                               '''of EEG/MEG data file formats -- including those imported using EEGLAB native menus'',' ...
                               '''(see www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:dataformat (FILEIO) and'',' ...
                               '''biosig.sourceforge.net/SupportedSystems.html (BIOSIG) for supported file formats).'',' ...
                               ''' ''));' ];
if exist('ft_chantype')
    uimenu( import_m, 'Label', 'Using the FILE-IO interface', 'CallBack', cb_fileio, 'separator', 'on'); 
end;
if biosigflag
    uimenu( import_m, 'Label', 'Using the BIOSIG interface' , 'CallBack', cb_biosig); 
end;
uimenu( import_m, 'Label', 'Troubleshooting data formats...', 'CallBack', cb_others);    

% changing plugin menu color
% --------------------------
fourthsub_m = findobj('parent', tools_m);
plotsub_m   = findobj('parent', plot_m);
importsub_m = findobj('parent', neuro_m);
epochsub_m  = findobj('parent', epoch_m);
eventsub_m  = findobj('parent', event_m);    
editsub_m   = findobj('parent', edit_m);
exportsub_m = findobj('parent', exportm);
filter_m    = findobj('parent', filter_m);
icadefs; % containing PLUGINMENUCOLOR
if length(fourthsub_m) > 11, set(fourthsub_m(1:end-11), 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(plotsub_m)   > 17, set(plotsub_m  (1:end-17), 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(importsub_m) > 9,  set(importsub_m(1:end-9) , 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(epochsub_m ) > 3 , set(epochsub_m (1:end-3 ), 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(eventsub_m ) > 4 , set(eventsub_m (1:end-4 ), 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(exportsub_m) > 4 , set(exportsub_m(1:end-4 ), 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(editsub_m)   > 10, set(editsub_m(  1:end-10), 'foregroundcolor', PLUGINMENUCOLOR); end;
if length(filter_m)    > 3 , set(filter_m   (1:end-1 ), 'foregroundcolor', PLUGINMENUCOLOR); end;

EEGMENU = uimenu( set_m, 'Label', '------', 'Enable', 'off');
eval('set(W_MAIN, ''userdat'', { EEGUSERDAT{1} EEGMENU javaobj });');
eeglab('redraw');
if nargout < 1
    clear ALLEEG;
end;


%% automatic updater
try
    [dummy eeglabVersionNumber currentReleaseDateString] = eeg_getversion;
    eeglabUpdater = up.updater(eeglabVersionNumber, 'http://sccn.ucsd.edu/eeglab/updater/latest_version.php', 'EEGLAB', currentReleaseDateString);
    
    % create a new GUI item (e.g. under Help)
    %newerVersionMenu = uimenu(help_m, 'Label', 'Upgrade to the Latest Version', 'visible', 'off', 'userdata', 'startup:on;study:on');
    eeglabUpdater.menuItemHandle = newerVersionMenu;
    
    % set the callback to bring up the updater GUI
    icadefs; % for getting background color
    eeglabFolder = fileparts(which('eeglab.m'));
    eeglabUpdater.menuItemCallback = {@command_on_update_menu_click, eeglabUpdater, eeglabFolder, true, BACKEEGLABCOLOR};

    % place it in the base workspace.
    assignin('base', 'eeglabUpdater', eeglabUpdater);
    
    eeglabUpdater.checkForNewVersion({'eeglab_event' 'setup'});
    
    if eeglabUpdater.newerVersionIsAvailable
        fprintf('\nA newer version (%s, %d days newer) of EEGLAB is available to download at: \n%s\n\n', num2str(eeglabUpdater.latestVersionNumber), eeglabUpdater.daysCurrentVersionIsOlder, eeglabUpdater.downloadUrl);
        
        % make the Help menu item dark red
        set(help_m, 'foregroundColor', [0.6, 0 0]);
    elseif isempty(eeglabUpdater.lastTimeChecked)
        fprintf('Could not check for the latest EEGLAB version (internet may be disconnected).\n');
    else
        fprintf('You are using the latest version of EEGLAB.\n');
    end;    
    
catch
    fprintf('Updater could not be initialized.\n');
end;


% REMOVED MENUS
	%uimenu( tools_m, 'Label', 'Automatic comp. reject',  'enable', 'off', 'CallBack', '[EEG LASTCOM] = pop_rejcomp(EEG); eegh(LASTCOM); if ~isempty(LASTCOM), eeg_store(CURRENTSET); end;');
	%uimenu( tools_m, 'Label', 'Reject (synthesis)' , 'Separator', 'on', 'CallBack', '[EEG LASTCOM] = pop_rejall(EEG); eegh(LASTCOM); if ~isempty(LASTCOM), eeg_store; end; eeglab(''redraw'');');

     function command_on_update_menu_click(callerHandle, tmp, eeglabUpdater, installDirectory, goOneFolderLevelIn, backGroundColor)
         postInstallCallbackString = 'clear all function functions; eeglab';
         eeglabUpdater.launchGui(installDirectory, goOneFolderLevelIn, backGroundColor, postInstallCallbackString);
    
%     
% --------------------
% draw the main figure
% --------------------

function tb = eeg_mainfig(onearg);

icadefs;
COLOR = BACKEEGLABCOLOR;
WINMINX         = 17;
WINMAXX         = 260;
WINYDEC			= 13;
NBLINES         = 16;
WINY		    = WINYDEC*NBLINES;
javaChatFlag    = 1;

BORDERINT       = 4;
BORDEREXT       = 10;
comp = computer;
if strcmpi(comp(1:3), 'GLN') % Linux
    FONTNAME        = 'courier';
    FONTSIZE        = 8;
elseif strcmpi(comp(1:3), 'MAC') 
    FONTNAME        = 'courier';
    FONTSIZE        = 8;
else
    FONTNAME        = '';
    FONTSIZE        = 11;
end;    

hh = findobj('tag', 'EEGLAB');
if ~isempty(hh)
    disp('EEGLAB warning: there can be only one EEGLAB window, closing old one');
    close(hh);
end;
if strcmpi(onearg, 'remote')
    figure(	'name', [ 'EEGLAB v' eeg_getversion ], ... 
	'numbertitle', 'off', ...
	'Position',[200 100 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) 30 ], ...
	'color', COLOR, ...
	'Tag','EEGLAB', ...
	'Userdata', {[] []});
	%'resize', 'off', ...
    return;
end;

W_MAIN = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'name', [ 'EEGLAB v' eeg_getversion ], ... 
	'numbertitle', 'off', ...
	'Position',[200 100 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) (WINY+2*BORDERINT+2*BORDEREXT) ], ...
	'color', COLOR, ...
	'Tag','EEGLAB', ...
    'visible', 'off', ...   
	'Userdata', {[] []});
%	'resize', 'off', ...

% java chat
eeglab_options;
if option_chat == 1
    if is_sccn
        disp('Starting chat...');
        tmpp = fileparts(which('startpane.m'));
        if isempty(tmpp) || ~ismatlab
            disp('Cannot start chat');
            tb = [];
        else
            disp(' ----------------------------------- ');
            disp('| EEGLAB chat 0.9                   |');
            disp('| The chat currently only works     |'); 
            disp('| at the University of CA San Diego |');
            disp(' ----------------------------------- ');

            javaaddpath(fullfile(tmpp, 'Chat_with_pane.jar'));
            eval('import client.EEGLABchat.*;');
            eval('import client.VisualToolbar;');
            eval('import java.awt.*;');
            eval('import javax.swing.*;');

            try
                tb = VisualToolbar('137.110.244.26');
                F = W_MAIN;
                tb.setPreferredSize(Dimension(0, 75));

                javacomponent(tb,'South',F);
                javaclose = ['userdat = get(gcbf, ''userdata'');' ...
                             ' tb = userdat{3};' ...
                             'clear userdat; delete(gcbf); tb.close; clear tb' ];
                set(gcf, 'CloseRequestFcn',javaclose);

                refresh(F);
            catch,
                tb = [];
            end;
        end;
    else
        tb = [];
    end;
else
    tb = [];
end;

try,
    set(W_MAIN, 'NextPlot','new');
catch, end;

if ismatlab
    BackgroundColor = get(gcf, 'color'); %[0.701960784313725 0.701960784313725 0.701960784313725];
    H_MAIN(1) = uicontrol('Parent',W_MAIN, ...
        'Units','points', ...
        'BackgroundColor',COLOR, ...
        'ListboxTop',0, ...
        'HorizontalAlignment', 'left',...
        'Position',[BORDEREXT   BORDEREXT  (WINMINX+WINMAXX+2*BORDERINT)  (WINY)], ...
        'Style','frame', ...
       'Tag','Frame1');

    geometry = { [1] [1] [1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1] };
    listui = { { 'style', 'text', 'string', 'Parameters of the current set', 'tag', 'win0' } { } ...
               { 'style', 'text', 'tag', 'win1', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win2', 'string', 'Channels per frame', 'userdata', 'datinfo'} ...
               { 'style', 'text', 'tag', 'val2', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win3', 'string', 'Frames per epoch', 'userdata', 'datinfo'} ...
               { 'style', 'text', 'tag', 'val3', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win4', 'string', 'Epochs', 'userdata', 'datinfo'} ...
               { 'style', 'text', 'tag', 'val4', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win5', 'string', 'Events', 'userdata', 'datinfo'} ...
               { 'style', 'text', 'tag', 'val5', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win6', 'string', 'Sampling rate (Hz)', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'val6', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win7', 'string', 'Epoch start (sec)', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'val7', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win8', 'string', 'Epoch end (sec)', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'val8', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win9', 'string', 'Average reference', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'val9', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win10', 'string', 'Channel locations', 'userdata', 'datinfo'} ...
               { 'style', 'text', 'tag', 'val10', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win11', 'string', 'ICA weights', 'userdata', 'datinfo'  } ...
               { 'style', 'text', 'tag', 'val11', 'string', ' ', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'win12', 'string', 'Dataset size (Mb)', 'userdata', 'datinfo' } ...
               { 'style', 'text', 'tag', 'val12', 'string', ' ', 'userdata', 'datinfo' } {} };
    supergui(gcf, geometry, [], listui{:});
    geometry = { [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] };
    listui = { { } ...
               { } ...
               { 'style', 'text', 'tag', 'mainwin1', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin2', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin3', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin4', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin5', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin6', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin7', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin8', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin9', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin10', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin11', 'string', ' ', 'userdata', 'fullline' } ...
               { 'style', 'text', 'tag', 'mainwin12', 'string', ' ', 'userdata', 'fullline' }  ...
               { 'style', 'text', 'tag', 'mainwin13', 'string', ' ', 'userdata', 'fullline' } {} };
    supergui(gcf, geometry, [], listui{:});

    titleh   = findobj('parent', gcf, 'tag', 'win0');
    alltexth = findobj('parent', gcf, 'style', 'text');
    alltexth = setdiff(alltexth, titleh);

    set(gcf, 'Position',[200 100 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) (WINY+2*BORDERINT+2*BORDEREXT) ]);
    set(titleh, 'fontsize', 14, 'fontweight', 'bold');
    set(alltexth, 'fontname', FONTNAME, 'fontsize', FONTSIZE);
    set(W_MAIN, 'visible', 'on');
end;

return;

% eeglab(''redraw'')() - Update EEGLAB menus based on values of global variables.
%
% Usage: >> eeglab(''redraw'')( );
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_global(), eeglab()

% WHEN THIS FUNCTION WAS SEPARATED
% Revision 1.21  2002/04/23 19:09:25  arno
% adding automatic dataset search
% Revision 1.20  2002/04/18 20:02:23  arno
% retrIeve
% Revision 1.18  2002/04/18 16:28:28  scott
% EEG.averef printed as 'Yes' or 'No' -sm
% Revision 1.16  2002/04/18 16:03:15  scott
% editted "Events/epoch info (nb) -> Events  -sm
% Revision 1.14  2002/04/18 14:46:58  scott
% editted main window help msg -sm
% Revision 1.10  2002/04/18 03:02:17  scott
% edited opening instructions -sm
% Revision 1.9  2002/04/11 18:23:33  arno
% Oups, typo which crashed EEGLAB
% Revision 1.8  2002/04/11 18:07:59  arno
% adding average reference variable
% Revision 1.7  2002/04/11 17:49:40  arno
% corrected operator precedence problem
% Revision 1.6  2002/04/11 15:36:55  scott
% added parentheses to final ( - & - ), line 84. ARNO PLEASE CHECK -sm
% Revision 1.5  2002/04/11 15:34:50  scott
% put isempty(CURRENTSET) first in line ~80 -sm
% Revision 1.4  2002/04/11 15:31:47  scott
% added test isempty(CURRENTSET) line 78 -sm
% Revision 1.3  2002/04/11 01:41:27  arno
% checking dataset ... and inteligent menu update
% Revision 1.2  2002/04/09 20:47:41  arno
% introducing event number into gui

function updatemenu();
eeg_global;

W_MAIN = findobj('tag', 'EEGLAB');
EEGUSERDAT = get(W_MAIN, 'userdata');
H_MAIN  = EEGUSERDAT{1};
EEGMENU = EEGUSERDAT{2};
if length(EEGUSERDAT) > 2
     tb = EEGUSERDAT{3};
else tb = [];
end;
if ~isempty(tb) && ~isstr(tb)
    tb.RefreshToolbar();
end;
if exist('CURRENTSET') ~= 1, CURRENTSET = 0; end;
if isempty(ALLEEG), ALLEEG = []; end;
if isempty(EEG), EEG = []; end;

% test if the menu is present  
try
	figure(W_MAIN);
	set_m   = findobj( 'parent', W_MAIN, 'Label', 'Datasets');
catch, return; end;
index = 1;
indexmenu = 1;
MAX_SET = max(length( ALLEEG ), length(EEGMENU)-1);
	
clear functions;
eeglab_options;
if isempty(ALLEEG) && ~isempty(EEG) && ~isempty(EEG.data)
    ALLEEG = EEG;
end;

% setting the dataset menu
% ------------------------
while( index <= MAX_SET)
    try
        set( EEGMENU(index), 'Label', '------', 'checked', 'off');
    catch,
        if mod(index, 30) == 0
            tag = [ 'More (' int2str(index/30) ') ->' ];
            tmp_m = findobj('label', tag);
            if isempty(tmp_m)
                 set_m = uimenu( set_m, 'Label', tag, 'userdata', 'study:on'); 
            else set_m = tmp_m;
            end;	
        end;
        try
            set( EEGMENU(index), 'Label', '------', 'checked', 'off');
        catch, EEGMENU(index) = uimenu( set_m, 'Label', '------', 'Enable', 'on'); end;	
    end;        
	set( EEGMENU(index), 'Enable', 'on', 'separator', 'off' );
	try, ALLEEG(index).data;
		if ~isempty( ALLEEG(index).data)
            
            cb_retrieve = [ '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''retrieve'', ' int2str(index) ', ''study'', ~isempty(STUDY)+0);' ...
                            'if CURRENTSTUDY & ~isempty(LASTCOM), CURRENTSTUDY = 0; LASTCOM = [ ''CURRENTSTUDY = 0;'' LASTCOM ]; end; eegh(LASTCOM);' ...
                            'eeglab(''redraw'');' ];
            
       		menutitle   = sprintf('Dataset %d:%s', index, ALLEEG(index).setname);
			set( EEGMENU(index), 'Label', menutitle, 'userdata', 'study:on');
			set( EEGMENU(index), 'CallBack', cb_retrieve );
			set( EEGMENU(index), 'Enable', 'on' );
            if any(index == CURRENTSET), set( EEGMENU(index), 'checked', 'on' ); end;
		end;
	catch, end;	
	index = index+1;
end;
hh = findobj( 'parent', set_m, 'Label', '------');
set(hh, 'Enable', 'off');

% menu for selecting several datasets
% -----------------------------------
if index ~= 0
    cb_select = [ 'nonempty = find(~cellfun(''isempty'', { ALLEEG.data } ));' ...                  
                  'tmpind = pop_chansel({ ALLEEG(nonempty).setname }, ''withindex'', nonempty);' ... 
                  'if ~isempty(tmpind),' ...
                  '    [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''retrieve'', nonempty(tmpind), ''study'', ~isempty(STUDY)+0);' ...
                  '    eegh(LASTCOM);' ...
                  '    eeglab(''redraw'');' ...
                  'end;' ...
                  'clear tmpind nonempty;' ];
    if MAX_SET == length(EEGMENU), EEGMENU(end+1) = uimenu( set_m, 'Label', '------', 'Enable', 'on'); end;
    
    set(EEGMENU(end), 'enable', 'on', 'Label', 'Select multiple datasets', ...
                      'callback', cb_select, 'separator', 'on', 'userdata', 'study:on');
end;

% STUDY consistency
% -----------------
exist_study = 0;
if exist('STUDY') & exist('CURRENTSTUDY')

    % if study present, check study consistency with loaded datasets
    % --------------------------------------------------------------
    if ~isempty(STUDY)
        if length(ALLEEG) > length(STUDY.datasetinfo) || ~isfield(ALLEEG, 'data') || any(cellfun('isempty', {ALLEEG.data}))
            if strcmpi(STUDY.saved, 'no')
                res = questdlg2( strvcat('The study is not compatible with the datasets present in memory', ...
                                         'It is self consistent but EEGLAB is not be able to process it.', ...
                                         'Do you wish to save the study as it is (EEGLAB will prompt you to', ...
                                         'enter a file name) or do you wish to remove it'), 'Study inconsistency', 'Save and remove', 'Remove', 'Remove' );
                if strcmpi(res, 'Remove')
                    STUDY = [];
                    CURRENTSTUDY = 0;
                else
                    pop_savestudy(STUDY, ALLEEG);
                    STUDY = [];
                    CURRENTSTUDY = 0;
                end;
            else
                warndlg2( strvcat('The study was not compatible any more with the datasets present in memory.', ...
                                  'Since it had not changed since last saved, it was simply removed from', ...
                                  'memory.') );
                STUDY = [];
                CURRENTSTUDY = 0;
            end;
        end;
    end;
    
    if ~isempty(STUDY)
        exist_study = 1;
    end;
end;

% menu for selecting STUDY set
% ----------------------------
if exist_study
    cb_select = [ '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''retrieve'', [STUDY.datasetinfo.index], ''study'', 1);' ...
                  'if ~isempty(LASTCOM), CURRENTSTUDY = 1; LASTCOM = [ LASTCOM ''CURRENTSTUDY = 1;'' ]; end;' ...
                  'eegh(LASTCOM);' ...
                  'eeglab(''redraw'');' ];
    tmp_m = findobj('label', 'Select the study set');
    delete(tmp_m); % in case it is not at the end
    tmp_m = uimenu( set_m, 'Label', 'Select the study set', 'Enable', 'on', 'userdata', 'study:on');
    set(tmp_m, 'enable', 'on', 'callback', cb_select, 'separator', 'on');        
else 
    delete( findobj('label', 'Select the study set') );
end;

EEGUSERDAT{2} = EEGMENU;
set(W_MAIN, 'userdata', EEGUSERDAT);

if (isempty(CURRENTSET) | length(ALLEEG) < CURRENTSET(1) | CURRENTSET(1) == 0 | isempty(ALLEEG(CURRENTSET(1)).data))
	CURRENTSET = 0;
	for index = 1:length(ALLEEG)
		if ~isempty(ALLEEG(index).data)
			CURRENTSET = index;
			break;
		end;
	end;
	if CURRENTSET ~= 0
		eegh([ '[EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,' int2str(CURRENTSET) ');' ])
		[EEG ALLEEG] = eeg_retrieve(ALLEEG, CURRENTSET);	
	else 
		EEG = eeg_emptyset;
	end;
end;

if (isempty(EEG) | isempty(EEG(1).data)) & CURRENTSET(1) ~= 0
	eegh([ '[EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,' int2str(CURRENTSET) ');' ])
	[EEG ALLEEG] = eeg_retrieve(ALLEEG, CURRENTSET);	
end;

% test if dataset has changed
% ---------------------------
if length(EEG) == 1
    if ~isempty(ALLEEG) & CURRENTSET~= 0 & ~isequal(EEG.data, ALLEEG(CURRENTSET).data) & ~isnan(EEG.data(1))
        % the above comparison does not work for ome structures
        %tmpanswer = questdlg2(strvcat('The current EEG dataset has changed. What should eeglab do with the changes?', ' '), ...
        %                      'Dataset change detected', ...
        %                      'Keep changes', 'Delete changes', 'New dataset', 'Make new dataset');
        disp('Warning: for some reason, the backup dataset in EEGLAB memory does not');
        disp('         match the current dataset. The dataset in memory has been overwritten');
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eegh('[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);');
        
        %if tmpanswer(1) == 'D' % delete changes
        %    [EEG ALLEEG] = eeg_retrieve(ALLEEG, CURRENTSET);	
        %    eegh('[EEG ALLEEG] = eeg_retrieve( ALLEEG, CURRENTSET);');
        %elseif tmpanswer(1) == 'K' % keep changes
        %    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        %    eegh('[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);');
        %else % make new dataset
        %    [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); 
        %    eegh(LASTCOM);
        %    MAX_SET = max(length( ALLEEG ), length(EEGMENU));
        %end;
    end;
end;

% print some information on the main figure
% ------------------------------------------
g = myguihandles(gcf);
if ~isfield(g, 'win0') % no display
    return;
end;

study_selected = 0;
if exist('STUDY') & exist('CURRENTSTUDY')
    if CURRENTSTUDY == 1, study_selected = 1; end;
end;

menustatus = {};
if study_selected
    menustatus = { menustatus{:} 'study' };
    
    hh = findobj('parent', gcf, 'userdata', 'fullline'); set(hh, 'visible', 'off');
    hh = findobj('parent', gcf, 'userdata', 'datinfo');  set(hh, 'visible', 'on');

    % head string
    % -----------
    set( g.win0, 'String', sprintf('STUDY set: %s', STUDY.name) );
    
    % dataset type
    % ------------
    datasettype = unique( [ EEG.trials ] );
    if datasettype(1) == 1 & length(datasettype) == 1, datasettype = 'continuous';
    elseif datasettype(1) == 1,                        datasettype = 'epoched and continuous';
    else                                               datasettype = 'epoched';
    end;
    
    % number of channels and channel locations
    % ----------------------------------------
    chanlen    = unique( [ EEG.nbchan ] );
    chanlenstr = vararg2str( mattocell(chanlen) );
    anyempty    = unique( cellfun( 'isempty', { EEG.chanlocs }) );
    if length(anyempty) == 2,   chanlocs = 'mixed, yes and no';
    elseif anyempty == 0,       chanlocs = 'yes';
    else                        chanlocs = 'no';
    end;

    % ica weights
    % -----------
    anyempty    = unique( cellfun( 'isempty', { EEG.icaweights }) );
    if length(anyempty) == 2,   studystatus = 'Missing ICA dec.';
    elseif anyempty == 0,       studystatus = 'Ready to precluster';
    else                        studystatus = 'Missing ICA dec.';
    end;

    % consistency & other parameters
    % ------------------------------
    [EEG epochconsist] = eeg_checkset(EEG, 'epochconsist');        % epoch consistency
    [EEG chanconsist ] = eeg_checkset(EEG, 'chanconsist');         % channel consistency
    [EEG icaconsist  ] = eeg_checkset(EEG, 'icaconsist');          % ICA consistency
    totevents = num2str(sum( cellfun( 'length', { EEG.event }) )); % total number of events
    totsize   = whos('STUDY', 'ALLEEG');                              % total size
    if isempty(STUDY.session),   sessionstr = ''; else sessionstr = vararg2str(STUDY.session); end;
    if isempty(STUDY.condition), condstr    = ''; else condstr    = vararg2str(STUDY.condition); end;
    
    % determine study status
    % ----------------------
    if isfield(STUDY.etc, 'preclust')
        if ~isempty( STUDY.etc.preclust )
            studystatus = 'Pre-clustered';
        elseif length(STUDY.cluster) > 1
            studystatus = 'Clustered';
        end;
    elseif length(STUDY.cluster) > 1
        studystatus = 'Clustered';
    end;        
    
    % text
    % ----
    set( g.win2, 'String', 'Study task name');
    set( g.win3, 'String', 'Nb of subjects');
    set( g.win4, 'String', 'Nb of conditions');
    set( g.win5, 'String', 'Nb of sessions');
    set( g.win6, 'String', 'Nb of groups');
    set( g.win7, 'String', 'Epoch consistency');
    set( g.win8, 'String', 'Channels per frame');
    set( g.win9, 'String', 'Channel locations');
    set( g.win10, 'String', 'Clusters');
    set( g.win11, 'String', 'Status');
    set( g.win12, 'String', 'Total size (Mb)');
    
    % values
    % ------
    fullfilename = fullfile( STUDY.filepath, STUDY.filename);
    if length(fullfilename) > 26
        set( g.win1, 'String', sprintf('Study filename: ...%s\n', fullfilename(max(1,length(fullfilename)-26):end) ));
    else
        set( g.win1, 'String', sprintf('Study filename: %s\n'   , fullfilename));
    end;        	
    condconsist  = std_checkconsist(STUDY, 'uniform', 'condition');
    groupconsist = std_checkconsist(STUDY, 'uniform', 'group');
    sessconsist  = std_checkconsist(STUDY, 'uniform', 'session');
    txtcond  = fastif(condconsist , ' per subject', ' (some missing)');
    txtgroup = fastif(groupconsist, ' per subject', ' (some missing)');
    txtsess  = fastif(sessconsist , ' per subject', ' (some missing)');
    set( g.val2, 'String', STUDY.task);
    set( g.val3, 'String', int2str(max(1, length(STUDY.subject))));
    set( g.val4, 'String', [ int2str(max(1, length(STUDY.condition))) txtcond ]);
    set( g.val5, 'String', [ int2str(max(1, length(STUDY.session)))   txtsess ]);
    set( g.val6, 'String', [ int2str(max(1, length(STUDY.group)))    txtgroup ]);
    set( g.val7, 'String', epochconsist);
    set( g.val8, 'String', chanlenstr);
    set( g.val9, 'String', chanlocs);
    set( g.val10, 'String', length(STUDY.cluster));
    set( g.val11, 'String', studystatus);
    set( g.val12, 'String', num2str(round(sum( [ totsize.bytes] )/1E6*10)/10));        
    
elseif (exist('EEG') == 1) & ~isnumeric(EEG) & ~isempty(EEG(1).data) 

    hh = findobj('parent', gcf, 'userdata', 'fullline'); set(hh, 'visible', 'off');
    hh = findobj('parent', gcf, 'userdata', 'datinfo');  set(hh, 'visible', 'on');
    
    if length(EEG) > 1 % several datasets

        menustatus = { menustatus{:} 'multiple_datasets' };
        
        % head string
        % -----------
        strsetnum = 'Datasets ';
        for i = CURRENTSET
            strsetnum = [ strsetnum int2str(i) ',' ];
        end;
        strsetnum = strsetnum(1:end-1);
        set( g.win0, 'String', strsetnum);
        
        % dataset type
        % ------------
        datasettype = unique( [ EEG.trials ] );
        if datasettype(1) == 1 & length(datasettype) == 1, datasettype = 'continuous';
        elseif datasettype(1) == 1,                        datasettype = 'epoched and continuous';
        else                                               datasettype = 'epoched';
        end;
        
        % number of channels and channel locations
        % ----------------------------------------
        chanlen    = unique( [ EEG.nbchan ] );
        chanlenstr = vararg2str( mattocell(chanlen) );
        anyempty    = unique( cellfun( 'isempty', { EEG.chanlocs }) );
        if length(anyempty) == 2,   chanlocs = 'mixed, yes and no';
        elseif anyempty == 0,       chanlocs = 'yes';
        else                        chanlocs = 'no';
        end;

        % ica weights
        % -----------
        anyempty    = unique( cellfun( 'isempty', { EEG.icaweights }) );
        if length(anyempty) == 2,   icaweights = 'mixed, yes and no';
        elseif anyempty == 0,       icaweights = 'yes';
        else                        icaweights = 'no';
        end;

        % consistency & other parameters
        % ------------------------------
        [EEG epochconsist] = eeg_checkset(EEG, 'epochconsist');        % epoch consistency
        [EEG chanconsist ] = eeg_checkset(EEG, 'chanconsist');         % channel consistency
        [EEG icaconsist  ] = eeg_checkset(EEG, 'icaconsist');          % ICA consistency
        totevents = num2str(sum( cellfun( 'length', { EEG.event }) )); % total number of events
        srate     = vararg2str( mattocell( unique( [ EEG.srate ] ) )); % sampling rate
        totsize   = whos('EEG');                                       % total size
                
        % text
        % ----
        set( g.win2, 'String', 'Number of datasets');
        set( g.win3, 'String', 'Dataset type');
        set( g.win4, 'String', 'Epoch consistency');
        set( g.win5, 'String', 'Channels per frame');
        set( g.win6, 'String', 'Channel consistency');
        set( g.win7, 'String', 'Channel locations');
        set( g.win8, 'String', 'Events (total)');
        set( g.win9, 'String', 'Sampling rate (Hz)');
        set( g.win10, 'String', 'ICA weights');
        set( g.win11, 'String', 'Identical ICA');
        set( g.win12, 'String', 'Total size (Mb)');

        % values
        % ------
        set( g.win1, 'String', sprintf('Groupname: -(soon)-\n'));
        set( g.val2, 'String', int2str(length(EEG)));
        set( g.val3, 'String', datasettype);
        set( g.val4, 'String', epochconsist);
        set( g.val5, 'String', chanlenstr);
        set( g.val6, 'String', chanconsist);
        set( g.val7, 'String', chanlocs);
        set( g.val8, 'String', totevents);
        set( g.val9, 'String', srate);
        set( g.val10, 'String', icaweights);
        set( g.val11, 'String', icaconsist);
        set( g.val12, 'String', num2str(round(totsize.bytes/1E6*10)/10));        
        
    else % one continous dataset selected
        
        menustatus = { menustatus{:} 'continuous_dataset' };
        
        % text
        % ----
        set( g.win2, 'String', 'Channels per frame');
        set( g.win3, 'String', 'Frames per epoch');
        set( g.win4, 'String', 'Epochs');
        set( g.win5, 'String', 'Events');
        set( g.win6, 'String', 'Sampling rate (Hz)');
        set( g.win7, 'String', 'Epoch start (sec)');
        set( g.win8, 'String', 'Epoch end (sec)');
        set( g.win9, 'String', 'Reference');
        set( g.win10, 'String', 'Channel locations');
        set( g.win11, 'String', 'ICA weights');
        set( g.win12, 'String', 'Dataset size (Mb)');
        
        if CURRENTSET == 0, strsetnum = '';
        else                strsetnum = ['#' int2str(CURRENTSET) ': '];
        end;
        maxchar = 28;
        if ~isempty( EEG.setname )
            if length(EEG.setname) > maxchar+2
                set( g.win0, 'String', [strsetnum EEG.setname(1:min(maxchar,length(EEG.setname))) '...' ]);
            else set( g.win0, 'String', [strsetnum EEG.setname ]);
            end;
        else
            set( g.win0, 'String', [strsetnum '(no dataset name)' ] );
        end;

        fullfilename = fullfile(EEG.filepath, EEG.filename);
        if ~isempty(fullfilename)
            if length(fullfilename) > 26
                set( g.win1, 'String', sprintf('Filename: ...%s\n', fullfilename(max(1,length(fullfilename)-26):end) ));
            else
                set( g.win1, 'String', sprintf('Filename: %s\n', fullfilename));
            end;        	
        else
            set( g.win1, 'String', sprintf('Filename: none\n'));
        end;
        
        set( g.val2, 'String', int2str(fastif(isempty(EEG.data), 0, size(EEG.data,1))));
        set( g.val3, 'String', int2str(EEG.pnts));
        set( g.val4, 'String', int2str(EEG.trials));
        set( g.val5, 'String', fastif(isempty(EEG.event), 'none', int2str(length(EEG.event))));
        set( g.val6, 'String', int2str( round(EEG.srate)) );
        if round(EEG.xmin) == EEG.xmin & round(EEG.xmax) == EEG.xmax
            set( g.val7, 'String', sprintf('%d\n', EEG.xmin));
            set( g.val8, 'String', sprintf('%d\n', EEG.xmax));
        else 
            set( g.val7, 'String', sprintf('%6.3f\n', EEG.xmin));
            set( g.val8, 'String', sprintf('%6.3f\n', EEG.xmax));
        end;

        % reference
        if isfield(EEG(1).chanlocs, 'ref')
            [curref tmp allinds] = unique( { EEG(1).chanlocs.ref });
            maxind = 1;
            for ind = unique(allinds)
                if length(find(allinds == ind)) > length(find(allinds == maxind))
                    maxind = ind;
                end;
            end;
            curref = curref{maxind};
            if isempty(curref), curref = 'unknown'; end;
        else curref = 'unknown';
        end;
        set( g.val9, 'String', curref);
        if isempty(EEG.chanlocs)
            set( g.val10, 'String', 'No');
        else
            if ~isfield(EEG.chanlocs, 'theta') | all(cellfun('isempty', { EEG.chanlocs.theta }))
                set( g.val10, 'String', 'No (labels only)');           
            else
                set( g.val10, 'String', 'Yes');
            end;
        end;
        
        set( g.val11, 'String', fastif(isempty(EEG.icasphere), 'No', 'Yes'));
        tmp = whos('EEG');
        if ~isa(EEG.data, 'memmapdata')
            set( g.val12, 'String', num2str(round(tmp.bytes/1E6*10)/10));
        else
            set( g.val12, 'String', [ num2str(round(tmp.bytes/1E6*10)/10) ' (memory mapped)' ]);
        end;

        if EEG.trials > 1
            menustatus = { menustatus{:} 'epoched_dataset' };
        else
            menustatus = { menustatus{:} 'continuous_dataset' };
        end
        if ~isfield(EEG.chanlocs, 'theta')
            menustatus = { menustatus{:} 'chanloc_absent' };
        end;
        if isempty(EEG.icaweights)
            menustatus = { menustatus{:} 'ica_absent' };
        end;
    end;
else
    menustatus = { menustatus{:} 'startup' };
    
	hh = findobj('parent', gcf, 'userdata', 'fullline'); set(hh, 'visible', 'on');
	hh = findobj('parent', gcf, 'userdata', 'datinfo');  set(hh, 'visible', 'off');
	set( g.win0, 'String', 'No current dataset');
	set( g.mainwin1, 'String', '- Create a new or load an existing dataset:');
	set( g.mainwin2, 'String', '   Use "File > Import data"           (new)'); 
	set( g.mainwin3, 'String', '   Or  "File > Load existing dataset" (old)');
	set( g.mainwin4, 'String', '- If new,');
	set( g.mainwin5, 'String', '  "File > Import epoch info" (data epochs) else');
	set( g.mainwin6, 'String', '  "File > Import event info" (continuous data)');
	set( g.mainwin7, 'String',  '  "Edit > Dataset info" (add/edit dataset info)');
	set( g.mainwin8, 'String', '  "File > Save dataset" (save dataset)');
	set( g.mainwin9, 'String', '- Prune data: "Edit > Select data"');
	set( g.mainwin10,'String', '- Reject data: "Tools > Reject continuous data"');
	set( g.mainwin11,'String', '- Epoch data: "Tools > Extract epochs"');
	set( g.mainwin12,'String', '- Remove baseline: "Tools > Remove baseline"');
	set( g.mainwin13,'String', '- Run ICA:    "Tools > Run ICA"');
end;

% enable selected menu items
% --------------------------
allmenus = findobj( W_MAIN, 'type', 'uimenu');
allstrs  = get(allmenus, 'userdata');
if any(strcmp(menustatus, 'startup'))
    
    set(allmenus, 'enable', 'on');  
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'startup:off'))), allstrs);  
    set(allmenus(indmatch), 'enable', 'off');
    
elseif any(strcmp(menustatus, 'study'))
    
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'study:on'))), allstrs);            
    set(allmenus          , 'enable', 'off');  
    set(allmenus(indmatch), 'enable', 'on');
    
elseif any(strcmp(menustatus, 'multiple_datasets'))
    
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'study:on'))), allstrs);            
    set(allmenus          , 'enable', 'off');  
    set(allmenus(indmatch), 'enable', 'on');        
    set(findobj('parent', W_MAIN, 'label', 'Study'), 'enable', 'off');
    
elseif any(strcmp(menustatus, 'epoched_dataset'))

    set(allmenus, 'enable', 'on');  
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'epoch:off'))), allstrs);  
    set(allmenus(indmatch), 'enable', 'off');
    
elseif any(strcmp(menustatus, 'continuous_dataset'))
    
    set(allmenus, 'enable', 'on');  
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'continuous:off'))), allstrs);  
    set(allmenus(indmatch), 'enable', 'off');
    
end;
if any(strcmp(menustatus, 'chanloc_absent'))
    
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'chanloc:on'))), allstrs);  
    set(allmenus(indmatch), 'enable', 'off');
    
end;
if any(strcmp(menustatus, 'ica_absent'))
    
    indmatch = cellfun(@(x)(~isempty(findstr(num2str(x), 'ica:on'))), allstrs);  
    set(allmenus(indmatch), 'enable', 'off');
    
end;

% adjust title extent
% -------------------
poswin0 = get(g.win0, 'position');
extwin0 = get(g.win0, 'extent');
set(g.win0, 'position', [poswin0(1:2) extwin0(3) extwin0(4)]);

return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;

function g = myguihandles(fig)
	g = [];
	hh = findobj('parent', gcf);
	for index = 1:length(hh)
		if ~isempty(get(hh(index), 'tag'))
			g = setfield(g, get(hh(index), 'tag'), hh(index));
		end;
	end;

% find a function path and add path if not present
% ------------------------------------------------
function myaddpath(eeglabpath, functionname, pathtoadd);

    tmpp = which(functionname);
    tmpnewpath = [ eeglabpath pathtoadd ];
    if ~isempty(tmpp)
        tmpp = tmpp(1:end-length(functionname));
        if length(tmpp) > length(tmpnewpath), tmpp = tmpp(1:end-1); end; % remove trailing filesep
        if length(tmpp) > length(tmpnewpath), tmpp = tmpp(1:end-1); end; % remove trailing filesep
        %disp([ tmpp '     |        ' tmpnewpath '(' num2str(~strcmpi(tmpnewpath, tmpp)) ')' ]);
        if ~strcmpi(tmpnewpath, tmpp)
            warning('off', 'MATLAB:dispatcher:nameConflict');
            addpath(tmpnewpath);
            warning('on', 'MATLAB:dispatcher:nameConflict');
        end;
    else
        %disp([ 'Adding new path ' tmpnewpath ]);
        addpath(tmpnewpath);
    end;

function addpathexist(p)
    if exist(p) == 7
        addpath(p);
    end;

function val = iseeglabdeployed2;
%val = 1; return;
if exist('isdeployed')
     val = isdeployed;
else val = 0;
end;

function buildhelpmenu;
