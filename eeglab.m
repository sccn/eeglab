% eeglab() - EEGLAB (Version 3.0) - Graphic user interface by Arnaud Delorme 
%            with Scott Makeig for the ICA EEG toolbox of Makeig et al. 
%            (now the EEGLAB toolbox) hosted 1997-2001 by the Computational 
%            Neuroscience Laboratory of Terrence J. Sejnowski at The Salk 
%            Institute, La Jolla. EEGLAB is hosted (2002-) by the Swartz Center 
%            for Computatational Neuroscience, Institute for Neural Computation,
%            University of California San Diego (http://sccn.ucsd.edu/)
%
% Description:
%   EEGLAB is Matlab software for processing continuous or epoched EEG data
%   using tools from the ICA/EEG toolbox. 
%   It is built to be useful either to novice or expert Matlab users. 
%   In normal use, the EEGLAB graphic interface calls graphic functions via 
%   pop-up windows. For the more advanced user, a single structure, 'EEG,' 
%   contains all dataset parameters and may be accessed and modified 
%   directly from the Matlab command line. Computations are carried out
%   using functions from the ICA/EEG toolbox (now the EEGLAB toolbox). 
%
% Note: To redaw the graphic interface without erasing datasets, simply
%       type >> eeglab redraw
%
% See first:
%    eeg_checkset(), for the structure of EEG datasets
%    The full documentation at http://www.cnl.salk.edu/~arno/eeglab.html
%    The licence under which the software is distributed (licence.txt)
%
% Main files:
% ---------- 
% eeglab()         - main graphic interface
% license.txt      - GNU license
% 
% Lastest signal processing functions (not in the latest distribution):
% --------------------------------------------------------------------
% cell2mat()      - cell to matrix, overwrite Neural network toolbox function
% compvar()       - compute component variance
% convolve()      - smart conv2 (fewer boundary problems)
% del2map()       - compute a surface Laplacian transform of the data
% eegplot()       - scrolling multichannel data viewer
% eegplot2event() - process data rejection info from eegplot()
% eegplot2trial() - process eegplot() rejection info
% eegrej()        - reject portions of continuous eeg data
% eegthresh()     - simple thresholding method
% entropy()       - compute component entropy
% epoch()         - extract epochs from a continuous dataset
% fastif()        - fast if function
% gabor2d()       - 2D Gabor matrix
% gauss2d()       - 2D Gauss matrix
% getallmenus()   - retrieve all menus of a GUI
% gradmap()       - compute the gradient of a map
% h()             - EEGLAB history function
% help2html()     - help header to HTML file conversion
% inputgui()      - function to program GUI (replace inputdlg)
% jointprob()     - joint probability function
% loadcnt()       - load continous CNT neuroscan file
% laplac2d()      - generate a Laplacian matrix output
% loadavg()       - load neuroscan .AVG file (not in EEGLAB, only for ERPs)
% loaddat()       - load neuroscan .DAT file
% loadeeg()       - load neuroscan .EEG file
% loadtxt()       - load text file
% makehtml()      - generate html pages for directories (uses help2html)
% mat2cell()      - matrix to cell, overwrite Matlab neural network toolbox function
% pophelp()       - format the help header  !!!
% realproba()     - compute observed probability (used by entropy)
% rejepoch()      - reject trials based on a given statistical measure
% rejkurt()       - calculate and reject data based on kurtosis
% rejtrend()      - reject EEG showing linear trends  !!!
% slider()        - graphic slider function
% supergui()      - allow generation of advanced GUI
% readlocs()      - read location files .loc, .sph, .xyz, .elp (uses readelp)
% parsetxt()      - parse a line of text for 
% readelp()       - read Polhemus .ELP file
% textgui()       - create a text window with sliders (for help text)
%
% GUI Functions calling eponymous functions:
% -----------------------------------------
% pop_eegfilt()   - bandpass filter data (eegfilt())
% pop_eegplot()   - scrolling multichannel data viewer (eegplot())
% pop_eegthresh() - simple thresholding method (eegthresh())
% pop_envtopo()   - plot ERP data and component contributions (envtopo())
% pop_epoch()     - extract epochs from a continuous dataset (epoch())
% pop_erpimage()  - plot single epochs as an image (erpimage())
% pop_jointprob() - reject epochs using joint probability (jointprob())
% pop_loaddat()   - load Neuroscan .DAT info file (loaddat())
% pop_loadcnt()   - load Neuroscan .CNT data (lndcnt())
% pop_loadeeg()   - load Neuroscan .EEG data (loadeeg())
% pop_plotdata()  - plot data epochs in rectangular array (plotdata())
% pop_rejkurt()   - compute data kurtosis (rejkurt())
% pop_rejtrend()  - reject EEG epochs showing linear trends  (rejtrend())
% pop_resample()  - change data sampling rate (resample())
% pop_rmbase()    - remove epoch baseline (rmbase())
% pop_runica()    - run infomax ICA decomposition (runica())
% pop_timef()     - event-related time-frequency (timef())
% pop_timtopo()   - plot ERP and scalp maps  (timtopo())
% pop_topoplot()  - plot scalp maps (topoplot())
% pop_snapread()  - read Snapmaster .SMA files (snapread())
% pop_crossf()    - event-realted cross-coherence (crossf())
% pop_spectopo()  - plot all channel spectra and scalp maps (spectopo())
% pop_plottopo()  - plot a data epoch in a topographic array (plottopo())
% pop_readedf()   - read .EDF EEG data format (readedf())
% pop_headplot()  - plot a 3-D data scalp map (headplot())
% pop_averef()    - convert data to average reference (averef())
%
% Other GUI functions:
% -------------------
% pop_comments()       - edit dataset comment ('about') text
% pop_compareerps()    - compare two dataset ERPs using plottopo()
% pop_compprop()       - plot component properties (erpimage(), spectra, map)
% pop_copyset()        - copy dataset
% pop_dispcomp()       - display component scalp maps with reject buttons
% pop_editeventfield() - edit event fields
% pop_editeventvals()  - edit event values
% pop_editset()        - edit dataset information
% pop_icathresh()      - choose rejection thresholds (in development)
% pop_importepoch()    - import epoch info ASCII file
% pop_importevent()    - import event info ASCII file
% pop_importpres()     - import Presentation info file
% pop_loadset()        - load dataset
% pop_loadwks()        - load workspace
% pop_mergeset()       - merge two datasets
% pop_rejepoch()       - reject pre-identified epochs in a EEG dataset
% pop_rejspec()        - reject based on spectrum (computes spectrum -% eegthresh)
% pop_saveh()          - save EEGLAB command history
% pop_saveset()        - save dataset
% pop_savewks()        - save workspace
% pop_select()         - select data (epochs, time points, channels ...)
% pop_selectevent()    - select events
% pop_subcomp()        - subtract components from data
%
% Non-GUI functions use for handling the EEG structure:
% ----------------------------------------------------
% eeg_checkset()       - check dataset parameter consistency
% eeg_consts()         - constants definition
% pop_delset()         - delete dataset
% pop_editoptions()    - edit the option file
% eeg_emptyset()       - empty dataset
% eeg_epochformat()    - convert epoch array to structure
% eeg_eventformat()    - convert event array to structure
% eeg_getepochevent()  - return event values for a subset of event types
% eeg_global()         - global variables
% eeg_multieegplot()   - plot several rejections (using different colors)
% eeg_options()        - option file
% eeg_rejsuperpose()   - use by rejmenu to superpose all rejections
% eeg_rejmacro()       - used by all rejection functions
% pop_rejmenu()        - rejection menu (with all rejection methods visible)
% eeg_retrieve()       - retrieve dataset from ALLSET
% eeg_store()          - store dataset into ALLSET
% eeg_updatemenu()     - update eeglab() gui
%
% Help functions:
% --------------
% eeg_helpadmin()      - help on admin function
% eeg_helphelp()       - help on help
% eeg_helpmenu()       - EEG help menus
% eeg_helppop()        - help on pop_ and eeg_ functions
% eeg_helpsigproc()    - help on signal processing functions

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

% $Log: not supported by cvs2svn $
% Revision 1.22  2002/04/21 00:03:43  scott
% [same] -sm
%
% Revision 1.21  2002/04/21 00:02:40  scott
% [same] -sm
%
% Revision 1.20  2002/04/21 00:01:29  scott
% 'datafile' -> 'data file' -sm
%
% Revision 1.19  2002/04/21 00:00:10  scott
% 'Import Matlab data array' -> 'Import data file or Matlab array' -sm
%
% Revision 1.18  2002/04/20 17:19:56  arno
% adding Done to callbacks
%
% Revision 1.17  2002/04/18 20:01:13  arno
% retrIeve
%
% Revision 1.16  2002/04/18 03:10:14  scott
% changed Edit menu item names -sm
%
% Revision 1.15  2002/04/18 02:36:22  scott
% Load existing dataset(s) -> Load existing dataset -sm
%
% Revision 1.14  2002/04/18 02:21:01  arno
% adding pop-up errors
%
% Revision 1.13  2002/04/18 00:22:42  scott
% Load dataset(s) -> Load existing dataset(s) -sm
%
% Revision 1.12  2002/04/11 19:09:09  arno
% typo in LASTCOM
%
% Revision 1.11  2002/04/11 19:08:02  arno
% adding lastcom argument for erpimage
%
% Revision 1.10  2002/04/11 18:04:08  arno
% not saving a new dataset when average referencing
%
% Revision 1.9  2002/04/11 03:35:31  arno
% editing load and save set menus
%
% Revision 1.8  2002/04/11 02:33:33  arno
% editing menus
%
% Revision 1.7  2002/04/10 23:16:49  arno
% chaning save worksapce menu
%
% Revision 1.6  2002/04/10 02:12:25  arno
% testing version control
%
% Revision 1.5  2002/04/10 01:01:10  arno
% testing version control
%
% Revision 1.4  2002/04/08 22:59:17  arno
% changing rmbase dataset saving status
%
% Revision 1.3  2002/04/06 02:54:10  arno
% change comments call
%
% Revision 1.2  2002/04/06 02:03:58  arno
% editing menus
%
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 
% 02-15-02 text interface editing -sm & ad 
% 03-08-02 removed toolbar option (matlab 5.2 compatibility) -ad
% 03-13-02 updated event function calls -ad
% 03-16-02 text interface editing -sm & ad 
% 3/19/02 Help msg edited by sm 

function eeglab( onearg )
eeg_options; % eeglaboptions;  <-- wasnt found 3/30 am -Scott
eeg_global;
eeg_consts;

if nargin < 1
	clear global EEG ALLEEG CURRENTSET ALLCOM LASTCOM;
	eeg_global;
	EEG = eeg_emptyset;
	evalin('base', 'eeg_global;');
	h('eeglab;');
else
	h('eeglab redraw;');
end;
  
% checking strings
% ----------------
e_try = 'try,';
e_catch = 'catch, errordlg(lasterr, ''EEGLAB error''); end;';
nocheck           = e_try;
check             = ['[EEG LASTCOM] = eeg_checkset(EEG, ''data''); h(LASTCOM);' e_try];
checkica          = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica''); h(LASTCOM);' e_try];
checkepoch        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch''); h(LASTCOM);' e_try];
checkevent        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''event''); h(LASTCOM);' e_try];
checkepochica     = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica''); h(LASTCOM);' e_try];
checkplot         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''chanloc''); h(LASTCOM);' e_try];
checkicaplot      = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'', ''chanloc''); h(LASTCOM);' e_try];
checkepochplot    = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''chanloc''); h(LASTCOM);' e_try];
checkepochicaplot = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'', ''chanloc''); h(LASTCOM);' e_try];
e_newnonempty     = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM), EEG = EEGTMP; clear EEGTMP; eeg_store; fprintf(''Done.''); end; eeg_updatemenu;'];
e_newset          = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM), eeg_store; fprintf(''Done.''); end; eeg_updatemenu;'];
e_store           = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM), eeg_store(CURRENTSET); fprintf(''Done.''); end; eeg_updatemenu;'];
e_hist            = [e_catch 'h(LASTCOM); eeg_updatemenu;'];

% menu definition
% --------------- 
eeg_mainfig;
W_MAIN = findobj('tag', 'EEGLAB');
EEGUSERDAT = get(W_MAIN, 'userdata');
set(W_MAIN, 'MenuBar', 'none');
first_m = uimenu( W_MAIN, 'Label', 'File');
	neuromenu = uimenu( first_m, 'Label', 'Import data'); 
	uimenu( neuromenu, 'Label', 'Read ascii/float data file or Matlab array'              ,     'CallBack', [ nocheck '[EEGTMP LASTCOM] = pop_importdata;' e_newnonempty ]);
	uimenu( neuromenu, 'Label', 'Read .SMA data file (Snapmaster)'       ,     'CallBack', [ nocheck 'eeg_global; [EEGTMP LASTCOM]= pop_snapread;' e_newnonempty ],  'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'Read .CNT data file (Neuroscan continuous)',  'CallBack', [ nocheck 'eeg_global; [EEGTMP LASTCOM]= pop_loadcnt;' e_newnonempty ], 'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'Read .EEG data file (Neuroscan epochs)'  ,    'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadeeg;' e_newnonempty ]); 
	importepoch = uimenu( first_m, 'Label', 'Import epoch info'); 
    uimenu( importepoch, 'Label', 'Import Matlab array or ASCII file',        'CallBack', [ check   '[EEG LASTCOM] = pop_importepoch(EEG);' e_store ]);
	uimenu( importepoch, 'Label', 'Import .DAT info file (Neuroscan epochs)', 'CallBack', [ check   'eeg_global; [EEG LASTCOM]= pop_loaddat(EEG);' e_store]); 
	importevent = uimenu( first_m, 'Label', 'Import event info'); 
	uimenu( importevent, 'Label', 'Import Matlab array or ASCII file',        'CallBack', [ check   '[EEG LASTCOM] = pop_importevent(EEG);' e_store]);
	uimenu( importevent, 'Label', 'Import .LOG event file (Presentation)'   , 'CallBack', [ check   'eeg_global; [EEG LASTCOM]= pop_importpres(EEG);' e_store]); 

	uimenu( first_m, 'Label', 'Load existing dataset' , 'Separator', 'on', 'CallBack', [ nocheck 'eeg_global; [TMPVAR LASTCOM]= pop_loadset; clear TMPVAR;' e_hist]); 
	uimenu( first_m, 'Label', 'Save current dataset'     , 'Separator', 'on', 'CallBack', [ check   'LASTCOM = pop_saveset(EEG);' e_hist]);
	uimenu( first_m, 'Label', 'Save datasets'                               , 'CallBack', [ check   'LASTCOM = pop_saveset;' e_hist ]);
	uimenu( first_m, 'Label', 'Clear dataset(s)'                            , 'CallBack', [ nocheck 'LASTCOM = pop_delset;' e_hist ]);
	uimenu( first_m, 'Label', 'Maximize memory'  , 'Separator', 'on'        , 'CallBack', [ nocheck 'LASTCOM = pop_editoptions;' e_hist]);
	uimenu( first_m, 'Label', 'Save history'     , 'Separator', 'on'        , 'CallBack', [ nocheck 'LASTCOM = pop_saveh(ALLCOM);' e_hist]);
	uimenu( first_m, 'Label', 'Quit'             , 'Separator', 'on'        , 'CallBack', ...
	       [ 'close(gcf); disp(''Do not forget to save the history, use >> pop_saveh(ALLCOM);''); clear global EEG ALLEEG LASTCOM CURRENTSET;']);

second_m = uimenu( W_MAIN, 'Label', 'Edit');
	uimenu( second_m, 'Label', 'Dataset info'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_editset(EEG);' e_store]);
	uimenu( second_m, 'Label', 'Event fields'     , 'CallBack', [ checkevent '[EEG LASTCOM] = pop_editeventfield(EEG);' e_store]);
	uimenu( second_m, 'Label', 'Event values'     , 'CallBack', [ checkevent '[EEG LASTCOM] = pop_editeventvals(EEG);' e_store]);
	uimenu( second_m, 'Label', 'About this dataset'    , 'CallBack', [ check      '[EEG.comments LASTCOM] =pop_comments(EEG.comments, ''About this dataset'');' e_store]);
	uimenu( second_m, 'Label', 'Select data'           , 'CallBack', [ check      '[EEG LASTCOM] = pop_select(EEG);' e_newset], 'Separator', 'on');
	uimenu( second_m, 'Label', 'Select events'         , 'CallBack', [ checkevent '[EEG TMP LASTCOM] = pop_selectevent(EEG); clear TMP' e_newset ]);
	uimenu( second_m, 'Label', 'Copy current dataset'  , 'CallBack', [ check      '[LASTCOM] = pop_copyset; h(LASTCOM); eeg_updatemenu;' e_hist], 'Separator', 'on');
	uimenu( second_m, 'Label', 'Append another dataset', 'CallBack', [ check      '[EEG LASTCOM] = pop_mergeset(EEG);' e_newset]);
	uimenu( second_m, 'Label', 'Delete dataset(s)'     , 'CallBack', [ nocheck    'LASTCOM = pop_delset;' e_hist]);
		
fourth_m  = uimenu( W_MAIN, 'Label', 'Tools');
	uimenu( fourth_m, 'Label', 'Change sampling rate', 'CallBack', [ check      '[EEG LASTCOM] = pop_resample(EEG);' e_newset]);
	uimenu( fourth_m, 'Label', 'Filter the data'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_eegfilt(EEG);' e_newset]);
	uimenu( fourth_m, 'Label', 'Average reference'   , 'CallBack', [ check      '[EEG LASTCOM] = pop_averef(EEG,1);' e_store]);
	uimenu( fourth_m, 'Label', 'Reject continuous data','CallBack',[ check      '[LASTCOM] = pop_eegplot(EEG, 1);' e_hist]);
	uimenu( fourth_m, 'Label', 'Extract epochs'      , 'CallBack', [ check      '[EEG tmp LASTCOM] = pop_epoch(EEG); clear tmp;' e_newset], 'Separator', 'on');
	uimenu( fourth_m, 'Label', 'Remove baseline'     , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_rmbase(EEG);' e_store]);
	fourth_sub1 = uimenu( fourth_m, 'Label', 'Reject data epochs');
	uimenu( fourth_m, 'Label', 'Run ICA'             , 'CallBack', [ check      '[EEG LASTCOM] = pop_runica(EEG);' e_store], 'foregroundcolor', 'b', 'Separator', 'on');
	uimenu( fourth_m, 'Label', 'Remove components'   , 'CallBack', [ checkica   '[EEG LASTCOM] = pop_subcomp(EEG);' e_newset]);
	fourth_sub2 = uimenu( fourth_m, 'Label', 'Reject using ICA');

	uimenu( fourth_sub1, 'Label', 'Reject data (all methods)', 'CallBack', [ check      'pop_rejmenu(EEG, 1); LASTCOM = '''';' e_hist]);
	uimenu( fourth_sub1, 'Label', 'Reject by inspection'     , 'CallBack', [ check      '[LASTCOM] = pop_eegplot(EEG, 1);' e_hist]);
	uimenu( fourth_sub1, 'Label', 'Reject extreme values'    , 'CallBack', [ checkepoch '[TMP LASTCOM] = pop_eegthresh(EEG, 1); clear TMP;' e_hist]);
	uimenu( fourth_sub1, 'Label', 'Reject flat line data'    , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_rejtrend(EEG, 1);' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject by probability'    , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_jointprob(EEG, 1);' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject by kurtosis'       , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_rejkurt(EEG, 1);' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject by spectra'        , 'CallBack', [ checkepoch '[EEG Itmp LASTCOM] = pop_rejspec(EEG, 1); clear Itmp;' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject labeled epochs', 'separator', 'on', 'CallBack', [ checkepoch ...
	     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 1,1,1,1,1,1,1,1); h(LASTCOM);' ...
	     '[EEG LASTCOM] = pop_rejepoch(EEG, EEG.reject.rejglobal,1);' e_newset]);
	   
	uimenu( fourth_sub2, 'Label', 'Reject components by map', 'CallBack', [ checkicaplot  '[EEG LASTCOM] = pop_selectcomps(EEG);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject data (all methods)','CallBack', [ checkepochica 'pop_rejmenu(EEG, 0); LASTCOM ='''';' e_hist], 'Separator', 'on');
	uimenu( fourth_sub2, 'Label', 'Reject by inspection',     'CallBack', [ checkica      '[LASTCOM] = pop_eegplot(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject extreme values',    'CallBack', [ checkepochica '[TMP LASTCOM] = pop_eegthresh(EEG, 0); clear TMP;' e_hist]);
	uimenu( fourth_sub2, 'Label', 'Reject flat line activity','CallBack', [ checkepochica '[EEG LASTCOM] = pop_rejtrend(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject by probability',    'CallBack', [ checkepochica '[EEG LASTCOM] = pop_jointprob(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject by kurtosis',       'CallBack', [ checkepochica '[EEG LASTCOM] = pop_rejkurt(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject by spectra',        'CallBack', [ checkepochica '[EEG LASTCOM] = pop_rejspec(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject labeled epochs', 'separator', 'on', 'CallBack', [ checkepochica ...
	     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); h(LASTCOM);' ...
	     '[EEG LASTCOM] = pop_rejepoch(EEG, EEG.reject.rejglobal,1); eeg_store; h(LASTCOM); eeg_updatemenu;' ]);
   
third_m = uimenu( W_MAIN, 'Label', 'Plot');
	loc_m = uimenu( third_m, 'Label', 'Channel locations'   );
       uimenu( loc_m, 'Label', 'By name'   , 'CallBack'  , [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''', ''''electrodes'''', ''''labelpoint'''');'']; eval(LASTCOM);' e_hist]);
	   uimenu( loc_m, 'Label', 'By number'   , 'CallBack', [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''', ''''electrodes'''', ''''numpoint'''');'']; eval(LASTCOM);' e_hist]);
    uimenu( third_m, 'Label', 'EEG data (scroll)'        , 'CallBack', [ check          'LASTCOM = pop_eegplot(EEG, 1, 0, 0);' e_hist], 'Separator', 'on');
	uimenu( third_m, 'Label', 'Channel spectra and maps' , 'CallBack', [ checkplot      'LASTCOM = pop_spectopo(EEG, 1);' e_hist]);
	uimenu( third_m, 'Label', 'Channel ERP image'        , 'CallBack', [ checkepoch     'LASTCOM = pop_erpimage(EEG, 1, LASTCOM);' e_hist]);
	ERP_m = uimenu( third_m, 'Label', 'ERP plots');
		uimenu( ERP_m, 'Label', 'ERP and scalp maps'     , 'CallBack', [ checkepochplot 'LASTCOM = pop_timtopo(EEG);' e_hist]);
		uimenu( ERP_m, 'Label', 'ERP in scalp array'     , 'CallBack', [ checkplot      'LASTCOM = pop_plottopo(EEG);' e_hist]);
		uimenu( ERP_m, 'Label', 'ERP in rect. array'     , 'CallBack', [ checkepoch     '[tmpeeg LASTCOM] = pop_plotdata(EEG, 1); clear tmpeeg;' e_hist]);
	topo_m = uimenu( third_m, 'Label', 'ERP maps');
		uimenu( topo_m, 'Label', 'As 2-D scalp maps'     , 'CallBack', [ checkplot      'LASTCOM = pop_topoplot(EEG, 1);' e_hist]);
		uimenu( topo_m, 'Label', 'As 3-D head plots'     , 'CallBack', [ checkplot      '[EEG LASTCOM] = pop_headplot(EEG, 1);' e_store]);
	uimenu( third_m, 'Label', 'Compare ERPs'             , 'CallBack', [ checkepoch     'LASTCOM = pop_compareerps(ALLEEG);' e_hist]);

    uimenu( third_m, 'Label', 'Component activations (scroll)', 'CallBack', [ checkica  '[LASTCOM] = pop_eegplot(EEG, 0, 0, 0);' e_hist],'Separator', 'on');
	topoica_m = uimenu( third_m, 'Label', 'Component maps');
		uimenu( topoica_m, 'Label', 'As 2-D scalp maps'  , 'CallBack', [ checkicaplot   'LASTCOM = pop_topoplot(EEG, 0);' e_hist]);
		uimenu( topoica_m, 'Label', 'As 3-D head plots'  , 'CallBack', [ checkicaplot   '[EEG LASTCOM] = pop_headplot(EEG, 0);' e_store]);
	uimenu( third_m, 'Label', 'Component properties'     , 'CallBack', [ checkicaplot   'LASTCOM = pop_compprop(EEG);' e_hist]);
	uimenu( third_m, 'Label', 'Component ERP image'      , 'CallBack', [ checkepochica  'LASTCOM = pop_erpimage(EEG, 0, LASTCOM);' e_hist]);
	ERPC_m = uimenu( third_m, 'Label', 'Component ERPs');
	   uimenu( ERPC_m, 'Label', 'On the same axis with maps', 'CallBack', [ checkepochicaplot 'LASTCOM = pop_envtopo(EEG);' e_hist]);
	   uimenu( ERPC_m, 'Label', 'In rectangular array'      , 'CallBack', [ checkepochica     '[tmpeeg LASTCOM] = pop_plotdata(EEG, 0); clear tmpeeg;' e_hist]);
	spec_m = uimenu( third_m, 'Label', 'Time-frequency', 'Separator', 'on');
		uimenu( spec_m, 'Label', 'Channel time-frequency'   , 'CallBack', [ checkepoch    'LASTCOM = pop_timef(EEG, 1);' e_hist]);
		uimenu( spec_m, 'Label', 'Channel cross-coherence'  , 'CallBack', [ checkepoch    'LASTCOM = pop_crossf(EEG, 1);' e_hist]);
		uimenu( spec_m, 'Label', 'Component time-frequency' , 'CallBack', [ checkepochica 'LASTCOM = pop_timef(EEG, 0);' e_hist],'Separator', 'on');
		uimenu( spec_m, 'Label', 'Component cross-coherence', 'CallBack', [ checkepochica 'LASTCOM = pop_crossf(EEG, 0);' e_hist]);
		
set_m   = uimenu( W_MAIN, 'Label', 'Datasets');
help_m  = uimenu( W_MAIN, 'Label', 'Help');
uimenu( help_m, 'Label', 'About EEGLAB', 'CallBack', 'pophelp(''eeglab'');');
uimenu( help_m, 'Label', 'About EEGLAB help', 'CallBack', 'pophelp(''eeg_helphelp'');');
uimenu( help_m, 'Label', 'EEGLAB menus', 'CallBack', 'eeg_helpmenu;');
help_subm1 = uimenu( help_m, 'Label', 'EEGLAB functions');
    uimenu( help_subm1, 'Label', 'Signal processing overview', 'CallBack', 'pophelp(''ica'');');
	uimenu( help_subm1, 'Label', 'Signal processing functions', 'callback', 'eeg_helpsigproc;');	
	uimenu( help_subm1, 'Label', 'Interactive (pop_) functions', 'callback', 'eeg_helppop;');	
help_subm2 = uimenu( help_m, 'Label', 'EEGLAB advanced');
    uimenu( help_subm2, 'Label', 'Dataset structure', 'CallBack', 'pophelp(''eeg_checkset'');');
	uimenu( help_subm2, 'Label', 'Admin functions', 'callback', 'eeg_helpadmin;');	

EEGMENU = uimenu( set_m, 'Label', '------', 'Enable', 'off');
set(W_MAIN, 'userdat', { EEGUSERDAT{1} EEGMENU });
eeg_updatemenu;

% REMOVED MENUS
	%uimenu( fourth_m, 'Label', 'Automatic comp. reject',  'enable', 'off', 'CallBack', '[EEG LASTCOM] = pop_rejcomp(EEG); h(LASTCOM); if ~isempty(LASTCOM), eeg_store(CURRENTSET); end;');
	%uimenu( fourth_m, 'Label', 'Reject (synthesis)' , 'Separator', 'on', 'CallBack', '[EEG LASTCOM] = pop_rejall(EEG); h(LASTCOM); if ~isempty(LASTCOM), eeg_store; end; eeg_updatemenu;');

% --------------------
% draw the main figure
% --------------------
function eeg_mainfig;

colormap('gray');
mat0 = colormap;
close;
WINMINX         = 17;
WINMAXX         = 260;

WINYDEC			= 13;
NBLINES         = 16;
WINY		    = WINYDEC*NBLINES;

BORDERINT       = 4;
BORDEREXT       = 10;
FONTNAME        = 'courrier';
FONTSIZE        = 11;

h = findobj('tag', 'EEGLAB');
if ~isempty(h)
    disp('EEGLAB warning: there can be only one EEGLAB window, closing old one');
    close(h);
end;

W_MAIN = figure('Units','points', ...
	'Colormap',mat0, ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'name', 'EEGLAB', ... 
	'numbertitle', 'off', ...
	'resize', 'off', ...
	'Position',[545.6028543307087 192.1136811023622 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) (WINY+2*BORDERINT+2*BORDEREXT) ], ...
	'Tag','EEGLAB', ...
	'Userdata', {[] []});
BackgroundColor = get(gcf, 'color'); %[0.701960784313725 0.701960784313725 0.701960784313725];
H_MAIN(1) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[BORDEREXT   BORDEREXT  (WINMINX+WINMAXX+2*BORDERINT)  (WINY)], ...
	'Style','frame', ...
   'Tag','Frame1');
H_MAIN(2) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'FontSize',14, ...
	'HorizontalAlignment', 'left',...
	'ListboxTop',0, ...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'String','Parameters of the current set', ...
	'Style','text', ...
   'Tag','1', ...
   'fontweight', 'bold', ...
   'fontsize', 12);
WINY = WINY-5/2*WINYDEC;
H_MAIN(3) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','2');
WINY = WINY-WINYDEC;
H_MAIN(4) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','3');
WINY = WINY-WINYDEC;
H_MAIN(5) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','4');
WINY = WINY-WINYDEC;
H_MAIN(6) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','5');
WINY = WINY-WINYDEC;
H_MAIN(7) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','6');
WINY = WINY-WINYDEC;
H_MAIN(8) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','7');
WINY = WINY-WINYDEC;
H_MAIN(9) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'HorizontalAlignment', 'left',...
	'ListboxTop',0, ...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','8');
WINY = WINY-WINYDEC;
H_MAIN(10) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','9');
WINY = WINY-WINYDEC;
H_MAIN(11) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'HorizontalAlignment', 'left',...
	'ListboxTop',0, ...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','10');
WINY = WINY-WINYDEC;
H_MAIN(12) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'HorizontalAlignment', 'left',...
	'ListboxTop',0, ...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','11');
WINY = WINY-WINYDEC;
H_MAIN(13) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','12');
WINY = WINY-WINYDEC;
H_MAIN(14) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','13');
WINY = WINY-WINYDEC;
H_MAIN(15) = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'ListboxTop',0, ...
	'HorizontalAlignment', 'left',...
	'Position',[WINMINX WINY WINMAXX WINYDEC], ...
	'Style','text', ...
	'Tag','14');

set(H_MAIN(3:end), 'fontname', FONTNAME, 'fontsize', FONTSIZE);
set(W_MAIN, 'userdata', { H_MAIN [] });
return;

h19 = uicontrol('Parent',W_MAIN, ...
	'Units','points', ...
	'BackgroundColor',BackgroundColor, ...
	'HorizontalAlignment', 'left',...
	'ListboxTop',0, ...
	'Position',[ BORDEREXT  (WINYDEC*NBLINES-3*WINYDEC/4)  (WINMINX+WINMAXX+2*BORDERINT) BORDERINT], ...
	'Style','frame', ...
	'Tag','Frame2');

