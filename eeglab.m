% eeglab() - Matlab graphic user interface environment for 
%   electrophysiological data analysis incorporating the ICA/EEG toolbox 
%   (Makeig et al.) developed at CNL / The Salk Institute, 1997-2001. 
%   Released 11/2002- as EEGLAB (Delorme, Makeig, et al.) at the Swartz Center 
%   for Computational Neuroscience, Institute for Neural Computation, 
%   University of California San Diego (http://sccn.ucsd.edu/). 
%   User feedback welcome: email eeglab@sccn.ucsd.edu
%
% Authors: Arnaud Delorme, Scott Makeig, with substantial contributions
%   from Colin Humphries, Sigurd Enghoff, Tzyy-Ping Jung, plus contributions 
%   from Tony Bell, Te-Won Lee, Luca Finelli and other contributors. 
%
% Description:
%   EEGLAB is Matlab software for processing continuous or event-related 
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
%
%   >> type "license.txt" % the GNU public license
%   >> web http://sccn.ucsd.edu/eeglab/tutorial/ % the EEGLAB tutorial
%   >> help eeg_checkset  % the EEG dataset structure
%
% Main files:
% ---------- 
% eeglab()        - main graphic interface
% license.txt     - GNU license
% 
% Functions added to EEGLAB: 
% --------------------------------------------------------------------
% cell2mat()      - cell to matrix, overwrites neural network toolbox function
% compvar()       - compute component variance
% convolve()      - smart conv2 (fewer boundary problems)
% del2map()       - compute a surface Laplacian transform of the data
% eegplot()       - scrolling multichannel data viewer (with data rejection)
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
% mat2cell()      - matrix to cell (local)
% pophelp()       - format the help header  !!!
% readedf()       - read binary EEG EDF file
% readegi()       - read binary EEG EGI file 
% read_erpss()    - read ERPSS data 
% readegihdr()    - read binary EEG EGI file header
% rejkurt()       - calculate and reject data based on kurtosis
% rejtrend()      - reject EEG showing linear trends  !!!
% reref()         - re-reference data
% slider()        - graphic slider function
% supergui()      - allow generation of advanced GUI
% readlocs()      - read location files .loc, .sph, .xyz, .elp (uses readelp)
% parsetxt()      - parse a line of text for 
% readelp()       - read Polhemus .ELP file
% textgui()       - create a text window with sliders (for help text)
%
% GUI Functions calling eponymous processing and plotting functions:
% ------------------------------------------------------------------
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
% pop_loadbva()   - load Brain Vision Analyser matlab files
% pop_plotdata()  - plot data epochs in rectangular array (plotdata())
% pop_read_erpss() - read ERPSS data (read_erpss())
% pop_readegi()   - load binary EGI data file (readegi())
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
% pop_reref()     - re-reference data (reref())
% pop_signalstat() - plot signal or component statistic (signalstat())
%
% Other GUI functions:
% -------------------
% pop_chanevent()      - import events stored in data channel(s)
% pop_comments()       - edit dataset comment ('about') text
% pop_compareerps()    - compare two dataset ERPs using plottopo()
% pop_prop()           - plot channel or component properties (erpimage, spectra, map)
% pop_copyset()        - copy dataset
% pop_dispcomp()       - display component scalp maps with reject buttons
% pop_editeventfield() - edit event fields
% pop_editeventvals()  - edit event values
% pop_editset()        - edit dataset information
% pop_export()         - export data or ica activity to ASCII file
% pop_expica()         - export ica weights or inverse matrix to ASCII file
% pop_icathresh()      - choose rejection thresholds (in development)
% pop_importepoch()    - import epoch info ASCII file
% pop_importevent()    - import event info ASCII file
% pop_importpres()     - import Presentation info file
% pop_loadset()        - load dataset
% pop_mergeset()       - merge two datasets
% pop_rejepoch()       - reject pre-identified epochs in a EEG dataset
% pop_rejspec()        - reject based on spectrum (computes spectrum -% eegthresh)
% pop_saveh()          - save EEGLAB command history
% pop_saveset()        - save dataset
% pop_select()         - select data (epochs, time points, channels ...)
% pop_selectevent()    - select events
% pop_subcomp()        - subtract components from data
%
% Non-GUI functions use for handling the EEG structure:
% ----------------------------------------------------
% eeg_checkset()       - check dataset parameter consistency
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
% eeg_retrieve()       - retrieve dataset from ALLEEG
% eeg_store()          - store dataset into ALLEEG
%
% Help functions:
% --------------
% eeg_helpadmin()      - help on admin function
% eeg_helphelp()       - help on help
% eeg_helpmenu()       - EEG help menus
% eeg_helppop()        - help on pop_ and eeg_ functions
% eeg_helpsigproc()    - help on signal processing functions

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.370  2005/07/28 17:38:55  arno
% same
%
% Revision 1.369  2005/07/28 17:36:53  arno
% enable/diable menus
%
% Revision 1.368  2005/07/27 19:27:32  arno
% allow selecting different dataset
% not different -> multiple
%
% Revision 1.367  2005/07/21 17:11:29  arno
% add check data
%
% Revision 1.366  2005/05/12 15:59:25  arno
% typo
%
% Revision 1.365  2005/05/12 15:57:24  arno
% changing for Matlab 7.2 compatibility thanks to Andreas Widmann
%
% Revision 1.364  2005/03/08 16:16:35  arno
% adding eeglab path
%
% Revision 1.363  2005/03/05 02:08:02  arno
% same
%
% Revision 1.362  2005/03/05 02:07:13  arno
% adding chaninfo
%
% Revision 1.361  2005/03/04 23:42:06  arno
% make function compatible with channel info
%
% Revision 1.360  2005/02/03 18:46:36  arno
% fixing link to license
%
% Revision 1.359  2005/02/02 20:25:37  arno
% menu
%
% Revision 1.358  2004/11/17 18:29:43  arno
% put BVA import in plugin
%
% Revision 1.357  2004/11/12 18:43:42  arno
% add other path for biosig
%
% Revision 1.356  2004/11/12 18:40:23  arno
% same
%
% Revision 1.355  2004/11/12 18:37:30  arno
% same
%
% Revision 1.354  2004/11/12 18:36:57  arno
% same
%
% Revision 1.353  2004/11/12 18:35:57  arno
% same
%
% Revision 1.352  2004/11/12 18:34:16  arno
% debug biosig menus ...
% /
%
% Revision 1.351  2004/11/12 18:27:17  arno
% adding biosig menus
%
% Revision 1.350  2004/11/12 18:09:59  arno
% debug last
%
% Revision 1.349  2004/11/12 18:09:20  arno
% use EEGLAB_VERSION instead of LATEST
%
% Revision 1.348  2004/11/05 17:59:43  arno
% put latest string for revision number (SCOTT: DO NOT CHANGE!)
%
% Revision 1.347  2004/10/07 18:37:51  hilit
% changed from && -> &
%
% Revision 1.346  2004/10/04 21:55:05  hilit
% plugin paths are now added to the end of the path
%
% Revision 1.345  2004/09/14 16:33:28  arno
% debug delete dataset
%
% Revision 1.344  2004/09/13 01:31:46  arno
% same
%
% Revision 1.343  2004/09/13 01:29:31  arno
% same
%
% Revision 1.342  2004/09/13 01:26:57  arno
% BIOSIG not in plugin folder
%
% Revision 1.341  2004/09/13 01:17:33  arno
% same
%
% Revision 1.340  2004/09/13 01:17:07  arno
% debug last
%
% Revision 1.339  2004/09/12 02:12:03  arno
% eeglab subpath for plugin
%
% Revision 1.338  2004/08/31 22:05:32  arno
% add comment
%
% Revision 1.337  2004/08/31 22:04:55  arno
% isequal compare data array
%
% Revision 1.336  2004/08/31 21:46:40  arno
% same
%
% Revision 1.335  2004/08/31 21:45:24  arno
% menu color
%
% Revision 1.334  2004/08/31 21:43:56  arno
% removing EDF/BDF menu
%
% Revision 1.333  2004/08/31 21:40:52  arno
% adding new command string
%
% Revision 1.332  2004/08/31 18:52:50  arno
% fix typo
%
% Revision 1.331  2004/08/31 18:00:07  arno
% move checks
%
% Revision 1.330  2004/08/31 16:30:46  scott
% dataset changes edits -sm
%
% Revision 1.329  2004/08/31 13:32:02  scott
% changed the 'dataset inconsistency' message -sm
% PS. Shouldnt the default be 'keep the changes'?
%
% Revision 1.328  2004/08/31 01:11:14  arno
% debug last
%
% Revision 1.327  2004/08/31 01:01:29  arno
% update MAX_SET
%
% Revision 1.326  2004/08/31 00:59:47  arno
% same
%
% Revision 1.325  2004/08/31 00:58:30  arno
% command line warning
%
% Revision 1.324  2004/08/25 17:51:21  arno
% debug last change
%
% Revision 1.323  2004/08/25 17:27:16  arno
% better history
%
% Revision 1.322  2004/07/10 00:09:51  arno
% same
%
% Revision 1.321  2004/07/10 00:09:09  arno
% same
%
% Revision 1.320  2004/07/10 00:07:55  arno
% same
%
% Revision 1.319  2004/07/10 00:04:16  arno
% debug msg
%
% Revision 1.318  2004/07/09 23:58:59  arno
% debug plugin
%
% Revision 1.317  2004/07/09 23:42:46  arno
% same
%
% Revision 1.316  2004/07/09 23:41:55  arno
% update plugin
%
% Revision 1.315  2004/07/02 18:17:23  arno
% saving a dataset is recorded in dataset history
%
% Revision 1.314  2004/06/01 21:25:54  arno
% history save
%
% Revision 1.313  2004/06/01 21:19:29  arno
% history for read binary
%
% Revision 1.312  2004/05/27 23:14:17  arno
% allowing to define EEG then call eeglab redraw
%
% Revision 1.311  2004/03/25 15:44:08  arno
% version option
%
% Revision 1.310  2004/03/13 03:04:21  arno
% history
%
% Revision 1.309  2004/03/04 18:48:28  arno
% debug last
%
% Revision 1.308  2004/03/04 18:46:39  arno
% update text for channel location
%
% Revision 1.307  2004/01/31 01:50:43  arno
% create new dataset after re-referencing
%
% Revision 1.306  2003/12/18 00:33:58  arno
% same
%
% Revision 1.305  2003/12/18 00:33:21  arno
% debug last
%
% Revision 1.304  2003/12/18 00:32:34  arno
% channel location detection
%
% Revision 1.303  2003/12/17 01:41:22  arno
% plugin color
%
% Revision 1.302  2003/12/16 23:31:36  arno
% adding brain vision analyser import function
%
% Revision 1.301  2003/12/12 01:17:55  arno
% debug history for channel locations
%
% Revision 1.300  2003/12/11 17:29:19  arno
% fixing lastcom for pop_chaneditset and pop_copyset
%
% Revision 1.299  2003/12/05 18:06:39  arno
% debug adding path
%
% Revision 1.298  2003/12/05 00:49:20  arno
% no local history when changing comments
%
% Revision 1.297  2003/12/05 00:17:03  arno
% history in datasets
%
% Revision 1.296  2003/12/04 03:08:27  arno
% plugin msg
%
% Revision 1.295  2003/12/04 02:05:59  arno
% resources
%
% Revision 1.294  2003/12/04 01:47:08  arno
% still adding paths
%
% Revision 1.293  2003/12/04 01:34:02  arno
% plugins
%
% Revision 1.292  2003/12/04 01:26:47  arno
% debug path
%
% Revision 1.291  2003/12/04 01:24:12  arno
% adding paths
%
% Revision 1.290  2003/12/04 01:14:17  arno
% adding paths
%
% Revision 1.289  2003/12/03 18:36:49  arno
% same
%
% Revision 1.288  2003/12/03 18:35:48  arno
% removing de-activated menus
%
% Revision 1.287  2003/12/03 02:48:35  arno
% enable spectra in all cases
%
% Revision 1.286  2003/12/03 00:20:11  arno
% removing version number once more (why?)
%
% Revision 1.285  2003/12/02 17:31:49  arno
% path debugging
%
% Revision 1.284  2003/12/02 17:27:36  arno
% still debuging path
%
% Revision 1.283  2003/12/02 17:25:09  arno
% debuging path add
%
% Revision 1.282  2003/12/01 00:44:06  arno
% undoing revision number
%
% Revision 1.281  2003/11/29 23:36:34  scott
% 4.0 -> 4.3
%
% Revision 1.280  2003/11/27 20:18:39  arno
% adding path datafile
%
% Revision 1.279  2003/11/27 20:10:11  arno
% remove trailing delimitor for path
%
% Revision 1.278  2003/11/27 20:05:42  arno
% debug path
%
% Revision 1.277  2003/11/27 20:00:57  arno
% still working on path
%
% Revision 1.276  2003/11/27 01:26:43  arno
% path
%
% Revision 1.275  2003/11/27 01:20:30  arno
% debug path
%
% Revision 1.274  2003/11/27 01:18:00  arno
% eeglab path
%
% Revision 1.273  2003/11/27 01:03:42  arno
% automatically adding paths
%
% Revision 1.272  2003/11/26 01:23:15  arno
% color for plugin in plot menu
%
% Revision 1.271  2003/11/25 22:58:23  arno
% same
%
% Revision 1.270  2003/11/25 22:57:21  arno
% pluginmenucolor
%
% Revision 1.269  2003/11/25 22:55:08  arno
% same
%
% Revision 1.268  2003/11/25 22:54:05  arno
% polugin menu color
%
% Revision 1.267  2003/11/25 19:40:58  arno
% removing BIOSIG
%
% Revision 1.266  2003/11/25 19:30:01  arno
% removing ant menus
%
% Revision 1.265  2003/11/25 19:10:11  arno
% smart plugin
%
% Revision 1.264  2003/11/25 19:00:18  arno
% same
%
% Revision 1.263  2003/11/25 18:59:31  arno
% typo last
%
% Revision 1.262  2003/11/25 18:58:38  arno
% adding ant file format
%
% Revision 1.261  2003/11/17 20:23:20  arno
% plot/compare menu
%
% Revision 1.260  2003/10/29 19:25:12  arno
% detecting ICALAB
%
% Revision 1.259  2003/10/29 19:06:26  arno
% text typo
%
% Revision 1.258  2003/10/29 19:02:59  arno
% clearing variable h
%
% Revision 1.257  2003/10/29 18:31:48  arno
% polishing message
%
% Revision 1.256  2003/10/29 18:29:11  arno
% same
%
% Revision 1.255  2003/10/29 18:28:35  arno
% debug last
%
% Revision 1.254  2003/10/29 18:27:33  arno
% adding biosign optional menu
%
% Revision 1.253  2003/10/29 00:50:57  arno
% check epoch when importing epoch info
%
% Revision 1.252  2003/08/11 00:55:17  arno
% remote option
%
% Revision 1.251  2003/07/31 22:39:22  arno
% debug last
%
% Revision 1.250  2003/07/31 22:32:35  arno
% same
%
% Revision 1.249  2003/07/31 22:30:23  arno
% debug storeallcall
%
% Revision 1.248  2003/07/30 02:14:10  arno
% envtopo compare error message
%
% Revision 1.247  2003/07/28 17:55:44  arno
% changing test for reference
%
% Revision 1.246  2003/07/28 15:21:16  arno
% remove averef
%
% Revision 1.245  2003/07/22 17:14:09  arno
% automatic reject button for channel scroll
%
% Revision 1.244  2003/06/05 15:28:09  arno
% update import menu
%
% Revision 1.243  2003/05/14 15:25:17  arno
% typo
%
% Revision 1.242  2003/05/14 15:21:45  arno
% disp done
% for save
%
% Revision 1.241  2003/05/14 15:15:50  arno
% adding export menu
%
% Revision 1.240  2003/05/12 16:23:01  scott
% edited text fprintfs
%
% Revision 1.239  2003/05/12 16:16:11  scott
% plugin "executing" -> plugin "adding"
%
% Revision 1.238  2003/05/12 15:58:15  arno
% updated store call
%
% Revision 1.237  2003/04/18 00:54:45  arno
% update menu label for compare ERPs
%
% Revision 1.236  2003/04/11 01:52:27  arno
% same
%
% Revision 1.235  2003/04/11 00:43:37  arno
% adding menu for reading segmented EGI files
%
% Revision 1.234  2003/03/18 00:31:19  arno
% debug last
%
% Revision 1.233  2003/03/18 00:30:26  arno
% adding menu for ICA erps
%
% Revision 1.232  2003/03/17 23:43:30  arno
% including new erp function
%
% Revision 1.231  2003/03/14 22:39:41  arno
% auto baseline removal after epoch extraction
%
% Revision 1.230  2003/03/14 03:20:34  arno
% adding menu to compare 2 envtopo
%
% Revision 1.229  2003/03/13 19:35:57  arno
% removing besa option
%
% Revision 1.228  2003/03/03 17:09:38  arno
% debug last
%
% Revision 1.227  2003/03/03 17:08:20  arno
% adding path detection for plugins
%
% Revision 1.226  2003/02/28 17:48:34  arno
% do not return output var for eeglab redraw
%
% Revision 1.225  2003/02/25 01:20:07  scott
% header edit -sm
%
% Revision 1.224  2003/02/24 19:56:40  arno
% debuginh last
%
% Revision 1.223  2003/02/24 19:54:53  arno
% do not crash eeglab if bad plugin
%
% Revision 1.222  2003/02/24 19:35:58  arno
% adding plugin functionality
%
% Revision 1.221  2003/02/21 22:25:52  arno
% small problem with 2 import menus
%
% Revision 1.220  2003/02/13 00:57:32  arno
% showing warning if low screen color depth
%
% Revision 1.219  2003/02/03 22:30:29  arno
% adding output command call to history
%
% Revision 1.218  2003/02/03 22:24:45  arno
% removing output if no output
%
% Revision 1.217  2003/02/03 17:15:02  arno
% adding ALLCOM output
%
% Revision 1.216  2003/02/03 17:09:59  arno
% adding outputs to EEGLAB
%
% Revision 1.215  2003/01/29 23:07:53  arno
% typo in readerpss in header
%
% Revision 1.214  2003/01/28 18:03:27  arno
% removing channel location file constraint for spectopo
%
% Revision 1.213  2003/01/24 01:28:27  arno
% adding the function to read ERPSS data
%
% Revision 1.212  2003/01/17 16:03:17  arno
% fixing window update problem for delete
%
% Revision 1.211  2003/01/10 01:15:28  arno
% do not redaw eeglab when only plotting
%
% Revision 1.210  2003/01/10 00:53:07  arno
% changing position once more
%
% Revision 1.209  2003/01/10 00:51:19  arno
% same
%
% Revision 1.208  2003/01/10 00:49:54  arno
% change default position
%
% Revision 1.207  2003/01/09 18:35:33  arno
% putting dataset name in title
%
% Revision 1.206  2003/01/03 22:09:15  arno
% nothing
%
% Revision 1.205  2002/11/18 02:39:43  arno
% allowing pop_rmbase to baseline portion of continuous data
%
% Revision 1.204  2002/11/15 18:33:07  arno
% updating besa checks
%
% Revision 1.203  2002/11/15 18:18:41  arno
% adding tutorial to the help menu
%
% Revision 1.202  2002/11/15 03:17:20  arno
% remove old functions from header
%
% Revision 1.201  2002/11/14 23:33:44  arno
% merging .EDF and .BDF menu lines
%
% Revision 1.200  2002/11/14 18:15:10  arno
% Average reference -> Re-reference
%
% Revision 1.199  2002/11/14 18:12:27  arno
% updating average reference flag
%
% Revision 1.198  2002/11/14 18:08:47  arno
% adding readegi
%
% Revision 1.197  2002/11/13 17:03:06  arno
% averef -> reref
%
% Revision 1.196  2002/11/13 02:40:35  arno
% nothing
%
% Revision 1.195  2002/11/12 22:57:16  arno
% adding check for extra channel
%
% Revision 1.194  2002/11/12 16:54:42  scott
% TIP text edit
%
% Revision 1.193  2002/11/12 02:15:49  scott
% Biosemi BDF -> standard BDF
%
% Revision 1.192  2002/11/12 02:12:34  arno
% add BDF function to file header
%
% Revision 1.191  2002/11/12 02:12:24  scott
% edit help msg
%
% Revision 1.190  2002/11/12 01:32:35  arno
% adding read bdf format
%
% Revision 1.189  2002/11/11 15:33:22  arno
% adding checks for BESA
%
% Revision 1.188  2002/11/09 18:33:08  scott
% BESA menu item edit
%
% Revision 1.187  2002/11/09 18:24:45  scott
% BESA menu item titles
%
% Revision 1.186  2002/11/09 17:53:28  scott
% edit dipole menu items
%
% Revision 1.185  2002/11/09 17:52:04  scott
% Edit dipole menu items
%
% Revision 1.184  2002/10/29 17:05:00  arno
% menu select event
%
% Revision 1.183  2002/10/28 02:25:16  arno
% debug h -> hh
%
% Revision 1.182  2002/10/25 23:45:17  arno
% debug last
%
% Revision 1.181  2002/10/25 23:40:05  arno
% h -> hh
%
% Revision 1.180  2002/10/23 22:28:26  arno
% releasing the timtopo ploting for continuous data
%
% Revision 1.179  2002/10/23 15:44:26  arno
% disable function warning
%
% Revision 1.178  2002/10/23 15:41:25  arno
% typo
%
% Revision 1.177  2002/10/23 15:33:56  arno
% put warning for chan edit
%
% Revision 1.176  2002/10/15 17:18:46  arno
% visible off and on
%
% Revision 1.175  2002/10/15 16:15:52  arno
% message
%
% Revision 1.174  2002/10/13 23:23:05  arno
% also clear EEGTMP if error
%
% Revision 1.173  2002/10/09 21:23:45  arno
% updating check for pop_envtopo
%
% Revision 1.172  2002/10/04 01:39:10  arno
% same
%
% Revision 1.171  2002/10/04 01:38:36  arno
% mri dipole summary menu
%
% Revision 1.170  2002/10/03 23:49:46  arno
% besaplot options
%
% Revision 1.169  2002/10/03 14:55:36  arno
% change besaplot summary
%
% Revision 1.168  2002/09/27 15:33:48  arno
% reverting old time range text displayed
%
% Revision 1.167  2002/09/23 23:32:34  arno
% updating display of time range
%
% Revision 1.166  2002/09/11 19:43:26  arno
% debug besa menu
%
% Revision 1.165  2002/09/08 00:45:03  scott
% pop_envtopo() menu label -scott
%
% Revision 1.164  2002/09/05 00:04:29  scott
% Plus -> With in menu -scott
%
% Revision 1.163  2002/09/05 00:03:06  scott
% Averaged ERP -> Channel ERPs in menu, for symmetry with comps -scott
%
% Revision 1.162  2002/09/04 23:20:04  arno
% swaping menu order
%
% Revision 1.161  2002/09/04 23:17:33  arno
% pop_compprop -> pop_prop
%
% Revision 1.160  2002/09/03 16:32:47  arno
% same, debug
%
% Revision 1.159  2002/09/03 16:31:55  arno
% dataset size
%
% Revision 1.158  2002/08/22 21:12:36  arno
% typo
%
% Revision 1.157  2002/08/22 17:18:21  arno
% more checks
%
% Revision 1.156  2002/08/22 17:13:48  arno
% debug eeglab besa
%
% Revision 1.155  2002/08/21 17:56:21  arno
% add menu to export rejections
%
% Revision 1.154  2002/08/21 02:22:00  arno
% add continuous data check for rejection continuous data menu
%
% Revision 1.153  2002/08/19 22:00:13  arno
% remove 10 as CR
%
% Revision 1.152  2002/08/19 19:50:23  arno
% hiding statistic menu if no ksstat func
%
% Revision 1.151  2002/08/17 22:23:23  scott
% As 2-D scalp series -> In 2-D
% As 3-D head plots -> In 3-D
%
% Revision 1.150  2002/08/17 21:59:13  scott
% As 2-D scalp maps -> In 2-D
% As 3-D head plots -> In 3-D
%
% Revision 1.149  2002/08/17 20:40:30  scott
% Plot > ERP maps -> Plot > ERP map series
%
% Revision 1.148  2002/08/17 20:38:49  scott
% Plot > ERP plots -> Plot > Averaged ERP
%
% Revision 1.147  2002/08/17 18:52:55  scott
% Reject epochs using ICA -> Reject data using ICA
%
% Revision 1.146  2002/08/17 18:51:06  scott
% Reject using ICA -> Reject epochs using ICA
%
% Revision 1.145  2002/08/17 02:08:21  arno
% debugging error catch
%
% Revision 1.144  2002/08/17 01:17:55  scott
% Plot > EEG data (scroll) -> Plot > Channel data (scroll)
%
% Revision 1.143  2002/08/17 00:44:15  arno
% remove rejepoch from header
%
% Revision 1.142  2002/08/17 00:35:17  arno
% editing header for internet
%
% Revision 1.141  2002/08/15 22:15:43  arno
% EDF
%
% Revision 1.140  2002/08/15 16:57:30  arno
% update text
%
% Revision 1.139  2002/08/15 16:38:00  arno
% debug stat menu
%
% Revision 1.138  2002/08/15 16:34:30  arno
% new statistic menu
%
% Revision 1.137  2002/08/15 15:48:59  arno
% debugging filename display
%
% Revision 1.136  2002/08/15 00:51:26  arno
% editing text
%
% Revision 1.135  2002/08/14 23:54:36  arno
% new format for dataset parameter ui
%
% Revision 1.134  2002/08/14 21:43:05  arno
% besa menu toggle
%
% Revision 1.133  2002/08/14 20:43:46  arno
% disable menus if function not available
%
% Revision 1.132  2002/08/14 00:15:32  arno
% debug unique dataset
%
% Revision 1.131  2002/08/14 00:03:13  arno
% *** empty log message ***
%
% Revision 1.130  2002/08/13 23:51:33  arno
% debug CURRENSET (T missing)
%
% Revision 1.129  2002/08/13 23:17:54  arno
% error message update
% .,
%
% Revision 1.128  2002/08/13 22:37:38  arno
% Dataset hide in low mem option
%
% Revision 1.127  2002/08/13 18:12:57  scott
% allow longer datset name to print on main window
%
% Revision 1.126  2002/08/13 18:06:54  scott
% title
%
% Revision 1.125  2002/08/13 18:06:21  scott
% title
%
% Revision 1.124  2002/08/13 18:02:49  scott
% help message and window title
%
% Revision 1.123  2002/08/13 17:12:15  scott
% menu
%
% Revision 1.122  2002/08/13 17:10:37  scott
% menu
%
% Revision 1.121  2002/08/13 17:09:10  scott
% menu
%
% Revision 1.120  2002/08/13 17:05:45  scott
% menu
%
% Revision 1.119  2002/08/13 17:04:13  scott
% menu
%
% Revision 1.118  2002/08/13 16:22:32  scott
% menu
%
% Revision 1.117  2002/08/13 00:30:14  scott
% text
%
% Revision 1.116  2002/08/13 00:14:46  scott
% menu text
%
% Revision 1.115  2002/08/13 00:13:37  scott
% menu text
% ,
% ,
%
% Revision 1.114  2002/08/13 00:09:36  scott
% menu names
%
% Revision 1.113  2002/08/12 22:31:46  arno
% menu text
%
% Revision 1.112  2002/08/12 22:04:35  arno
% text
%
% Revision 1.111  2002/08/12 21:49:20  arno
% menus
%
% Revision 1.110  2002/08/12 18:53:51  arno
% errordlg2
%
% Revision 1.109  2002/08/12 18:37:13  arno
% questdlg2
%
% Revision 1.108  2002/08/12 16:13:24  arno
% same
%
% Revision 1.107  2002/08/12 16:09:54  arno
% debug
%
% Revision 1.106  2002/08/12 01:15:40  arno
% update
%
% Revision 1.105  2002/08/11 22:50:38  arno
% change background
%
% Revision 1.104  2002/08/11 19:17:38  arno
% removing eeg_const eeg_updatemenu
%
% Revision 1.103  2002/08/11 17:37:47  arno
% same
%
% Revision 1.102  2002/08/11 17:37:23  arno
% header
%
% Revision 1.101  2002/08/09 17:59:38  arno
% adding signal statistics
%
% Revision 1.100  2002/08/01 00:35:15  arno
% debugging history for erpimage
%
% Revision 1.99  2002/07/31 16:13:50  arno
% updating check for spectopo
%
% Revision 1.98  2002/07/30 18:06:56  arno
% new spectopo menu
%
% Revision 1.97  2002/07/30 00:43:08  arno
% adding history for pop_editchan
%
% Revision 1.96  2002/07/30 00:08:00  arno
% adding history for crossf
%
% Revision 1.95  2002/07/29 23:26:47  arno
% adding history to timef
%
% Revision 1.94  2002/07/29 18:30:38  arno
% typo
%
% Revision 1.93  2002/07/29 17:55:37  arno
% add import pop_chanevent
%
% Revision 1.92  2002/07/29 16:07:12  arno
% adding help menus
%
% Revision 1.91  2002/07/29 15:29:55  arno
% updating message
%
% Revision 1.90  2002/07/27 00:46:31  arno
% debugging erp_image last modif
%
% Revision 1.89  2002/07/26 17:13:43  arno
% debugging
%
% Revision 1.88  2002/07/25 17:41:39  arno
% debugging
%
% Revision 1.87  2002/07/25 17:27:03  arno
% check all datasets after editing options
%
% Revision 1.86  2002/07/25 14:33:34  arno
% pop_editchan -> redraw eeglab
%
% Revision 1.85  2002/07/25 00:53:22  arno
% change default precision display
%
% Revision 1.84  2002/07/24 23:33:20  scott
% editing "On same axis (with maps)" -sm
%
% Revision 1.83  2002/07/24 16:47:29  arno
% debugging
%
% Revision 1.82  2002/07/24 01:17:58  arno
% edf debuging
%
% Revision 1.81  2002/07/24 01:16:48  arno
% adding pop_readedf
%
% Revision 1.80  2002/07/22 18:12:33  arno
% debugging eeg_retreive
%
% Revision 1.79  2002/07/18 18:06:20  arno
% same
%
% Revision 1.78  2002/07/18 18:05:17  arno
% same
%
% Revision 1.77  2002/07/18 18:04:15  arno
% same
%
% Revision 1.76  2002/07/18 18:02:55  arno
% correct dipoles
%
% Revision 1.75  2002/07/18 18:01:54  arno
% typo
%
% Revision 1.74  2002/07/18 17:58:01  arno
% add dipole menu
%
% Revision 1.73  2002/07/18 02:25:51  arno
% debug ALLSET->ALLEEG
%
% Revision 1.72  2002/07/13 00:12:26  arno
% debugging close
%
% Revision 1.71  2002/07/12 23:59:42  arno
% removing close satatement when creating a new eeglab
%
% Revision 1.70  2002/07/08 19:23:31  arno
% adding the import BCI menu command
%
% Revision 1.69  2002/06/25 14:36:31  arno
% debugging no dataset error
%
% Revision 1.68  2002/05/03 00:10:42  arno
% special case for channel editing
%
% Revision 1.67  2002/05/02 23:09:04  arno
% channel menu
%
% Revision 1.66  2002/05/02 23:07:53  arno
% eeglab redraw spetial case
%
% Revision 1.65  2002/05/01 03:21:10  arno
% adding channel editor
%
% Revision 1.64  2002/04/30 19:07:44  scott
% *** empty log message ***
%
% Revision 1.63  2002/04/30 19:07:10  scott
% *** empty log message ***
%
% Revision 1.62  2002/04/30 19:06:36  scott
% *** empty log message ***
%
% Revision 1.61  2002/04/30 19:04:21  scott
% *** empty log message ***
%
% Revision 1.60  2002/04/30 19:03:27  scott
% trying global COLOR -sm
%
% Revision 1.59  2002/04/30 19:01:21  scott
% *** empty log message ***
%
% Revision 1.58  2002/04/30 18:59:59  scott
% *** empty log message ***
%
% Revision 1.57  2002/04/30 18:56:33  scott
% *** empty log message ***
%
% Revision 1.56  2002/04/30 18:52:04  scott
% *** empty log message ***
%
% Revision 1.55  2002/04/30 18:41:05  scott
% editting help msg -sm
%
% Revision 1.54  2002/04/30 18:34:30  scott
% *** empty log message ***
%
% Revision 1.53  2002/04/30 18:33:50  scott
% same
%
% Revision 1.52  2002/04/30 18:33:17  scott
% same
%
% Revision 1.51  2002/04/30 18:32:40  scott
% same -sm
%
% Revision 1.50  2002/04/30 18:31:54  scott
% same -sm
% .,
%
% Revision 1.49  2002/04/30 18:30:38  scott
% same
%
% Revision 1.48  2002/04/30 18:30:07  scott
% editting color -sm
%
% Revision 1.47  2002/04/30 18:29:20  scott
% editting colors -sm
%
% Revision 1.46  2002/04/30 18:07:49  arno
% multiple colors
%
% Revision 1.45  2002/04/30 15:16:31  scott
% edited help msg -sm
%
% Revision 1.44  2002/04/26 21:18:20  arno
% updating pop_copyset call
%
% Revision 1.43  2002/04/26 20:46:50  arno
% removing old pop_erpimage
%
% Revision 1.42  2002/04/26 02:58:07  arno
% adding pop_newset.m
%
% Revision 1.41  2002/04/26 02:48:32  scott
% editing message
%
% Revision 1.40  2002/04/25 19:10:26  arno
% debugged typp
%
% Revision 1.39  2002/04/25 19:03:39  arno
% single dataset option
%
% Revision 1.38  2002/04/25 18:51:15  arno
% adding extra checks
%
% Revision 1.37  2002/04/25 18:50:02  arno
% updating memory options
%
% Revision 1.36  2002/04/24 17:05:46  arno
% CURRENTSET existance
%
% Revision 1.35  2002/04/24 15:12:03  scott
% [same] -sm
%
% Revision 1.34  2002/04/24 15:09:41  scott
% [same] -sm
%
% Revision 1.33  2002/04/24 15:09:07  scott
% [same] -sm
%
% Revision 1.32  2002/04/24 15:07:20  scott
% trying background color change -sm
%
% Revision 1.31  2002/04/23 23:49:09  arno
% new embeded update
%
% Revision 1.30  2002/04/23 21:34:20  arno
% updating pop_savesets.m
%
% Revision 1.29  2002/04/23 20:56:39  arno
% full reprogramming of the functionfurther updates of pop_delset call
%
% Revision 1.28  2002/04/23 19:08:55  arno
% debuging for pop_delset standalone call
%
% Revision 1.27  2002/04/23 17:58:16  arno
% modifying pop_loadset
%
% Revision 1.26  2002/04/23 02:03:17  arno
% updating menus for old and new pop_erpimage.m
%
% Revision 1.25  2002/04/22 23:34:40  arno
% correcting typo
%
% Revision 1.24  2002/04/22 01:00:51  arno
% fprintf->disp
%
% Revision 1.23  2002/04/21 01:00:14  scott
% edited help msg -sm
%
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

function [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab( onearg )

% add the paths
% -------------
eeglabpath = which('eeglab.m');
eeglabpath = eeglabpath(1:end-length('eeglab.m'));

% test for local SCCN copy
% ------------------------
addpath(eeglabpath);
if (strcmpi(computer, 'GLNX86') & exist( [ eeglabpath 'functions/adminfunc' ] ) == 7)
    myaddpath( eeglabpath, 'readeetraklocs.m', 'functions/sigprocfunc');
    myaddpath( eeglabpath, 'eeg_checkset.m',   'functions/adminfunc');
    myaddpath( eeglabpath, 'pop_loadbci.m',    'functions/popfunc');
    myaddpath( eeglabpath, 'icademo.m',        'functions/miscfunc');
    myaddpath( eeglabpath, 'VolumeMNI.bin',    'functions/resources');
elseif (strcmpi(computer, 'pcwin') & exist( [ eeglabpath 'functions\adminfunc' ] ) == 7)
    myaddpath( eeglabpath, 'readeetraklocs.m', 'functions\sigprocfunc');
    myaddpath( eeglabpath, 'eeg_checkset.m',   'functions\adminfunc');
    myaddpath( eeglabpath, 'pop_loadbci.m',    'functions\popfunc');
    myaddpath( eeglabpath, 'icademo.m',        'functions\miscfunc');
    myaddpath( eeglabpath, 'eeglab1020.ced',   'functions\resources');
else
    myaddpath( eeglabpath, 'readeetraklocs.m', 'functions');    
    funcpath = which('readeetraklocs.m');
    funcpath = funcpath(1:end-length('readeetraklocs.m'));
    myaddpath( funcpath , 'eeglab1020.ced', 'resources');    
end;
myaddpath( eeglabpath, 'eegplugin_dipfit', 'plugins');

eeg_options; 
if nargin == 1 &  strcmp(onearg, 'redraw')
    if evalin('base', 'exist(''EEG'')', '0') == 1
        evalin('base', 'warning off; eeg_global; warning on;');
    end;
end;
eeg_global;

% for the history function
% ------------------------
evalin('base'  , 'if exist(''h'') == 1, clear h; disp(''EEGLAB note: variable h cleared''); end;');

comtmp = 'warning off MATLAB:mir_warning_variable_used_as_function';
evalin('base'  , comtmp, '');
evalin('caller', comtmp, '');
    
evalin('base', 'eeg_global;');
if nargin < 1 | exist('EEG') ~= 1
	clear global EEG ALLEEG CURRENTSET ALLCOM LASTCOM;
	eeg_global;
	EEG = eeg_emptyset;
	h('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;');
    if get(0, 'screendepth') <= 8
        disp('Warning: screen color depth too low, some colors will be inaccurate in time-frequency plots');
    end;
end;

if nargin == 1
	if strcmp(onearg, 'verions')
        disp( [ 'EEGLAB v' EEGLAB_VERSION ] );
	elseif strcmp(onearg, 'redraw')
		W_MAIN = findobj('tag', 'EEGLAB');
		if ~isempty(W_MAIN)
			updatemenu;
            if nargout < 1, clear ALLEEG; end; % do not return output var
			return;
		else
			h('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab(''rebuild'');');
		end;
	elseif strcmp(onearg, 'besa');
		disp('Besa option deprecated. Download the BESA plugin to add the BESA menu.');
        h('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;');
	else
        h('[ALLEEG EEG CURRENTSET ALLCOM] = eeglab(''rebuild'');');
	end;
else 
    onearg = 'rebuild';
end;
ALLCOM = ALLCOM;
colordef white

% checking strings
% ----------------
e_try = 'try, if exist(''h'') == 1, clear h; disp(''EEGLAB note: variable h cleared''); end;';
e_catch = 'catch, errordlg2(lasterr, ''EEGLAB error''); LASTCOM= ''''; clear EEGTMP; end;';
nocheck           = e_try;
check             = ['[EEG LASTCOM] = eeg_checkset(EEG, ''data''); h(LASTCOM);' e_try];
checkcont         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''contdata''); h(LASTCOM);' e_try];
checkica          = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica''); h(LASTCOM);' e_try];
checkepoch        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch''); h(LASTCOM);' e_try];
checkevent        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''event''); h(LASTCOM);' e_try];
checkbesa         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''besa''); h(''% no history yet for BESA dipole localization'');' e_try];
checkepochica     = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica''); h(LASTCOM);' e_try];
checkplot         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''chanloc''); h(LASTCOM);' e_try];
checkicaplot      = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'', ''chanloc''); h(LASTCOM);' e_try];
checkepochplot    = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''chanloc''); h(LASTCOM);' e_try];
checkepochicaplot = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'', ''chanloc''); h(LASTCOM);' e_try];

storecall    = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); h(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);'');';
storeload    = '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); h(''[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);'');';
storenewcall = '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); h(LASTCOM);';
storeallcall = 'if ~isempty(ALLEEG) & ~isempty(ALLEEG(1).data), ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_checkset(EEG); h(''ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_checkset(EEG);''); end;';

e_newnonempty_nh     = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEGTMP), EEG = EEGTMP;' storenewcall 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_load_nh            = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEGTMP), EEG = EEGTMP;' storeload 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_newset_nh          = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG),' storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store_nh           = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG) & ~isempty(findstr(''='',LASTCOM)),' storecall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_storeall_nh        = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG),' storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist_nh            = [e_catch 'h(LASTCOM);'];
e_histdone_nh        = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM), disp(''Done.''); end;' ];

% same as above but also save history in dataset
% ----------------------------------------------
e_newnonempty   = [e_catch 'EEGTMP = h(LASTCOM, EEGTMP); if ~isempty(LASTCOM) & ~isempty(EEGTMP), EEG = EEGTMP;' storenewcall 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_load          = [e_catch 'EEGTMP = h(LASTCOM, EEGTMP); if ~isempty(LASTCOM) & ~isempty(EEGTMP), EEG = EEGTMP;' storeload 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_newset        = [e_catch 'EEG = h(LASTCOM, EEG); if ~isempty(LASTCOM) & ~isempty(EEG),' storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store         = [e_catch 'EEG = h(LASTCOM, EEG); if ~isempty(LASTCOM) & ~isempty(EEG) & ~isempty(findstr(''='',LASTCOM)),' storecall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_storeall      = [e_catch 'EEG = h(LASTCOM, EEG); if ~isempty(LASTCOM) & ~isempty(EEG),' storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist          = [e_catch 'EEG = h(LASTCOM, EEG);'];
e_histdone      = [e_catch 'EEG = h(LASTCOM, EEG); if ~isempty(LASTCOM), disp(''Done.''); end;' ];

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
catchstrs.new_non_empty          = e_newnonempty;

    % create eeglab figure
    % --------------------
    eeg_mainfig(onearg);

    % detecting icalab
    % ----------------
    if exist('icalab')
        disp('ICALAB toolbox detected (algo. added to "run ICA" interface)');
    end;
    
    % BIOSIG plugin (not in plugin folder)
    % ------------------------------------
    p = which('eeglab.m');
    p = p(1:findstr(p,'eeglab.m')-1);
    delimiter = p(end); if strcmpi(delimiter, ':'), delmiter = '::'; end;
    path_biosig = [ p '..' delimiter 'biosig' delimiter 't200' ];
    biosigflag = 0;
    if exist(path_biosig) == 7
        try,
            addpath(path_biosig);
            addpath([ p '..' delimiter 'biosig' ]); % for str2double
            version = [ p '..' delimiter 'biosig' delimiter 'VERSION' ];
            version = loadtxt(version, 'convert', 'off', 'verbose', 'off');
            version = [ version{2,3}(1) '.' version{2,3}(2:end) ];
            disp(['eeglab: adding "BIOSIGv' version '" plugin' ]);
            biosigflag = 1;
        catch
            disp([ 'eeglab: cannot find BIOSIG plugin' ] ); 
            disp([ '   ' lasterr] );
        end;
    end;            
    
    % menu definition
    % --------------- 
    W_MAIN = findobj('tag', 'EEGLAB');
    EEGUSERDAT = get(W_MAIN, 'userdata');
    set(W_MAIN, 'MenuBar', 'none');
    first_m = uimenu( W_MAIN, 'Label', 'File');
	neuromenu = uimenu( first_m, 'Label', 'Import data', 'tag', 'import data'); 
	uimenu( neuromenu, 'Label', 'From ASCII/float file or Matlab array'              ,     'CallBack', [ nocheck '[EEGTMP LASTCOM] = pop_importdata;' e_newnonempty ]);
	uimenu( neuromenu, 'Label', 'From continuous or seg. EGI .RAW file'           , 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_readegi;' e_newnonempty ],  'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Multiple seg. EGI .RAW files', 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_readsegegi;' e_newnonempty ]); 
	uimenu( neuromenu, 'Label', 'From BCI2000 ASCII file'      , 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadbci;' e_newnonempty ],  'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Snapmaster .SMA file'    , 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_snapread;' e_newnonempty ],  'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Neuroscan .CNT file'     , 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadcnt;' e_newnonempty ], 'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Neuroscan .EEG file'     , 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadeeg;' e_newnonempty ]); 
	uimenu( neuromenu, 'Label', 'From ERPSS .RAW or .RDF file',  'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_read_erpss;' e_newnonempty ], 'Separator', 'on'); 
    
    % BIOSIG MENUS
    % ------------
    if biosigflag
        combdf = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_readbdf;' catchstrs.new_non_empty ]; 
        comedf = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_readedf;' catchstrs.new_non_empty ]; 
        combio = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_biosig;'  catchstrs.new_non_empty ]; 
        uimenu( neuromenu, 'Label', 'From Biosemi .BDF file using BIOSIG', 'CallBack', combdf, 'Separator', 'on'); 
        uimenu( neuromenu, 'Label', 'From .EDF file using BIOSIG'        , 'CallBack', comedf, 'tag', 'biosig'); 
        uimenu( neuromenu, 'Label', 'From other formats using BIOSIG'    , 'CallBack', combio); 
    end;
    
    importepoch = uimenu( first_m, 'Label', 'Import epoch info', 'tag', 'import epoch'); 
    uimenu( importepoch, 'Label', 'From Matlab array or ASCII file',        'CallBack', [ checkepoch   '[EEG LASTCOM] = pop_importepoch(EEG);' e_store ]);
	uimenu( importepoch, 'Label', 'From Neuroscan .DAT file', 'CallBack', [ checkepoch   '[EEG LASTCOM]= pop_loaddat(EEG);' e_store]); 
	importevent = uimenu( first_m, 'Label', 'Import event info', 'tag', 'import event'); 
	uimenu( importevent, 'Label', 'From Matlab array or ASCII file',        'CallBack', [ check   '[EEG LASTCOM] = pop_importevent(EEG);' e_store]);
	uimenu( importevent, 'Label', 'From data channel'          , 'CallBack', [ check   '[EEG LASTCOM]= pop_chanevent(EEG);' e_store]); 
	uimenu( importevent, 'Label', 'From Presentation .LOG file'   , 'CallBack', [ check   '[EEG LASTCOM]= pop_importpres(EEG);' e_store]); 
	exportm = uimenu( first_m, 'Label', 'Export', 'tag', 'export'); 
	uimenu( exportm, 'Label', 'Data and ICA activity to text file', 'CallBack', [ check   'LASTCOM = pop_export(EEG);' e_histdone ]);
	uimenu( exportm, 'Label', 'Weight matrix to text file'        , 'CallBack', [ check   'LASTCOM = pop_expica(EEG, ''weights'');' e_histdone ]); 
	uimenu( exportm, 'Label', 'Inverse weight matrix to text file', 'CallBack', [ check   'LASTCOM = pop_expica(EEG, ''inv'');' e_histdone ]); 

	uimenu( first_m, 'Label', 'Load existing dataset' , 'Separator', 'on'   , 'CallBack', [ nocheck '[EEGTMP LASTCOM] = pop_loadset;' e_load_nh]); 
	uimenu( first_m, 'Label', 'Save current dataset'     , 'Separator', 'on', 'CallBack', [ check   '[EEG LASTCOM]    = pop_saveset(EEG);' e_store]);
	uimenu( first_m, 'Label', 'Save datasets'                               , 'CallBack', [ check   '[ALLEEG LASTCOM] = pop_saveset(ALLEEG);' e_hist_nh ]);
	uimenu( first_m, 'Label', 'Clear dataset(s)'                            , 'CallBack', [ nocheck '[ALLEEG LASTCOM] = pop_delset(ALLEEG, -CURRENTSET);' e_hist_nh 'eeglab redraw;' ]);
	uimenu( first_m, 'Label', 'Maximize memory'  , 'Separator', 'on'        , 'CallBack', [ nocheck 'LASTCOM = pop_editoptions;' e_storeall_nh]);
    
	hist_m = uimenu( first_m, 'Label', 'Save history'     , 'Separator', 'on');
	uimenu( hist_m, 'Label', 'Dataset history', 'CallBack', [ nocheck 'LASTCOM = pop_saveh(EEG.history);' e_hist_nh]);
	uimenu( hist_m, 'Label', 'Session history', 'CallBack', [ nocheck 'LASTCOM = pop_saveh(ALLCOM);' e_hist_nh]);
    
	uimenu( first_m, 'Label', 'Quit'             , 'Separator', 'on'        , 'CallBack', ...
	       [ 'close(gcf); disp(''To save the EEGLAB command history  >> pop_saveh(ALLCOM);''); clear global EEG ALLEEG LASTCOM CURRENTSET;']);

second_m = uimenu( W_MAIN, 'Label', 'Edit');
	uimenu( second_m, 'Label', 'Dataset info'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_editset(EEG);' e_store]);
	uimenu( second_m, 'Label', 'Event fields'     , 'CallBack', [ checkevent '[EEG LASTCOM] = pop_editeventfield(EEG);' e_store]);
	uimenu( second_m, 'Label', 'Event values'     , 'CallBack', [ checkevent '[EEG LASTCOM] = pop_editeventvals(EEG);' e_store]);
	uimenu( second_m, 'Label', 'About this dataset', 'CallBack',[ check      '[EEG.comments LASTCOM] =pop_comments(EEG.comments, ''About this dataset'');' e_store]);
	uimenu( second_m, 'Label', 'Channel locations'   , 'CallBack', [ ...
                        'disp(''IMPORTANT: After importing/modifying data channels, you must close'');' ...
                        'disp(''the channel editing window for the changes to take effect in EEGLAB.'');' ...
                        'disp(''TIP: call this function directy from the prompt, ">> pop_chanedit([]);"'');' ...
                        'disp(''     to convert between channel location file formats'');' ...
                        '[TMPCHAN TMPINFO LASTCOM] = pop_chanedit(EEG.chanlocs, EEG.chaninfo); if ~isempty(LASTCOM), EEG.chanlocs = TMPCHAN; EEG.chaninfo = TMPINFO; eeg_checkset(EEG, ''chanlocsize'');' ...
                        'clear TMPCHAN TMPINFO; EEG = h(LASTCOM, EEG);' storecall 'end; eeglab(''redraw'');']);
	uimenu( second_m, 'Label', 'Select data'           , 'CallBack', [ check      '[EEG LASTCOM] = pop_select(EEG);' e_newset], 'Separator', 'on');
	uimenu( second_m, 'Label', 'Select epochs/events'         , 'CallBack', [ checkevent '[EEG TMP LASTCOM] = pop_selectevent(EEG); clear TMP;' e_newset ]);
	uimenu( second_m, 'Label', 'Copy current dataset'  , 'CallBack', [ check      '[ALLEEG LASTCOM] = pop_copyset(ALLEEG, CURRENTSET); eeglab(''redraw'');' e_hist_nh], 'Separator', 'on');
	uimenu( second_m, 'Label', 'Append datasets', 'CallBack', [ check      '[EEGTMP LASTCOM] = pop_mergeset(ALLEEG);' e_newnonempty_nh]);
	uimenu( second_m, 'Label', 'Delete dataset(s)'     , 'CallBack', [ nocheck    '[ALLEEG LASTCOM] = pop_delset(ALLEEG, -CURRENTSET);' e_hist_nh 'eeglab redraw;' ]);
		
fourth_m  = uimenu( W_MAIN, 'Label', 'Tools', 'tag', 'tools');
	uimenu( fourth_m, 'Label', 'Change sampling rate', 'CallBack', [ check      '[EEG LASTCOM] = pop_resample(EEG);' e_newset]);
	uimenu( fourth_m, 'Label', 'Filter the data'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_eegfilt(EEG);' e_newset]);
	uimenu( fourth_m, 'Label', 'Re-reference'   , 'CallBack', [ check      '[EEG LASTCOM] = pop_reref(EEG);' e_newset]);
	uimenu( fourth_m, 'Label', 'Reject continuous data','CallBack',[ checkcont  '[LASTCOM] = pop_eegplot(EEG, 1);' e_hist]);
	uimenu( fourth_m, 'Label', 'Extract epochs'      , 'CallBack', [ check      '[EEG tmp LASTCOM] = pop_epoch(EEG); clear tmp;' e_newset check '[EEG LASTCOM] = pop_rmbase(EEG);' e_store], 'Separator', 'on');
	uimenu( fourth_m, 'Label', 'Remove baseline'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_rmbase(EEG);' e_store]);
	fourth_sub1 = uimenu( fourth_m, 'Label', 'Reject data epochs');
	uimenu( fourth_m, 'Label', 'Run ICA'             , 'CallBack', [ check      '[EEG LASTCOM] = pop_runica(EEG);' e_store], 'foregroundcolor', 'b', 'Separator', 'on');
	uimenu( fourth_m, 'Label', 'Remove components'   , 'CallBack', [ checkica   '[EEG LASTCOM] = pop_subcomp(EEG);' e_newset]);
	fourth_sub2 = uimenu( fourth_m, 'Label', 'Reject data using ICA');

	uimenu( fourth_sub1, 'Label', 'Reject data (all methods)', 'CallBack', [ check      'pop_rejmenu(EEG, 1); LASTCOM = '''';' e_hist]);
	uimenu( fourth_sub1, 'Label', 'Reject by inspection'     , 'CallBack', [ check      '[LASTCOM] = pop_eegplot(EEG, 1);' e_hist]);
	uimenu( fourth_sub1, 'Label', 'Reject extreme values'    , 'CallBack', [ checkepoch '[TMP LASTCOM] = pop_eegthresh(EEG, 1); clear TMP;' e_hist]);
	uimenu( fourth_sub1, 'Label', 'Reject flat line data'    , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_rejtrend(EEG, 1);' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject by probability'    , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_jointprob(EEG, 1);' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject by kurtosis'       , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_rejkurt(EEG, 1);' e_store]);
	uimenu( fourth_sub1, 'Label', 'Reject by spectra'        , 'CallBack', [ checkepoch '[EEG Itmp LASTCOM] = pop_rejspec(EEG, 1); clear Itmp;' e_store]);
	uimenu( fourth_sub1, 'Label', 'Export marks to ICA reject', 'separator', 'on', 'CallBack', [ checkepochica ...
					     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 1,1,1,1,1,1,1,1); h(LASTCOM);' ...
					     'LASTCOM = ''EEG.reject.icarejmanual = EEG.reject.rejglobal;''; eval(LASTCOM);' e_store ]);
	uimenu( fourth_sub1, 'Label', 'Reject marked epochs', 'separator', 'on', 'foregroundcolor', 'b', 'CallBack', [ checkepoch ...
	     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 1,1,1,1,1,1,1,1); EEG = h(LASTCOM, EEG);' ...
	     '[EEG LASTCOM] = pop_rejepoch(EEG);' e_newset]);
	   
	uimenu( fourth_sub2, 'Label', 'Reject components by map', 'CallBack', [ checkicaplot  '[EEG LASTCOM] = pop_selectcomps(EEG);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject data (all methods)','CallBack', [ checkepochica 'pop_rejmenu(EEG, 0); LASTCOM ='''';' e_hist], 'Separator', 'on');
	uimenu( fourth_sub2, 'Label', 'Reject by inspection',     'CallBack', [ checkica      '[LASTCOM] = pop_eegplot(EEG, 0);' e_hist]);
	uimenu( fourth_sub2, 'Label', 'Reject extreme values',    'CallBack', [ checkepochica '[TMP LASTCOM] = pop_eegthresh(EEG, 0); clear TMP;' e_hist]);
	uimenu( fourth_sub2, 'Label', 'Reject flat line activity','CallBack', [ checkepochica '[EEG LASTCOM] = pop_rejtrend(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject by probability',    'CallBack', [ checkepochica '[EEG LASTCOM] = pop_jointprob(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject by kurtosis',       'CallBack', [ checkepochica '[EEG LASTCOM] = pop_rejkurt(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Reject by spectra',        'CallBack', [ checkepochica '[EEG LASTCOM] = pop_rejspec(EEG, 0);' e_store]);
	uimenu( fourth_sub2, 'Label', 'Export marks to data reject', 'separator', 'on', 'CallBack', [ checkepochica ...
					     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); h(LASTCOM);' ...
					     'LASTCOM = ''EEG.reject.rejmanual = EEG.reject.rejglobal;''; eval(LASTCOM);' e_store ]);
	uimenu( fourth_sub2, 'Label', 'Reject marked epochs', 'separator', 'on', 'foregroundcolor', 'b', 'CallBack', [ checkepochica ...
	     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); EEG = h(LASTCOM, EEG);' ...
	     '[EEG LASTCOM] = pop_rejepoch(EEG);' e_newset ]);
   
third_m = uimenu( W_MAIN, 'Label', 'Plot', 'tag', 'plot');
	loc_m = uimenu( third_m, 'Label', 'Channel locations'   );
       uimenu( loc_m, 'Label', 'By name'   , 'CallBack'  , [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''', ''''electrodes'''', ''''labelpoint'''', ''''chaninfo'''', EEG.chaninfo);'']; eval(LASTCOM);' e_hist]);
	   uimenu( loc_m, 'Label', 'By number'   , 'CallBack', [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''', ''''electrodes'''', ''''numpoint'''', ''''chaninfo'''', EEG.chaninfo);'']; eval(LASTCOM);' e_hist]);
    uimenu( third_m, 'Label', 'Channel data (scroll)'        , 'CallBack', [ check          'LASTCOM = pop_eegplot(EEG, 1, 1, 1);' e_hist], 'Separator', 'on');
	uimenu( third_m, 'Label', 'Channel spectra and maps' , 'CallBack', [ check      'LASTCOM = pop_spectopo(EEG, 1);' e_hist]);
	uimenu( third_m, 'Label', 'Channel properties'       , 'CallBack', [ checkplot   'LASTCOM = pop_prop(EEG,1);' e_hist]);
	uimenu( third_m, 'Label', 'Channel ERP image'        , 'CallBack', [ checkepoch  'LASTCOM = pop_erpimage(EEG, 1, h(''find'',''pop_erpimage(EEG,1''));' e_hist]);
	ERP_m = uimenu( third_m, 'Label', 'Channel ERPs');
		uimenu( ERP_m, 'Label', 'With scalp maps'     , 'CallBack', [ checkplot     'LASTCOM = pop_timtopo(EEG);' e_hist]);
		uimenu( ERP_m, 'Label', 'In scalp array'     , 'CallBack', [ checkplot      'LASTCOM = pop_plottopo(EEG);' e_hist]);
		uimenu( ERP_m, 'Label', 'In rect. array'     , 'CallBack', [ checkepoch     '[tmpeeg LASTCOM] = pop_plotdata(EEG, 1); clear tmpeeg;' e_hist]);
	topo_m = uimenu( third_m, 'Label', 'ERP map series');
		uimenu( topo_m, 'Label', 'In 2-D'     , 'CallBack', [ checkplot      'LASTCOM = pop_topoplot(EEG, 1);' e_hist]);
		uimenu( topo_m, 'Label', 'In 3-D'     , 'CallBack', [ checkplot      '[EEG LASTCOM] = pop_headplot(EEG, 1);' e_store]);
	uimenu( third_m, 'Label', 'Sum/Compare ERPs'             , 'CallBack', [ checkepoch     'LASTCOM = pop_comperp(ALLEEG);' e_hist]);

    uimenu( third_m, 'Label', 'Component activations (scroll)', 'CallBack', [ checkica  '[LASTCOM] = pop_eegplot(EEG, 0, 1, 1);' e_hist],'Separator', 'on');
	uimenu( third_m, 'Label', 'Component spectra and maps' , 'CallBack', [ checkicaplot   'LASTCOM = pop_spectopo(EEG, 0);' e_hist]);
	topoica_m = uimenu( third_m, 'Label', 'Component maps');
		uimenu( topoica_m, 'Label', 'In 2-D', 'CallBack', [ checkicaplot   'LASTCOM = pop_topoplot(EEG, 0);' e_hist]);
		uimenu( topoica_m, 'Label', 'In 3-D', 'CallBack', [ checkicaplot   '[EEG LASTCOM] = pop_headplot(EEG, 0);' e_store]);
	uimenu( third_m, 'Label', 'Component properties'     , 'CallBack', [ checkicaplot   'LASTCOM = pop_prop(EEG,0);' e_hist]);
	uimenu( third_m, 'Label', 'Component ERP image'      , 'CallBack', [ checkepochica  'LASTCOM = pop_erpimage(EEG, 0, h(''find'',''pop_erpimage(EEG,0''));' e_hist]);
	ERPC_m = uimenu( third_m, 'Label', 'Component ERPs');
	   uimenu( ERPC_m, 'Label', 'With component maps', 'CallBack', [ checkica 'LASTCOM = pop_envtopo(EEG);' e_hist]);
	   uimenu( ERPC_m, 'Label', 'With comp. maps (compare)', 'CallBack', [ checkica 'if length(ALLEEG) == 1, error(''Need at least 2 datasets''); end; LASTCOM = pop_envtopo(ALLEEG);' e_hist]);
	   uimenu( ERPC_m, 'Label', 'In rectangular array'      , 'CallBack', [ checkepochica     '[tmpeeg LASTCOM] = pop_plotdata(EEG, 0); clear tmpeeg;' e_hist]);
	uimenu( third_m, 'Label', 'Sum/Compare comp. ERPs'      , 'CallBack', [ checkepochica     'LASTCOM = pop_comperp(ALLEEG, 0);' e_hist]);

	stat_m = uimenu( third_m, 'Label', 'Data statistics', 'Separator', 'on'); %, 'enable', fastif(exist('kstest'), 'on', 'off'));
	uimenu( stat_m, 'Label', 'Channel statistics'       , 'CallBack', [ check          'LASTCOM = pop_signalstat(EEG, 1);' e_hist]);
	uimenu( stat_m, 'Label', 'Component statistics'     , 'CallBack', [ checkica       'LASTCOM = pop_signalstat(EEG, 0);' e_hist]);
	uimenu( stat_m, 'Label', 'Event statistics'         , 'CallBack', [ checkevent     'LASTCOM = pop_eventstat(EEG);' e_hist]);
	spec_m = uimenu( third_m, 'Label', 'Time-frequency transforms', 'Separator', 'on');
		uimenu( spec_m, 'Label', 'Channel time-frequency'   , 'CallBack', [ checkepoch    'LASTCOM = pop_timef(EEG, 1, h(''find'',''pop_timef(EEG,1''));' e_hist]);
		uimenu( spec_m, 'Label', 'Channel cross-coherence'  , 'CallBack', [ checkepoch    'LASTCOM = pop_crossf(EEG, 1,h(''find'',''pop_crossf(EEG,1''));' e_hist]);
		uimenu( spec_m, 'Label', 'Component time-frequency' , 'CallBack', [ checkepochica 'LASTCOM = pop_timef(EEG, 0, h(''find'',''pop_timef(EEG,0''));' e_hist],'Separator', 'on');
		uimenu( spec_m, 'Label', 'Component cross-coherence', 'CallBack', [ checkepochica 'LASTCOM = pop_crossf(EEG, 0,h(''find'',''pop_crossf(EEG,0''));' e_hist]);
		
    set_m   = uimenu( W_MAIN, 'Label', 'Datasets');
    help_m  = uimenu( W_MAIN, 'Label', 'Help');
    uimenu( help_m, 'Label', 'About EEGLAB', 'CallBack', 'pophelp(''eeglab'');');
    uimenu( help_m, 'Label', 'About EEGLAB help', 'CallBack', 'pophelp(''eeg_helphelp'');');
    uimenu( help_m, 'Label', 'EEGLAB license', 'CallBack', 'pophelp(''eeglablicense.txt'', 1);');
    uimenu( help_m, 'Label', 'EEGLAB menus', 'CallBack', 'eeg_helpmenu;','separator','on');
    help_subm1 = uimenu( help_m, 'Label', 'EEGLAB functions');
    uimenu( help_subm1, 'Label', 'Toolbox functions', 'CallBack', 'pophelp(''ica'');');
	uimenu( help_subm1, 'Label', 'Signal processing functions', 'callback', 'eeg_helpsigproc;');	
	uimenu( help_subm1, 'Label', 'Interactive pop_ functions', 'callback', 'eeg_helppop;');	
    help_subm2 = uimenu( help_m, 'Label', 'EEGLAB advanced');
    uimenu( help_subm2, 'Label', 'Dataset structure', 'CallBack', 'pophelp(''eeg_checkset'');');
	uimenu( help_subm2, 'Label', 'Admin functions', 'callback', 'eeg_helpadmin;');	
    uimenu( help_m, 'Label', 'Tutorial (web)', 'CallBack', 'tutorial;');
    uimenu( help_m, 'Label', 'Contact us (email)', 'CallBack', 'web(''mailto:eeglab@sccn.ucsd.edu'');');

    % looking for eeglab plugins
    % --------------------------
    p = which('eeglab.m');
    p = p(1:findstr(p,'eeglab.m')-1);
    dircontent1 = what(p);
    dircontent  = dir([ p 'plugins' ]);
    delimiter = p(end); if strcmpi(delimiter, ':'), delimiter = '::'; end;
    dircontent  = { dircontent1.m{:} dircontent.name };

    % scan plugin folder
    % ------------------
    for index = 1:length(dircontent)

        % find function
        % -------------
        funcname = '';
        if exist([p 'plugins' delimiter dircontent{index}]) == 7
            if ~strcmpi(dircontent{index}, '.') & ~strcmpi(dircontent{index}, '..')
                addpath([ p 'plugins' delimiter dircontent{index} ]);
                tmpdir = dir([ p 'plugins' delimiter dircontent{index}]);
                for tmpind = 1:length(tmpdir)
                    % find plugin function in subfolder
                    % ---------------------------------
                    if ~isempty(findstr(tmpdir(tmpind).name, 'eegplugin')) & tmpdir(tmpind).name(end) == 'm'
                        funcname = tmpdir(tmpind).name(1:end-2);
                        tmpind = length(tmpdir);
                    end;
                    
                    % spetial case of eeglab subfolder (for BIOSIG)
                    % --------------------------------
                    if strcmpi(tmpdir(tmpind).name, 'eeglab')
                        addpath([ p 'plugins' delimiter dircontent{index} delimiter 'eeglab' ],'-end');
                        tmpdir2 = dir([ p 'plugins' delimiter dircontent{index} delimiter 'eeglab' ]);
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
                disp(['eeglab: adding "' vers '" plugin (see >> help ' funcname ')' ]);    
           catch
                try,
                    eval( [ funcname '(gcf, trystrs, catchstrs)' ]);
                    disp(['eeglab: adding plugin function "' funcname '"' ]);    
                catch
                    disp([ 'eeglab: error while adding plugin "' funcname '"' ] ); 
                    disp([ '   ' lasterr] );
                end;
            end;
        end;
    end;

    % changing plugin menu color
    % --------------------------
    fourthsub_m = findobj('parent', fourth_m);
    plotsub_m   = findobj('parent', third_m);
    importsub_m = findobj('parent', neuromenu);
    epochsub_m  = findobj('parent', importepoch);
    eventsub_m  = findobj('parent', importevent);
    exportsub_m = findobj('parent', exportm);
    icadefs; % containing PLUGINMENUCOLOR
    if length(fourthsub_m) > 10, set(fourthsub_m(1:end-10), 'foregroundcolor', PLUGINMENUCOLOR); end;
    if length(plotsub_m)   > 17, set(plotsub_m  (1:end-17), 'foregroundcolor', PLUGINMENUCOLOR); end;
    if length(importsub_m) > 8,  set(importsub_m(1:end-8) , 'foregroundcolor', PLUGINMENUCOLOR); end;
    if length(epochsub_m ) > 2 , set(epochsub_m (1:end-2 ), 'foregroundcolor', PLUGINMENUCOLOR); end;
    if length(eventsub_m ) > 3 , set(eventsub_m (1:end-3 ), 'foregroundcolor', PLUGINMENUCOLOR); end;
    if length(exportsub_m) > 3 , set(exportsub_m(1:end-3 ), 'foregroundcolor', PLUGINMENUCOLOR); end;

EEGMENU = uimenu( set_m, 'Label', '------', 'Enable', 'off');
set(W_MAIN, 'userdat', { EEGUSERDAT{1} EEGMENU });
eeglab('redraw');
if nargout < 1
    clear ALLEEG;
end;

% REMOVED MENUS
	%uimenu( fourth_m, 'Label', 'Automatic comp. reject',  'enable', 'off', 'CallBack', '[EEG LASTCOM] = pop_rejcomp(EEG); h(LASTCOM); if ~isempty(LASTCOM), eeg_store(CURRENTSET); end;');
	%uimenu( fourth_m, 'Label', 'Reject (synthesis)' , 'Separator', 'on', 'CallBack', '[EEG LASTCOM] = pop_rejall(EEG); h(LASTCOM); if ~isempty(LASTCOM), eeg_store; end; eeglab(''redraw'');');

% --------------------
% draw the main figure
% --------------------
function eeg_mainfig(onearg);

icadefs;
COLOR = BACKEEGLABCOLOR;
WINMINX         = 17;
WINMAXX         = 260;
WINYDEC			= 13;
NBLINES         = 16;
WINY		    = WINYDEC*NBLINES;

BORDERINT       = 4;
BORDEREXT       = 10;
FONTNAME        = 'courrier';
FONTSIZE        = 11;

hh = findobj('tag', 'EEGLAB');
if ~isempty(hh)
    disp('EEGLAB warning: there can be only one EEGLAB window, closing old one');
    close(hh);
end;
if strcmpi(onearg, 'remote')
    figure(	'name', [ 'EEGLAB v' EEGLAB_VERSION ], ... 
	'numbertitle', 'off', ...
	'resize', 'off', ...
	'Position',[200 100 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) 30 ], ...
	'color', COLOR, ...
	'Tag','EEGLAB', ...
	'Userdata', {[] []});
    return;
end;

W_MAIN = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'name', [ 'EEGLAB v' EEGLAB_VERSION ], ... 
	'numbertitle', 'off', ...
	'resize', 'off', ...
	'Position',[200 100 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) (WINY+2*BORDERINT+2*BORDEREXT) ], ...
	'color', COLOR, ...
	'Tag','EEGLAB', ...
    'visible', 'off', ...   
	'Userdata', {[] []});
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
set(W_MAIN, 'userdata', {[] []}, 'visible', 'on');
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
eeg_options;
if ~option_keepdataset
	if ~isempty(ALLEEG)
		if popask( strvcat('Remove all datasets except the present one ?', ...
						   'Otherwise, under "File > Maximize memory" clear dataset overwrite'))
			ALLEEG = []; CURRENTSET = 0;
			h('ALLEEG = []; CURRENTSET = 0;');
		else 
			return;
		end;
	end;
	set(findobj('parent', gcf, 'label', 'Datasets'), 'enable', 'off');
	CURRENTSET = 0;
else
	if isempty(ALLEEG) & ~isempty(EEG) & ~isempty(EEG.data)
		ALLEEG = EEG;
	end;
	set(findobj('parent', gcf, 'label', 'Datasets'), 'enable', 'on');
end;

% setting the dataset menu
% ------------------------
while( index <= MAX_SET)
	try
		set( EEGMENU(index), 'Label', '------', 'checked', 'off');
	catch, EEGMENU(index) = uimenu( set_m, 'Label', '------', 'Enable', 'on'); end;	
	set( EEGMENU(index), 'Enable', 'on' );
	try, ALLEEG(index).data;
		if ~isempty( ALLEEG(index).data)
       		menutitle   = sprintf('Dataset %d:%s', index, ALLEEG(index).setname);
			set( EEGMENU(index), 'Label', menutitle);
			set( EEGMENU(index), 'CallBack', ['com = ''EEG = eeg_retrieve(ALLEEG, ' int2str(index) ...
					'); CURRENTSET = ' int2str(index) ';''; eval(com); h(com); eeglab(''redraw'');' ]);
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
                  'CURRENTSET = nonempty(tmpind);' ...
                  'EEG = ALLEEG(CURRENTSET);' ...
                  'eeglab(''redraw'');' ];
    if MAX_SET == length(EEGMENU), EEGMENU(end+1) = uimenu( set_m, 'Label', '------', 'Enable', 'on'); end;
    
    set(EEGMENU(end), 'enable', 'on', 'Label', 'Select multiple datasets', ...
                      'callback', cb_select, 'separator', 'on');
end;

EEGUSERDAT{2} = EEGMENU;
set(W_MAIN, 'userdata', EEGUSERDAT);

if option_keepdataset & (isempty(CURRENTSET) | length(ALLEEG) < CURRENTSET(1) | CURRENTSET(1) == 0 | isempty(ALLEEG(CURRENTSET(1)).data))
	CURRENTSET = 0;
	for index = 1:length(ALLEEG)
		if ~isempty(ALLEEG(index).data)
			CURRENTSET = index;
			break;
		end;
	end;
	if CURRENTSET ~= 0
		h([ 'EEG = eeg_retrieve(ALLEEG,' int2str(CURRENTSET) '); CURRENTSET = ' int2str(CURRENTSET) ';'])
		EEG = eeg_retrieve(ALLEEG, CURRENTSET);	
	else 
		EEG = eeg_emptyset;
	end;
end;

if (isempty(EEG) | isempty(EEG(1).data)) & CURRENTSET(1) ~= 0 & option_keepdataset
	h([ 'EEG = eeg_retrieve(ALLEEG,' int2str(CURRENTSET) '); CURRENTSET = ' int2str(CURRENTSET) ';'])
	EEG = eeg_retrieve(ALLEEG, CURRENTSET);	
end;
% test if dataset has changed
% ---------------------------
if option_keepdataset
    if ~isempty(ALLEEG) & CURRENTSET~= 0 & ~isequal(EEG.data, ALLEEG(CURRENTSET).data) 
        % the above comparison does not work for ome structures
        tmpanswer = questdlg2(strvcat('The current EEG dataset has changed. What should eeglab do with the changes?', ' '), ...
                              'Dataset change detected', ...
                              'Keep changes', 'Delete changes', 'New dataset', 'Make new dataset');
        
        if tmpanswer(1) == 'D' % delete changes
            EEG = eeg_retrieve( ALLEEG, CURRENTSET);
            h('EEG = eeg_retrieve( ALLEEG, CURRENTSET);');
        elseif tmpanswer(1) == 'K' % keep changes
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            h('[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);');
        else % make new dataset
            [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); 
            h(LASTCOM);
            MAX_SET = max(length( ALLEEG ), length(EEGMENU));
        end;
    end;
end;

% print some informations on the main figure
% ------------------------------------------
g = myguihandles(gcf);
if ~isfield(g, 'win0') % no display
    return;
end;

if (exist('EEG') == 1) & isstruct(EEG) & ~isempty(EEG(1).data)        
    hh = findobj('parent', gcf, 'userdata', 'fullline'); set(hh, 'visible', 'off');
    hh = findobj('parent', gcf, 'userdata', 'datinfo');  set(hh, 'visible', 'on');
    
    if length(EEG) > 1 % several datasets

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
        
        % epoch consistency
        % -----------------
        epochconsist = 'no';
        if strcmpi(datasettype, 'epoched');
            allpnts = unique( [ EEG.pnts ] );
            allxmin = unique( [ EEG.xmin ] );
            if length(allpnts) == 1 & length(allxmin) == 1, epochconsist = 'yes'; end;
        end;
        
        % channels
        % --------
        chanlen    = unique( [ EEG.nbchan ] );
        chanlenstr = vararg2str( mattocell(chanlen) );
        
        % channel consistency
        % -------------------
        chanconsist = 'yes';
        anyempty    = unique( cellfun( 'isempty', { EEG.chanlocs }) );
        if length(chanlen) == 1 & all(anyempty == 0)
            channame1 = { EEG(1).chanlocs.labels };
            for i = 2:length(EEG)
                channame2 = { EEG(i).chanlocs.labels };
                if length(intersect(channame1, channame2)) ~= length(channame1)
                    chanconsist = 'no';
                end;
            end;
        else
            chanconsist = 'no';
        end;

        % channel location
        % ----------------
        if length(anyempty) == 2,   chanlocs = 'mixed, yes and no';
        elseif anyempty == 0,       chanlocs = 'yes';
        else                        chanlocs = 'no';
        end;
        
        % total number of events
        % ----------------------
        totevents    = num2str(sum( cellfun( 'length', { EEG.event }) ));
 
        % sampling rate
        % -------------
        srate = vararg2str( mattocell( unique( [ EEG.srate ] ) ));
        
        % ica weights
        % -----------
        anyempty    = unique( cellfun( 'isempty', { EEG.icaweights }) );
        if length(anyempty) == 2,   icaweights = 'mixed, yes and no';
        elseif anyempty == 0,       icaweights = 'yes';
        else                        icaweights = 'no';
        end;
        
        % ica consistency
        % ---------------
        icaconsist = 'yes';
        if strcmpi(icaweights, 'yes') & strcmpi(chanconsist, 'yes')
            ica1 = EEG(1).icaweights;
            for i = 2:length(EEG)
                if ~isequal(EEG(1).icaweights, EEG(i).icaweights)
                    icaconsist = 'no';
                end;
            end;
        else
            icaconsist = 'no';
        end;

        % total size
        % ----------
        totsize = whos('EEG');
        
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
        set( g.val5, 'String', chanlen);
        set( g.val6, 'String', chanconsist);
        set( g.val7, 'String', chanlocs);
        set( g.val8, 'String', totevents);
        set( g.val9, 'String', srate);
        set( g.val10, 'String', icaweights);
        set( g.val11, 'String', icaconsist);
        set( g.val12, 'String', num2str(round(totsize.bytes/1E6*10)/10));        
        
        % disable menus
        % -------------
        set( findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'File')    , 'enable', 'off');
        set( findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Edit')    , 'enable', 'off');
        set( findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Plot')    , 'enable', 'off');
        set( findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Datasets'), 'enable', 'on');
        tool_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Tools');
        set( findobj('parent', tool_m, 'type', 'uimenu'), 'enable', 'off');
        set( findobj('parent', tool_m, 'type', 'uimenu', 'label', 'Run ICA'), 'enable', 'on');
       
    else % one dataset selected
        
        % text
        % ----
        set( g.win2, 'String', 'Channels per frame');
        set( g.win3, 'String', 'Frames per epoch');
        set( g.win4, 'String', 'Epochs');
        set( g.win5, 'String', 'Events');
        set( g.win6, 'String', 'Sampling rate (Hz)');
        set( g.win7, 'String', 'Epoch start (sec)');
        set( g.win8, 'String', 'Epoch end (sec)');
        set( g.win9, 'String', 'Average reference');
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

        fullfilename = [ EEG.filepath EEG.filename];
        if ~isempty(fullfilename)
            if length(fullfilename) > 30
                set( g.win1, 'String', sprintf('Filename: ...%s\n', fullfilename(max(1,length(fullfilename)-30):end) ));
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

        set( g.val9, 'String', fastif(strcmpi(EEG.ref, 'averef'), 'Yes', 'No'));
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
        set( g.val12, 'String', num2str(round(tmp.bytes/1E6*10)/10));
    
        % enable menus
        % ------------
        file_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'File');  set(file_m, 'enable', 'on');
        edit_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Edit');  set(edit_m, 'enable', 'on');
        tool_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Tools'); set(tool_m, 'enable', 'on');
        plot_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Plot');  set(plot_m, 'enable', 'on');
        data_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Plot');  set(data_m, 'enable', 'on');
        set( findobj('parent', file_m, 'type', 'uimenu'), 'enable', 'on');
        set( findobj('parent', edit_m, 'type', 'uimenu'), 'enable', 'on');
        set( findobj('parent', plot_m, 'type', 'uimenu'), 'enable', 'on');
        set( findobj('parent', tool_m, 'type', 'uimenu'), 'enable', 'on');

    end;
else
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
    
    % disable menus
    % -------------
    file_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'File');  set(file_m, 'enable', 'on');
    set( findobj('parent', file_m, 'type', 'uimenu'), 'enable', 'off');
    set( findobj('parent', file_m, 'type', 'uimenu', 'label', 'Import data')          , 'enable', 'on');
    set( findobj('parent', file_m, 'type', 'uimenu', 'label', 'Load existing dataset'), 'enable', 'on');
    set( findobj('parent', file_m, 'type', 'uimenu', 'label', 'Maximize memory')      , 'enable', 'on');
    set( findobj('parent', file_m, 'type', 'uimenu', 'label', 'Quit')                 , 'enable', 'on');
    edit_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Edit');      set(edit_m, 'enable', 'off');
    tool_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Tools');     set(tool_m, 'enable', 'off');
    plot_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Plot');      set(plot_m, 'enable', 'off');
    data_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'Datasets');  set(data_m, 'enable', 'off');
    
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
        if length(tmpp) > length(tmpnewpath), tmpp = tmpp(1:end-1); end; % remove trailing delimiter
        if length(tmpp) > length(tmpnewpath), tmpp = tmpp(1:end-1); end; % remove trailing delimiter
        %disp([ tmpp '     |        ' tmpnewpath '(' num2str(~strcmpi(tmpnewpath, tmpp)) ')' ]);
        if ~strcmpi(tmpnewpath, tmpp)
            addpath(tmpnewpath);
        end;
    else
        %disp([ 'Adding new path ' tmpnewpath ]);
        addpath(tmpnewpath);
    end;
    
