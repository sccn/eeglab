% eeglab() - EEGLAB (v4.0) Matlab graphic user interface environment for 
%   electrophysiological data analysis incorporating the ICA EEG toolbox 
%   (Makeig et al.) developed at CNL / The Salk Institute, 1997-2001. 
%   Released 2002- as EEGLAB (Delorme, Makeig, et al.) at the Swartz Center 
%   for Computational Neuroscience, Institute for Neural Computation, 
%   University of California San Diego (http://sccn.ucsd.edu/). 
%   User feedback welcome: email eeglab@sccn.ucsd.edu
%
% Authors: Arnaud Delorme, Scott Makeig, et al.
%
% Description:
%   EEGLAB is Matlab software for processing continuous or epoched event-related 
%   EEG or other physiological data. It is designed for use by both novice and 
%   expert Matlab users. In normal use, the EEGLAB graphic interface calls 
%   graphic functions via pop-up function windows. the EEGLAB history mechanism 
%   can save the resulting Matlab calls to disk for later incorporation into 
%   Matlab scripts.  A single data structure ('EEG') containing all dataset 
%   parameters may be accessed and modified directly from the Matlab commandline. 
%
% Usage: 1) To (re)start EEGLAB from scratch, type
%            >> eeglab           % Ignores any loaded datasets
%        2) To redaw and update the EEGLAB interface, type
%            >> eeglab redraw    % Scans for non-empty datasets
%            >> eeglab rebuild   % Closes and rebuilds the EEGLAB window
%
% See: license.txt: The software distribution license 
%      http://sccn.ucsd.edu/eeglab/tutorial/: The EEGLAB tutorial
%      >> help eeg_checkset(), % Shows the structure of the EEG dataset 
%
% Main files:
% ---------- 
% eeglab()        - main graphic interface
% license.tx      - GNU license
% 
% Functions recently added to EEGLAB: 
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
% realproba()     - compute observed probability (used by entropy)
% rejkurt()       - calculate and reject data based on kurtosis
% rejtrend()      - reject EEG showing linear trends  !!!
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

function eeglab( onearg )
eeg_options; 
eeg_global;

evalin('base', 'eeg_global;');
if nargin < 1 | exist('EEG') ~= 1
	clear global EEG ALLEEG CURRENTSET ALLCOM LASTCOM;
	eeg_global;
	EEG = eeg_emptyset;
	h('eeglab;');
end;

besamenu = 0;
if nargin == 1
	if strcmp(onearg, 'redraw')
		W_MAIN = findobj('tag', 'EEGLAB');
		if ~isempty(W_MAIN)
			updatemenu;
			return;
		else
			h('eeglab rebuild;');
		end;
	elseif strcmp(onearg, 'besa');
		besamenu = 1;
		if exist('ALLEEG') == 0
			clear global EEG ALLEEG CURRENTSET ALLCOM LASTCOM;
			EEG = eeg_emptyset;
		end;
		eeg_global;
		h('eeglab besa;');
		disp('Besa menu activated');
	else
		h('eeglab rebuild;');
	end;
end;

colordef white

% checking strings
% ----------------
e_try = 'try,';
e_catch = 'catch, errordlg2(lasterr, ''EEGLAB error''); LASTCOM= ''''; clear EEGTMP; end;';
nocheck           = e_try;
check             = ['[EEG LASTCOM] = eeg_checkset(EEG, ''data''); h(LASTCOM);' e_try];
checkcont         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''contdata''); h(LASTCOM);' e_try];
checkica          = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica''); h(LASTCOM);' e_try];
checkepoch        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch''); h(LASTCOM);' e_try];
checkevent        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''event''); h(LASTCOM);' e_try];
checkepochica     = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica''); h(LASTCOM);' e_try];
checkplot         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''chanloc''); h(LASTCOM);' e_try];
checkicaplot      = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'', ''chanloc''); h(LASTCOM);' e_try];
checkepochplot    = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''chanloc''); h(LASTCOM);' e_try];
checkepochicaplot = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'', ''chanloc''); h(LASTCOM);' e_try];

storecall    = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); h(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);'');';
storeload    = '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); h(''[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);'');';
storenewcall = '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET); h(LASTCOM);';
storeallcall = 'ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_checkset(EEG); h(''ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_checkset(EEG);'');';

e_newnonempty     = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG), EEG = EEGTMP;' storenewcall 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_load            = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG), EEG = EEGTMP;' storeload 'disp(''Done.''); end;  clear EEGTMP; eeglab(''redraw'');'];
e_newset          = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG),' storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store           = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG),' storecall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_storeall        = [e_catch 'h(LASTCOM); if ~isempty(LASTCOM) & ~isempty(EEG),' storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist            = [e_catch 'h(LASTCOM); eeglab(''redraw'');'];

% menu definition
% --------------- 
eeg_mainfig;
W_MAIN = findobj('tag', 'EEGLAB');
EEGUSERDAT = get(W_MAIN, 'userdata');
set(W_MAIN, 'MenuBar', 'none');
first_m = uimenu( W_MAIN, 'Label', 'File');
	neuromenu = uimenu( first_m, 'Label', 'Import data'); 
	uimenu( neuromenu, 'Label', 'From ASCII/float file or Matlab array'              ,     'CallBack', [ nocheck '[EEGTMP LASTCOM] = pop_importdata;' e_newnonempty ]);
	uimenu( neuromenu, 'Label', 'From BCI2000 ASCII file'    ,     'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadbci;' e_newnonempty ],  'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Snapmaster .SMA file'       ,     'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_snapread;' e_newnonempty ],  'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Biosemi .EDF file'             ,  'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_readedf;' e_newnonempty ], 'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Neuroscan .CNT file',  'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadcnt;' e_newnonempty ], 'Separator', 'on'); 
	uimenu( neuromenu, 'Label', 'From Neuroscan .EEG file'  ,    'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadeeg;' e_newnonempty ]); 

	importepoch = uimenu( first_m, 'Label', 'Import epoch info'); 
    uimenu( importepoch, 'Label', 'From Matlab array or ASCII file',        'CallBack', [ check   '[EEG LASTCOM] = pop_importepoch(EEG);' e_store ]);
	uimenu( importepoch, 'Label', 'From Neuroscan .DAT file', 'CallBack', [ check   '[EEG LASTCOM]= pop_loaddat(EEG);' e_store]); 
	importevent = uimenu( first_m, 'Label', 'Import event info'); 
	uimenu( importevent, 'Label', 'From Matlab array or ASCII file',        'CallBack', [ check   '[EEG LASTCOM] = pop_importevent(EEG);' e_store]);
	uimenu( importevent, 'Label', 'From data channel'          , 'CallBack', [ check   '[EEG LASTCOM]= pop_chanevent(EEG);' e_store]); 
	uimenu( importevent, 'Label', 'From Presentation .LOG file'   , 'CallBack', [ check   '[EEG LASTCOM]= pop_importpres(EEG);' e_store]); 

	uimenu( first_m, 'Label', 'Load existing dataset' , 'Separator', 'on'   , 'CallBack', [ nocheck '[EEGTMP LASTCOM]= pop_loadset;' e_load]); 
	uimenu( first_m, 'Label', 'Save current dataset'     , 'Separator', 'on', 'CallBack', [ check   '[EEG LASTCOM] = pop_saveset(EEG);' e_store]);
	uimenu( first_m, 'Label', 'Save datasets'                               , 'CallBack', [ check   '[ALLEEG LASTCOM] = pop_saveset(ALLEEG);' e_hist ]);
	uimenu( first_m, 'Label', 'Clear dataset(s)'                            , 'CallBack', [ nocheck '[ALLEEG LASTCOM] = pop_delset(ALLEEG, -CURRENTSET);' e_hist ]);
	uimenu( first_m, 'Label', 'Maximize memory'  , 'Separator', 'on'        , 'CallBack', [ nocheck 'LASTCOM = pop_editoptions;' e_storeall]);
	uimenu( first_m, 'Label', 'Save history'     , 'Separator', 'on'        , 'CallBack', [ nocheck 'LASTCOM = pop_saveh(ALLCOM);' e_hist]);
	uimenu( first_m, 'Label', 'Quit'             , 'Separator', 'on'        , 'CallBack', ...
	       [ 'close(gcf); disp(''To save the EEGLAB command history  >> pop_saveh(ALLCOM);''); clear global EEG ALLEEG LASTCOM CURRENTSET;']);

second_m = uimenu( W_MAIN, 'Label', 'Edit');
	uimenu( second_m, 'Label', 'Dataset info'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_editset(EEG);' e_store]);
	uimenu( second_m, 'Label', 'Event fields'     , 'CallBack', [ checkevent '[EEG LASTCOM] = pop_editeventfield(EEG);' e_store]);
	uimenu( second_m, 'Label', 'Event values'     , 'CallBack', [ checkevent '[EEG LASTCOM] = pop_editeventvals(EEG);' e_store]);
	uimenu( second_m, 'Label', 'About this dataset', 'CallBack', [ check      '[EEG.comments LASTCOM] =pop_comments(EEG.comments, ''About this dataset'');' e_store]);
	uimenu( second_m, 'Label', 'Channel locations'   , 'CallBack', [ 'disp(''IMPORTANT: Close the channel editing window to import channels'');' ...
                        'disp(''WARNING: the number of channel information must match the number of'');' ... 
                        'disp(''         data channels (otherwise channel''s info is ignored)'');' ...
                        'disp(''TIP: to edit channel info only type "chanlocs = pop_chanedit([]);" from the command line'');' ...
    '[TMPCHAN LASTCOM] =pop_chanedit(EEG.chanlocs); if ~isempty(LASTCOM), EEG.chanlocs = TMPCHAN; clear TMPCHAN; h(LASTCOM);' storecall 'end; eeglab(''redraw'');']);
	uimenu( second_m, 'Label', 'Select data'           , 'CallBack', [ check      '[EEG LASTCOM] = pop_select(EEG);' e_newset], 'Separator', 'on');
	uimenu( second_m, 'Label', 'Select events'         , 'CallBack', [ checkevent '[EEG TMP LASTCOM] = pop_selectevent(EEG); clear TMP;' e_newset ]);
	uimenu( second_m, 'Label', 'Copy current dataset'  , 'CallBack', [ check      '[ALLEEG LASTCOM] = pop_copyset(ALLEEG, CURRENTSET); h(LASTCOM); eeglab(''redraw'');' e_hist], 'Separator', 'on');
	uimenu( second_m, 'Label', 'Append datasets', 'CallBack', [ check      '[EEGTMP LASTCOM] = pop_mergeset(ALLEEG);' e_newnonempty]);
	uimenu( second_m, 'Label', 'Delete dataset(s)'     , 'CallBack', [ nocheck    '[ALLEEG LASTCOM] = pop_delset(ALLEEG, -CURRENTSET);' e_hist]);
		
fourth_m  = uimenu( W_MAIN, 'Label', 'Tools');
	uimenu( fourth_m, 'Label', 'Change sampling rate', 'CallBack', [ check      '[EEG LASTCOM] = pop_resample(EEG);' e_newset], 'enable', fastif(exist('resample'), 'on', 'off'));
	uimenu( fourth_m, 'Label', 'Filter the data'     , 'CallBack', [ check      '[EEG LASTCOM] = pop_eegfilt(EEG);' e_newset], 'enable', fastif(exist('resample'), 'on', 'off'));
	uimenu( fourth_m, 'Label', 'Average reference'   , 'CallBack', [ check      '[EEG LASTCOM] = pop_averef(EEG,1);' e_store]);
	uimenu( fourth_m, 'Label', 'Reject continuous data','CallBack',[ checkcont  '[LASTCOM] = pop_eegplot(EEG, 1);' e_hist]);
	uimenu( fourth_m, 'Label', 'Extract epochs'      , 'CallBack', [ check      '[EEG tmp LASTCOM] = pop_epoch(EEG); clear tmp;' e_newset], 'Separator', 'on');
	uimenu( fourth_m, 'Label', 'Remove baseline'     , 'CallBack', [ checkepoch '[EEG LASTCOM] = pop_rmbase(EEG);' e_store]);
	fourth_sub1 = uimenu( fourth_m, 'Label', 'Reject data epochs');
	uimenu( fourth_m, 'Label', 'Run ICA'             , 'CallBack', [ check      '[EEG LASTCOM] = pop_runica(EEG);' e_store], 'foregroundcolor', 'b', 'Separator', 'on');
	uimenu( fourth_m, 'Label', 'Remove components'   , 'CallBack', [ checkica   '[EEG LASTCOM] = pop_subcomp(EEG);' e_newset]);
	fourth_sub2 = uimenu( fourth_m, 'Label', 'Reject data using ICA');

if besamenu
	fourth_sub3 = uimenu( fourth_m, 'Label', 'Localize dipoles using BESA');
	uimenu( fourth_sub3, 'Label', 'Export dipoles'   , 'CallBack', [ 'besaexport(EEG);']);
	uimenu( fourth_sub3, 'Label', 'Import dipoles'   , 'CallBack', [ check 'EEG = besaimport(EEG);' e_store]);
	uimenu( fourth_sub3, 'Label', 'Plot dipoles BESA'   , 'CallBack', [ 'besaplot(EEG.sources);']);
	uimenu( fourth_sub3, 'Label', 'Plot dipoles MRI'  , 'CallBack', [ 'besaplot(EEG.sources, ''image'', ''mri'');']);
	uimenu( fourth_sub3, 'Label', 'Plot dipoles summary BESA', 'CallBack', [ 'besaplot(EEG.sources, ''summary'', ''on'', ''dipolesize'', 30, ''mesh'', ''off'');']);
	uimenu( fourth_sub3, 'Label', 'Plot dipoles summary MRI', 'CallBack', [ 'besaplot(EEG.sources, ''image'', ''mri'', ''summary'', ''on'', ''dipolesize'', 30, ''mesh'', ''off'');']);
end;

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
	uimenu( fourth_sub2, 'Label', 'Export marks to data reject', 'separator', 'on', 'CallBack', [ checkepochica ...
					     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); h(LASTCOM);' ...
					     'LASTCOM = ''EEG.reject.rejmanual = EEG.reject.rejglobal;''; eval(LASTCOM);' e_store ]);
	uimenu( fourth_sub2, 'Label', 'Reject marked epochs', 'separator', 'on', 'foregroundcolor', 'b', 'CallBack', [ checkepochica ...
	     '[EEG LASTCOM] = eeg_rejsuperpose(EEG, 0,1,1,1,1,1,1,1); h(LASTCOM);' ...
	     '[EEG LASTCOM] = pop_rejepoch(EEG, EEG.reject.rejglobal,1);' e_newset ]);
   
third_m = uimenu( W_MAIN, 'Label', 'Plot');
	loc_m = uimenu( third_m, 'Label', 'Channel locations'   );
       uimenu( loc_m, 'Label', 'By name'   , 'CallBack'  , [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''', ''''electrodes'''', ''''labelpoint'''');'']; eval(LASTCOM);' e_hist]);
	   uimenu( loc_m, 'Label', 'By number'   , 'CallBack', [ checkplot 'LASTCOM = [''figure; topoplot([],EEG.chanlocs, ''''style'''', ''''blank'''', ''''electrodes'''', ''''numpoint'''');'']; eval(LASTCOM);' e_hist]);
    uimenu( third_m, 'Label', 'Channel data (scroll)'        , 'CallBack', [ check          'LASTCOM = pop_eegplot(EEG, 1, 0, 0);' e_hist], 'Separator', 'on');
	uimenu( third_m, 'Label', 'Channel spectra and maps' , 'CallBack', [ checkplot      'LASTCOM = pop_spectopo(EEG, 1);' e_hist], 'enable', fastif(exist('resample'), 'on', 'off'));
	uimenu( third_m, 'Label', 'Channel properties'       , 'CallBack', [ checkplot   'LASTCOM = pop_prop(EEG,1);' e_hist]);
	uimenu( third_m, 'Label', 'Channel ERP image'        , 'CallBack', [ checkepoch  'LASTCOM = pop_erpimage(EEG, 1, h(''find'',''pop_erpimage(EEG,1''));' e_hist]);
	ERP_m = uimenu( third_m, 'Label', 'Channel ERPs');
		uimenu( ERP_m, 'Label', 'With scalp maps'     , 'CallBack', [ checkepochplot 'LASTCOM = pop_timtopo(EEG);' e_hist]);
		uimenu( ERP_m, 'Label', 'In scalp array'     , 'CallBack', [ checkplot      'LASTCOM = pop_plottopo(EEG);' e_hist]);
		uimenu( ERP_m, 'Label', 'In rect. array'     , 'CallBack', [ checkepoch     '[tmpeeg LASTCOM] = pop_plotdata(EEG, 1); clear tmpeeg;' e_hist]);
	topo_m = uimenu( third_m, 'Label', 'ERP map series');
		uimenu( topo_m, 'Label', 'In 2-D'     , 'CallBack', [ checkplot      'LASTCOM = pop_topoplot(EEG, 1);' e_hist]);
		uimenu( topo_m, 'Label', 'In 3-D'     , 'CallBack', [ checkplot      '[EEG LASTCOM] = pop_headplot(EEG, 1);' e_store]);
	uimenu( third_m, 'Label', 'Compare ERPs'             , 'CallBack', [ checkepoch     'LASTCOM = pop_compareerps(ALLEEG);' e_hist]);

    uimenu( third_m, 'Label', 'Component activations (scroll)', 'CallBack', [ checkica  '[LASTCOM] = pop_eegplot(EEG, 0, 0, 0);' e_hist],'Separator', 'on');
	uimenu( third_m, 'Label', 'Component spectra and maps' , 'CallBack', [ checkicaplot   'LASTCOM = pop_spectopo(EEG, 0);' e_hist], 'enable', fastif(exist('resample'), 'on', 'off'));
	topoica_m = uimenu( third_m, 'Label', 'Component maps');
		uimenu( topoica_m, 'Label', 'In 2-D', 'CallBack', [ checkicaplot   'LASTCOM = pop_topoplot(EEG, 0);' e_hist]);
		uimenu( topoica_m, 'Label', 'In 3-D', 'CallBack', [ checkicaplot   '[EEG LASTCOM] = pop_headplot(EEG, 0);' e_store]);
	uimenu( third_m, 'Label', 'Component properties'     , 'CallBack', [ checkicaplot   'LASTCOM = pop_prop(EEG,0);' e_hist]);
	uimenu( third_m, 'Label', 'Component ERP image'      , 'CallBack', [ checkepochica  'LASTCOM = pop_erpimage(EEG, 0, h(''find'',''pop_erpimage(EEG,0''));' e_hist]);
	ERPC_m = uimenu( third_m, 'Label', 'Component ERPs');
	   uimenu( ERPC_m, 'Label', 'With component maps', 'CallBack', [ checkica 'LASTCOM = pop_envtopo(EEG);' e_hist]);
	   uimenu( ERPC_m, 'Label', 'In rectangular array'      , 'CallBack', [ checkepochica     '[tmpeeg LASTCOM] = pop_plotdata(EEG, 0); clear tmpeeg;' e_hist]);

	stat_m = uimenu( third_m, 'Label', 'Data statistics', 'Separator', 'on', 'enable', fastif(exist('kstest'), 'on', 'off'));
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
uimenu( help_m, 'Label', 'EEGLAB license', 'CallBack', 'pophelp(''license.txt'', 1);');
uimenu( help_m, 'Label', 'EEGLAB menus', 'CallBack', 'eeg_helpmenu;','separator','on');
help_subm1 = uimenu( help_m, 'Label', 'EEGLAB functions');
    uimenu( help_subm1, 'Label', 'Toolbox functions', 'CallBack', 'pophelp(''ica'');');
	uimenu( help_subm1, 'Label', 'Signal processing functions', 'callback', 'eeg_helpsigproc;');	
	uimenu( help_subm1, 'Label', 'Interactive pop_ functions', 'callback', 'eeg_helppop;');	
help_subm2 = uimenu( help_m, 'Label', 'EEGLAB advanced');
    uimenu( help_subm2, 'Label', 'Dataset structure', 'CallBack', 'pophelp(''eeg_checkset'');');
	uimenu( help_subm2, 'Label', 'Admin functions', 'callback', 'eeg_helpadmin;');	
uimenu( help_m, 'Label', 'Contact us (email)', 'CallBack', 'web(''mailto:eeglab@sccn.ucsd.edu'');');

EEGMENU = uimenu( set_m, 'Label', '------', 'Enable', 'off');
set(W_MAIN, 'userdat', { EEGUSERDAT{1} EEGMENU });
eeglab('redraw');

% REMOVED MENUS
	%uimenu( fourth_m, 'Label', 'Automatic comp. reject',  'enable', 'off', 'CallBack', '[EEG LASTCOM] = pop_rejcomp(EEG); h(LASTCOM); if ~isempty(LASTCOM), eeg_store(CURRENTSET); end;');
	%uimenu( fourth_m, 'Label', 'Reject (synthesis)' , 'Separator', 'on', 'CallBack', '[EEG LASTCOM] = pop_rejall(EEG); h(LASTCOM); if ~isempty(LASTCOM), eeg_store; end; eeglab(''redraw'');');

% --------------------
% draw the main figure
% --------------------
function eeg_mainfig;

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

h = findobj('tag', 'EEGLAB');
if ~isempty(h)
    disp('EEGLAB warning: there can be only one EEGLAB window, closing old one');
    close(h);
end;

W_MAIN = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'name', 'EEGLAB v4.0', ... 
	'numbertitle', 'off', ...
	'resize', 'off', ...
	'Position',[545.6028543307087 192.1136811023622 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) (WINY+2*BORDERINT+2*BORDEREXT) ], ...
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

geometry = { [1] [1] [1] [1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1] };
listui = { { 'style', 'text', 'string', 'Parameters of the current set', 'tag', 'win0' } { } ...
		   { 'style', 'text', 'tag', 'win1', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win2', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win3', 'string', 'Channels per frame'} ...
		   { 'style', 'text', 'tag', 'val3', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win4', 'string', 'Frames per epoch'} ...
		   { 'style', 'text', 'tag', 'val4', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win5', 'string', 'Epochs'} ...
		   { 'style', 'text', 'tag', 'val5', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win6', 'string', 'Events'} ...
		   { 'style', 'text', 'tag', 'val6', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win7', 'string', 'Sampling rate (Hz)' } ...
		   { 'style', 'text', 'tag', 'val7', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win8', 'string', 'Epoch start (sec)' } ...
		   { 'style', 'text', 'tag', 'val8', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win9', 'string', 'Epoch end (sec)' } ...
		   { 'style', 'text', 'tag', 'val9', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win10', 'string', 'Average reference' } ...
		   { 'style', 'text', 'tag', 'val10', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win11', 'string', 'Channel locations'} ...
		   { 'style', 'text', 'tag', 'val11', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win12', 'string', 'ICA weights'  } ...
		   { 'style', 'text', 'tag', 'val12', 'string', ' ' } ...
		   { 'style', 'text', 'tag', 'win13', 'string', 'Dataset size (Mb)' } ...
		   { 'style', 'text', 'tag', 'val13', 'string', ' ' } {} };
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

set(gcf, 'Position',[545.6028543307087 192.1136811023622 (WINMINX+WINMAXX+2*BORDERINT+2*BORDEREXT) (WINY+2*BORDERINT+2*BORDEREXT) ]);
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
MAX_SET = max(length( ALLEEG ), length(EEGMENU));
	
clear functions;
eeg_options;
if ~option_keepdataset
	if ~isempty(ALLEEG)
		if popask( strvcat('Remove all datasets except the present one ?', ...
						   'Otherwise go back to the memory menu to unset dataset overwrite'))
			ALLEEG = []; CURRENTSET = 0;
			h('ALLEEG = []; CURRENTSET = 0;');
		else 
			return;
		end;
	end;
	set(findobj('parent', gcf, 'label', 'Datasets'), 'enable', 'off');
	CURRENTSET = 0;
else
	if isempty(ALLEEG)
		ALLEEG = EEG;
	end;
	set(findobj('parent', gcf, 'label', 'Datasets'), 'enable', 'on');
end;

while( index <= MAX_SET)
	try
		set( EEGMENU(index), 'Label', '------');
	catch, EEGMENU(index) = uimenu( set_m, 'Label', '------', 'Enable', 'on'); end;	
	set( EEGMENU(index), 'Enable', 'on' );
	try, ALLEEG(index).data;
		if ~isempty( ALLEEG(index).data)
       		menutitle   = sprintf('Dataset %d:%s', index, ALLEEG(index).setname);
			set( EEGMENU(index), 'Label', menutitle);
			set( EEGMENU(index), 'CallBack', ['com = ''EEG = eeg_retrieve(ALLEEG, ' int2str(index) ...
					'); CURRENTSET = ' int2str(index) ';''; eval(com); h(com); eeglab(''redraw'');' ]);
			set( EEGMENU(index), 'Enable', 'on' );
		end;
	catch, end;	
	index = index+1;
end;
hh = findobj( 'parent', set_m, 'Label', '------');
set(hh, 'Enable', 'off');
EEGUSERDAT{2} = EEGMENU;
set(W_MAIN, 'userdata', EEGUSERDAT);

if option_keepdataset & (isempty(CURRENTSET) | length(ALLEEG) < CURRENTSET | CURRENTSET == 0 | isempty(ALLEEG(CURRENTSET).data))
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

if (isempty(EEG) | isempty(EEG.data)) & CURRENTSET ~= 0 & option_keepdataset
	h([ 'EEG = eeg_retrieve(ALLEEG,' int2str(CURRENTSET) '); CURRENTSET = ' int2str(CURRENTSET) ';'])
	EEG = eeg_retrieve(ALLEEG, CURRENTSET);	
end;	

% print some informations on the main figure
% ------------------------------------------
g = myguihandles(gcf);
if (exist('EEG') == 1) & isstruct(EEG) & ~isempty(EEG.data)
	h = findobj('parent', gcf, 'userdata', 'fullline');
	set(h, 'visible', 'off');
	
	if CURRENTSET == 0
		set( g.win0, 'String', sprintf('Parameters of %s dataset', ...
										  fastif(EEG.trials > 1, 'epoched', 'continuous')));	
	else  
		set( g.win0, 'String', sprintf('Parameters of %s dataset %d', ...
										  fastif(EEG.trials > 1, 'epoched', 'continuous'), CURRENTSET));	
	end;
	% set( g.win3, 'String', '');
	% set( g.win3, 'String', sprintf('Dataset name      \t\t%s\n', fastif(isempty(EEG.setname), 'none', EEG.setname)));
	set( g.win1, 'String', sprintf('Dataset name: %s\n', fastif(isempty(EEG.setname), 'none', EEG.setname)));
	fullfilename = [ EEG.filepath EEG.filename];
	if ~isempty(fullfilename)
		if length(fullfilename) > 30
			set( g.win2, 'String', sprintf('Filename: ...%s\n', fullfilename(max(1,length(fullfilename)-30):end) ));
		else
			set( g.win2, 'String', sprintf('Filename: %s\n', fullfilename));
		end;        	
	else
		set( g.win2, 'String', sprintf('Filename: none\n'));
	end;
	set( g.win3, 'String', 'Channels per frame');
	set( g.win4, 'String', 'Frames per epoch');
	set( g.win5, 'String', 'Epochs');
	set( g.win6, 'String', 'Events');
	set( g.win7, 'String', 'Sampling rate (Hz)');
	set( g.win8, 'String', 'Epoch start (sec)');
	set( g.win9, 'String', 'Epoch end (sec)');
	set( g.win10, 'String', 'Average reference');
	set( g.win11, 'String', 'Channel locations');
	set( g.win12, 'String', 'ICA weights');
	set( g.win13, 'String', 'Dataset size (Mb)');
	
	set( g.val3, 'String', int2str(fastif(isempty(EEG.data), 0, size(EEG.data,1))));
	set( g.val4, 'String', int2str(EEG.pnts));
	set( g.val5, 'String', int2str(EEG.trials));
	set( g.val6, 'String', fastif(isempty(EEG.event), 'none', int2str(length(EEG.event))));
	set( g.val7, 'String', int2str( round(EEG.srate)) );
	if round(EEG.xmin) == EEG.xmin & round(EEG.xmax) == EEG.xmax
		set( g.val8, 'String', sprintf('%d\n', EEG.xmin));
		set( g.val9, 'String', sprintf('%d\n', EEG.xmax));
	else 
		set( g.val8, 'String', sprintf('%6.3f\n', EEG.xmin));
		set( g.val9, 'String', sprintf('%6.3f\n', EEG.xmax));
	end;
	%set( g.val8, 'String', sprintf('%6.3f ±%1.3f\n', EEG.xmin+0.5/EEG.srate,0.5/EEG.srate));
	%set( g.val9, 'String', sprintf('%6.3f ±%1.3f\n', EEG.xmax+0.5/EEG.srate,0.5/EEG.srate));
	set( g.val10, 'String', EEG.averef);
	set( g.val11, 'String', fastif(isempty(EEG.chanlocs), 'No', 'Yes'));
	set( g.val12, 'String', fastif(isempty(EEG.icasphere), 'No', 'Yes'));
	tmp = whos('EEG');
	set( g.val13, 'String', num2str(round(tmp.bytes/1E6*10)/10));
else
	h = findobj('parent', gcf, 'userdata', 'fullline');
	set(h, 'visible', 'on');
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
	h = findobj('parent', gcf);
	for index = 1:length(h)
		if ~isempty(get(h(index), 'tag'))
			g = setfield(g, get(h(index), 'tag'), h(index));
		end;
	end;
