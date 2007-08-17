% eeg_checkset()   - check the consistency of the fields of an EEG dataset 
%                    Also: See EEG dataset structure field descriptions below.
%
% Usage: >> [EEGOUT,result] = eeg_checkset(EEG); % perform all checks
%                                                  except 'makeur'
%        >> [EEGOUT,result] = eeg_checkset(EEG, 'keyword'); % perform 'keyword' check(s)
%
% Inputs:
%       EEG        - EEGLAB dataset structure or (ALLEEG) array of EEG structures
%
% Optional keywords:
%   'icaconsist'   - if EEG contains several datasets, check whether they have 
%                    the same ICA decomposition
%   'epochconsist' - if EEG contains several datasets, check whether they have
%                    identical epoch lengths and time limits.
%   'chanconsist'  - if EEG contains several datasets, check whether they have 
%                    the same number of channela and channel labels.
%   'data'         - check whether EEG contains data (EEG.data)
%   'loaddata'     - load data array (if necessary)
%   'savedata'     - save data array (if necessary - see EEG.saved below)
%   'contdata'     - check whether EEG contains continuous data
%   'epoch'        - check whether EEG contains epoched or continuous data
%   'ica'          - check whether EEG contains an ICA decomposition
%   'besa'         - check whether EEG contains component dipole locations
%   'event'        - check whether EEG contains an event array
%   'makeur'       - remake the EEG.urevent structure
%   'checkur'      - check whether the EEG.urevent structure is consistent 
%                    with the EEG.event structure
%   'chanlocsize'  - check the EEG.chanlocs structure length; show warning if
%                    necessary.
%   'chanlocs_homogenous' - check whether EEG contains consistent channel
%                           information; if not, correct it.
%   'eventconsistency'    - check whether EEG.event information are consistent; 
%                           remake 'epoch' field (can be time consuming).
%
% Outputs:
%       EEGOUT     - output EEGLAB dataset or dataset array
%       result     - result code:  0 = OK; 1 = error; -1 = warning
%
% ===========================================================
% The structure of an EEG dataset under EEGLAB (as of v5.03):
%
% Basic dataset information:
%   EEG.setname      - descriptive name|title for the dataset
%   EEG.filename     - filename of the dataset file on disk
%   EEG.filepath     - filepath (directory/folder) of the dataset file(s)
%   EEG.trials       - number of epochs (or trials) in the dataset. 
%                      If data are continuous, this number is 1.
%   EEG.pnts         - number of time points (or data frames) per trial (epoch).
%                      If data are continuous (trials=1), the total number 
%                      of time points (frames) in the dataset
%   EEG.nbchan       - number of channels
%   EEG.srate        - data sampling rate (in Hz)
%   EEG.xmin         - epoch start latency|time (in sec. relative to the
%                      time-locking event at time 0)
%   EEG.xmax         - epoch end latency|time (in seconds)
%   EEG.times        - vector of latencies|times in seconds (one per time point)
%   EEG.ref          - ['common'|'averef'|integer] reference channel type or number
%   EEG.history      - cell array of ascii pop-window commands that created
%                      or modified the dataset
%   EEG.comments     - comments about the nature of the dataset (edit this via 
%                      menu selection Edit > About this dataset)
%   EEG.etc          - miscellaneous (technical or temporary) dataset information
%   EEG.saved        - ['yes'|'no'] 'no' flags need to save dataset changes before exit
%
% The data:
%   EEG.data         - two-dimensional continuous data array (chans, frames)
%                      ELSE, three-dim. epoched data array (chans, frames, epochs)
%
% The channel locations sub-structures:
%   EEG.chanlocs     - structure array containing names and locations 
%                      of the channels on the scalp
%   EEG.urchanlocs   - original (ur) dataset chanlocs structure containing
%                      all channels originally collected with these data
%                      (before channel rejection)
%   EEG.chaninfo     - structure containing additional channel info
%   EEG.ref          - type of channel reference ('common'|'averef'|+/-int]
%   EEG.splinefile   - location of the spline file used by headplot() to plot 
%                      data scalp maps in 3-D
%
% The event and epoch sub-structures:    
%   EEG.event        - event structure containing times and nature of experimental 
%                      events recorded as occurring at data time points
%   EEG.urevent      - original (ur) event structure containing all experimental
%                      events recorded as occurring at the original data time points
%                      (before data rejection)
%   EEG.epoch        - epoch event information structure array (one per epoch)
%   EEG.eventdescription - cell array of strings describing event fields.
%   EEG.epochdescription - cell array of strings describing epoch fields.
%   --> See the http://sccn.ucsd.edu/eeglab/maintut/eeglabscript.html for details
%
% ICA (or other linear) data components:
%   EEG.icasphere   - sphering array returned by linear (ICA) decomposition
%   EEG.icaweights  - unmixing weights array returned by linear (ICA) decomposition
%   EEG.icawinv     - inverse (ICA) weight matrix. Columns gives the projected
%                     topographies of the components to the electrodes.
%   EEG.icaact      - ICA activations matrix (components, frames, epochs)
%                     Note: [] here means that 'compute_ica' option has bee set 
%                     to 0 under 'File > Memory options' In this case, 
%                     component activations are computed only as needed.
%   EEG.icasplinefile - location of the spline file used by headplot() to plot 
%                     component scalp maps in 3-D
%   EEG.chaninfo.icachansind  - indices of channels used in the ICA decomposition
%   EEG.dipfit      - array of structures containing component map dipole models 
%
% Variables indicating membership of the dataset in a studyset:
%   EEG.subject     - studyset subject code 
%   EEG.group       - studyset group code 
%   EEG.condition   - studyset experimental condition code 
%   EEG.session     - studyset session number
%
% Variables used for manual and semi-automatic data rejection:
%   EEG.specdata           - data spectrum for every single trial
%   EEG.specica            - data spectrum for every single trial
%   EEG.stats              - statistics used for data rejection
%       EEG.stats.kurtc    - component kurtosis values
%       EEG.stats.kurtg    - global kurtosis of components      
%       EEG.stats.kurta    - kurtosis of accepted epochs      
%       EEG.stats.kurtr    - kurtosis of rejected epochs      
%       EEG.stats.kurtd    - kurtosis of spatial distribution      
%   EEG.reject            - statistics used for data rejection
%       EEG.reject.entropy - entropy of epochs  
%       EEG.reject.entropyc  - entropy of components
%       EEG.reject.threshold - rejection thresholds 
%       EEG.reject.icareject - epochs rejected by ICA criteria
%       EEG.reject.gcompreject - rejected ICA components
%       EEG.reject.sigreject  - epochs rejected by single-channel criteria
%       EEG.reject.elecreject - epochs rejected by raw data criteria
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

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
% Revision 1.219  2007/08/17 19:32:34  arno
% spped up by a 100 for studies
%
% Revision 1.218  2007/08/16 19:05:16  arno
% new template model format
%
% Revision 1.217  2007/08/09 16:37:14  arno
% debug last changes
%
% Revision 1.216  2007/08/09 16:10:55  arno
% same
%
% Revision 1.215  2007/08/09 16:08:30  arno
% coord_transform field
%
% Revision 1.214  2007/08/09 00:46:47  arno
% scaling ICA components to RMS
%
% Revision 1.213  2007/08/01 01:07:26  arno
% check data
%
% Revision 1.212  2007/07/27 23:09:16  arno
% rounding durations
%
% Revision 1.211  2007/03/21 22:47:05  arno
% computation of ica matrix, all matrices converted to double
%
% Revision 1.210  2007/03/08 02:55:14  toby
% Documentation edit
% Bug Fix courtesy Sebastian Jentschke
%
% Revision 1.209  2007/03/06 03:13:16  arno
% dipfitdefs check
%
% Revision 1.208  2007/02/23 19:03:38  scott
% completed help msg -- added 'result' output and the fact that EEG can also be ALLEEG
%
% Revision 1.207  2007/02/20 17:02:18  scott
% reworked help message - brought description of the EEG structure fields
% up to date -- Arno, please check.
% Toby -- please check this help message against the tutorial description
%  in the data structure appendix.
% -sm
%
% Revision 1.206  2007/02/20 14:16:50  arno
% Added memory warning
%
% Revision 1.205  2007/02/02 11:51:43  arno
% fix MEG file problem
%
% Revision 1.202  2006/12/08 19:36:11  arno
% chekcing labels
%
% Revision 1.201  2006/11/06 21:18:47  arno
% remove the +Y noze direction
%
% Revision 1.200  2006/05/10 14:37:10  arno
% nothing
%
% Revision 1.199  2006/04/14 18:04:42  arno
% checking eeg
%
% Revision 1.198  2006/04/14 17:54:25  arno
% fixing number of trials too
%
% Revision 1.197  2006/04/14 17:49:34  arno
% move com
%
% Revision 1.196  2006/04/14 17:46:10  arno
% checking data length
%
% Revision 1.195  2006/04/11 18:14:37  arno
% text message
%
% Revision 1.194  2006/04/10 21:07:14  arno
% nothing
%
% Revision 1.193  2006/03/12 03:12:06  arno
% typo
%
% Revision 1.192  2006/03/12 03:10:43  arno
% testing dipfit standard files
%
% Revision 1.191  2006/03/10 21:48:52  arno
% preserve saving status
%
% Revision 1.190  2006/03/10 20:30:41  arno
% do not recompute icawinv automatically
%
% Revision 1.189  2006/03/10 19:02:22  arno
% using pop_loadset to load data
%
% Revision 1.188  2006/03/07 17:47:39  arno
% automatically recomputing EEG.icawinv
%
% Revision 1.187  2006/03/03 22:55:50  arno
% nothing
%
% Revision 1.186  2006/02/16 19:29:08  arno
% nothing
%
% Revision 1.185  2006/02/16 00:18:42  arno
% test if icachansind empty
%
% Revision 1.184  2006/02/16 00:16:30  arno
% nothing
%
% Revision 1.183  2006/02/07 19:14:55  arno
% nothing
%
% Revision 1.182  2006/02/01 18:08:34  arno
% nothing
%
% Revision 1.181  2006/01/31 20:52:16  arno
% eeglab options
%
% Revision 1.180  2006/01/26 00:23:32  arno
% nothing
%
% Revision 1.179  2006/01/19 19:54:58  arno
% coord_tranfrom default field
%
% Revision 1.178  2005/12/01 03:32:42  toby
% *** empty log message ***
%
% Revision 1.177  2005/12/01 00:33:46  arno
% making icachansind empty
%
% Revision 1.176  2005/11/30 23:10:27  arno
% icasplinefile field
%
% Revision 1.175  2005/11/21 20:34:00  arno
% rewrite orderfields
%
% Revision 1.174  2005/11/10 23:48:27  arno
% copying icachansind
%
% Revision 1.173  2005/11/10 22:38:17  arno
% adding icachansind
%
% Revision 1.172  2005/11/07 19:40:39  arno
% nothing
%
% Revision 1.171  2005/11/04 23:02:37  arno
% fixing missing fields
%
% Revision 1.170  2005/11/04 22:55:27  arno
% changing field order
%
% Revision 1.169  2005/11/04 22:25:55  arno
% saved field
%
% Revision 1.168  2005/11/04 22:21:00  arno
% same
%
% Revision 1.167  2005/11/04 22:18:31  arno
% enforce new field order
%
% Revision 1.166  2005/10/11 23:02:32  arno
% convert data to single
%
% Revision 1.165  2005/10/11 17:17:04  arno
% move checks at the beginning
%
% Revision 1.164  2005/09/27 22:00:30  arno
% fix checking ALLEEG with some datasets missing; dealing with empty event structure
%
% Revision 1.163  2005/08/04 23:41:56  arno
% fix loaddata option
%
% Revision 1.162  2005/08/04 15:31:15  arno
% error if using savedata
%
% Revision 1.161  2005/08/03 17:37:25  arno
% updating loadset
%
% Revision 1.160  2005/08/03 01:41:10  arno
% also saving if no datfile
%
% Revision 1.159  2005/08/02 16:46:32  arno
% do not return immediately after saving data
%
% Revision 1.158  2005/08/02 16:27:57  arno
% text
%
% Revision 1.157  2005/08/02 16:22:16  arno
% test if datfile is empty
%
% Revision 1.156  2005/08/02 01:55:24  arno
% debug saving dataset
%
% Revision 1.155  2005/08/01 22:40:38  arno
% save also EEG structure
%
% Revision 1.154  2005/08/01 15:46:09  arno
% loading/writing dataset
%
% Revision 1.153  2005/07/30 01:52:11  arno
% debug loaddata savedata
%
% Revision 1.152  2005/07/30 01:45:52  arno
% loading and saving .dat array
%
% Revision 1.151  2005/07/29 15:48:20  arno
% document keywords
%
% Revision 1.150  2005/07/28 18:02:17  arno
% checking multiple dataset consistency
%
% Revision 1.149  2005/05/24 16:55:34  arno
% cell2mat and mat2cell
%
% Revision 1.148  2005/03/17 18:02:06  arno
% debug conversion
%
% Revision 1.147  2005/03/17 18:00:23  arno
% convert old dipfit structure to new one
%
% Revision 1.146  2005/03/05 02:08:46  arno
% chaninfo
%
% Revision 1.145  2005/03/04 23:45:15  arno
% implementing chaninfo
%
% Revision 1.144  2005/02/16 19:49:10  hilit
% typo
%
% Revision 1.143  2005/02/16 19:47:24  hilit
% Trying to solve compatible problems between Matlab 7 and 6
%
% Revision 1.142  2004/11/17 02:08:58  arno
% save filename
%
% Revision 1.141  2004/11/17 02:07:18  arno
% debug last
%
% Revision 1.140  2004/11/17 02:06:14  arno
% debug dat format
%
% Revision 1.139  2004/11/15 22:54:31  arno
% read data from .dat file
%
% Revision 1.138  2004/11/09 02:05:11  arno
% same
%
% Revision 1.137  2004/11/09 02:02:38  arno
% debug events for Matlab 7
%
% Revision 1.136  2004/09/22 16:50:23  hilit
% changed || -> |
%
% Revision 1.135  2004/09/21 16:54:51  hilit
% change && -> &
%
% Revision 1.134  2004/09/03 16:14:04  arno
% debug event duration
%
% Revision 1.133  2004/08/25 18:21:34  arno
% debug last change
%
% Revision 1.132  2004/08/25 17:58:52  arno
% check numeric format of events
%
% Revision 1.131  2004/08/03 01:22:56  arno
% debug cellfun for Matlab 7
%
% Revision 1.130  2004/07/30 16:59:02  arno
% new version detection
%
% Revision 1.129  2004/07/26 15:53:35  arno
% debug conversion
%
% Revision 1.128  2004/07/26 15:49:52  arno
% convert to single for Matlab 7
%
% Revision 1.127  2004/06/16 21:38:30  arno
% resorting urevents
%
% Revision 1.126  2004/06/16 18:58:58  arno
% same
%
% Revision 1.125  2004/06/16 18:57:45  arno
% resort by epoch first
%
% Revision 1.124  2004/06/14 16:10:04  arno
% resoring event latencies
%
% Revision 1.123  2004/06/04 01:05:02  arno
% allowing latencies of 0.5
%
% Revision 1.122  2004/05/26 23:12:02  arno
% epoch duration in ms
%
% Revision 1.121  2004/05/24 17:24:08  arno
% assign 0 to empty durations
%
% Revision 1.120  2004/05/14 17:47:15  arno
% convert history
%
% Revision 1.119  2004/05/06 21:54:55  arno
% message text
%
% Revision 1.118  2004/02/17 20:05:04  arno
% remove ICA weights if invalid
%
% Revision 1.117  2003/12/17 23:25:39  arno
% different check for chanlocs
%
% Revision 1.116  2003/12/17 00:45:57  arno
% adding E prefix to electrodes
%
% Revision 1.115  2003/12/11 17:59:41  arno
% msg
%
% Revision 1.114  2003/12/04 22:33:54  arno
% urevent remove
%
% Revision 1.113  2003/12/04 17:45:12  arno
% debug urevent (made Matlab 5.3 crash systematically)
%
% Revision 1.112  2003/12/04 02:40:56  arno
% msg
%
% Revision 1.111  2003/12/02 22:29:40  arno
% adding urchan
%
% Revision 1.110  2003/12/02 17:07:59  arno
% debug urchanlocs
%
% Revision 1.109  2003/12/02 03:39:06  arno
% urchanlocs
%
% Revision 1.108  2003/11/18 16:42:30  scott
% text labels
%
% Revision 1.107  2003/11/05 16:20:27  arno
% homogenous -> homogeneous
%
% Revision 1.106  2003/11/04 15:47:48  scott
% warning msg edits
%
% Revision 1.105  2003/11/04 01:15:06  arno
% makeur message made clearer
%
% Revision 1.104  2003/09/22 23:45:24  arno
% debug urevent
%
% Revision 1.103  2003/07/28 15:30:31  arno
% default for EEG.ref
%
% Revision 1.102  2003/07/28 15:19:09  arno
% detect reference electrode
%
% Revision 1.101  2003/07/21 14:32:17  arno
% convert single to double precision
%
% Revision 1.100  2003/07/16 20:53:43  arno
% auto creation of urevent table
%
% Revision 1.99  2003/06/27 16:59:11  arno
% updating ur
%
% Revision 1.98  2003/06/18 22:27:53  arno
% implementing makeur
%
% Revision 1.97  2003/06/13 16:41:51  arno
% adding chanlocs homogenous check
%
% Revision 1.96  2003/06/11 21:18:04  arno
% new limits for latencies
%
% Revision 1.95  2003/02/28 17:05:23  arno
% eeg_checkset() -> eeg_checkset in warnings
%
% Revision 1.94  2003/02/28 16:57:34  arno
% typo
%
% Revision 1.93  2003/02/28 15:35:41  scott
% header edits -sm
%
% Revision 1.92  2003/02/28 15:30:00  arno
% updating warning message
%
% Revision 1.91  2003/02/26 02:18:55  arno
% debugging if file has changed of location
%
% Revision 1.90  2003/02/03 20:07:45  arno
% error if no data
%
% Revision 1.89  2003/01/24 19:32:02  arno
% debugging ICA for NaN
%
% Revision 1.88  2003/01/02 17:13:01  scott
% edit header and msgs -sm
%
% Revision 1.87  2003/01/02 16:37:33  arno
% editing message - ad & sm
%
% Revision 1.86  2002/12/24 01:34:44  arno
% debug multiple checks
%
% Revision 1.85  2002/11/15 18:40:54  arno
% adding another test if chanlocs empty
%
% Revision 1.84  2002/11/15 02:11:04  arno
% debugging for single dataset
%
% Revision 1.83  2002/11/15 01:37:42  scott
% Can not -> Cannot
%
% Revision 1.82  2002/11/13 19:57:40  arno
% checkin shrink factor
%
% Revision 1.81  2002/11/13 19:21:09  arno
% updating average reference flag
%
% Revision 1.80  2002/11/13 17:41:11  arno
% editing chanlocs warning -sm
%
% Revision 1.79  2002/11/13 17:10:09  arno
% forcing channel labels to string
%
% Revision 1.78  2002/11/12 22:51:55  arno
% adding a warning for additional reference electrode location
%
% Revision 1.77  2002/11/11 15:28:54  arno
% besa check
%
% Revision 1.76  2002/10/29 01:17:37  arno
% implementing user abord
%
% Revision 1.75  2002/10/20 21:32:08  arno
% nan activation computation debug
%
% Revision 1.74  2002/10/16 22:44:04  arno
% ica recompute for NaNs
%
% Revision 1.73  2002/10/09 00:14:13  arno
% typo last
%
% Revision 1.72  2002/10/09 00:11:41  arno
% debug read float data file
%
% Revision 1.71  2002/09/23 16:42:24  arno
% adding comments
%
% Revision 1.70  2002/09/23 16:15:00  arno
% debug floatread
%
% Revision 1.69  2002/09/23 16:08:57  arno
% check for EEG.data empty
%
% Revision 1.68  2002/09/05 00:04:05  arno
% disp-> error
%
% Revision 1.67  2002/09/04 22:13:55  luca
% adding dataset name check -arno
%
% Revision 1.66  2002/08/28 01:02:48  arno
% changing error messages to disp
%
% Revision 1.65  2002/08/22 21:21:00  arno
% typo
%
% Revision 1.64  2002/08/21 17:56:20  arno
% debug checks
%
% Revision 1.63  2002/08/21 17:46:26  arno
% more reject field checks
%
% Revision 1.62  2002/08/21 02:24:27  arno
% change message
%
% Revision 1.61  2002/08/21 02:22:54  arno
% debug
%
% Revision 1.60  2002/08/21 02:19:38  arno
% add continuous data statement
%
% Revision 1.59  2002/08/21 00:15:20  arno
% debug
%
% Revision 1.58  2002/08/19 19:46:16  arno
% for non cellfun compatibility
%
% Revision 1.57  2002/08/14 02:01:08  arno
% debugging epoch info
%
% Revision 1.56  2002/08/12 18:53:14  arno
% errordlg2
%
% Revision 1.55  2002/08/12 18:51:58  arno
% errordlg2
%
% Revision 1.54  2002/08/12 18:39:23  arno
% questdlg2
%
% Revision 1.53  2002/08/12 00:16:46  arno
% same
%
% Revision 1.52  2002/08/12 00:13:44  arno
% same
%
% Revision 1.51  2002/08/12 00:05:44  arno
% changing manual color
%
% Revision 1.50  2002/08/09 00:58:38  arno
% text
%
% Revision 1.49  2002/08/09 00:39:18  arno
% debugging epoch
%
% Revision 1.48  2002/08/08 21:55:59  arno
% adding epoch creation
%
% Revision 1.47  2002/08/08 21:08:34  arno
% *** empty log message ***
%
% Revision 1.46  2002/07/30 22:05:24  arno
% adding disprej field
%
% Revision 1.45  2002/07/30 17:53:22  arno
% adding color for rejection
%
% Revision 1.44  2002/07/29 16:42:00  arno
% debugging
%
% Revision 1.43  2002/07/27 00:08:03  arno
% debugging
%
% Revision 1.42  2002/07/26 18:06:41  arno
% add warning
%
% Revision 1.41  2002/07/25 17:36:34  arno
% debugging
%
% Revision 1.40  2002/07/25 17:34:48  arno
% adding message when removing ICA array
%
% Revision 1.39  2002/07/25 17:12:55  arno
% debugging gcompreject
%
% Revision 1.38  2002/07/24 18:40:07  arno
% checking empty values in epochs
%
% Revision 1.37  2002/07/23 23:51:25  arno
% removing error
%
% Revision 1.36  2002/07/23 22:23:41  arno
% removing warning if icaact=[]
%
% Revision 1.35  2002/07/23 21:29:07  arno
% empty icaact
%
% Revision 1.34  2002/07/23 00:13:29  arno
% adding read float feature
%
% Revision 1.33  2002/06/25 13:40:00  arno
% adding EEG.times
%
% Revision 1.32  2002/06/25 02:31:06  arno
% gcompreject initialized to zeros for all components
%
% Revision 1.31  2002/06/25 00:45:59  arno
% removing epoch info unofor;isation
% ,
%
% Revision 1.30  2002/05/04 01:47:29  arno
% same
%
% Revision 1.29  2002/05/04 01:46:42  arno
% still correctin eventconsistency
%
% Revision 1.28  2002/05/04 01:45:13  arno
% typo
%
% Revision 1.27  2002/05/04 01:44:26  arno
% correcting typo
%
% Revision 1.26  2002/05/04 01:42:38  arno
% cellfun bug correction
%
% Revision 1.25  2002/05/03 01:33:58  luca
% debuging new event check
%
% Revision 1.24  2002/05/02 22:16:54  arno
% speeding up event checks
%
% Revision 1.23  2002/05/01 18:58:08  luca
% same
%
% Revision 1.22  2002/05/01 18:57:48  luca
% same
%
% Revision 1.21  2002/05/01 18:57:05  luca
% icaact reshape problem
%
% Revision 1.20  2002/04/30 15:21:44  scott
% editted help msg -sm
%
% Revision 1.19  2002/04/20 18:45:02  arno
% editing error message
%
% Revision 1.18  2002/04/18 16:14:25  scott
% EEG.ref = 'No' by default -sm
%
% Revision 1.17  2002/04/18 02:37:54  scott
% [same] -sm
%
% Revision 1.16  2002/04/18 02:35:42  scott
% [same] -sm
%
% Revision 1.15  2002/04/18 02:34:08  scott
% improved "empty dataset" msg -sm
%
% Revision 1.14  2002/04/11 23:34:49  arno
% adding event check in event consistency
%
% Revision 1.13  2002/04/11 18:21:43  arno
% add furhter check for EEG.ref
%
% Revision 1.12  2002/04/11 18:08:47  arno
% adding average reference variable check
% ,
%
% Revision 1.11  2002/04/10 00:42:55  arno
% reprograming eventconsistency for higher speed
%
% Revision 1.10  2002/04/10 00:08:28  arno
% debuging event consistency
%
% Revision 1.9  2002/04/10 00:03:03  arno
% debuging event consistency
%
% Revision 1.8  2002/04/09 21:02:05  arno
% adding further check to eventconsistency
%
% Revision 1.7  2002/04/09 20:11:33  arno
% eventdesciption advanced checking
%
% Revision 1.6  2002/04/09 02:38:12  arno
% bedugging epoch event consistency check
%
% Revision 1.5  2002/04/09 01:52:57  arno
% adding check of event latencies
%
% Revision 1.4  2002/04/08 21:52:51  arno
% checking event description consistency
%
% Revision 1.3  2002/04/08 20:49:26  arno
% add check for the comments field
%
% Revision 1.2  2002/04/08 02:13:09  scott
% improved wording of messages to user -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 01-26-02 chandeg events and trial condition format -ad
% 01-27-02 debug when trial condition is empty -ad
% 02-15-02 remove icawinv recompute for pop_epoch -ad & ja
% 02-16-02 remove last modification and test icawinv separatelly -ad
% 02-16-02 empty event and epoch check -ad
% 03-07-02 add the eeglab options -ad
% 03-07-02 corrected typos and rate/point calculation -ad & ja
% 03-15-02 add channel location reading & checking -ad
% 03-15-02 add checking of ICA and epochs with pop_up windows -ad
% 03-27-02 recorrected rate/point calculation -ad & sm

function [EEG, res] = eeg_checkset( EEG, varargin );
msg = '';
res = 0; % 0 = OK, 1 = error, -1=warning
com = sprintf('EEG = eeg_checkset( EEG );');

if nargin < 1
    help eeg_checkset;
    return;
end;

if isempty(EEG), return; end;
if ~isfield(EEG, 'data'), return; end;

% checking multiple datasets
% --------------------------
if length(EEG) > 1
 
    if nargin > 1
        switch varargin{1}
            case 'epochconsist', % test epoch consistency
                                 % ----------------------
            res = 'no';
            datasettype = unique( [ EEG.trials ] );
            if datasettype(1) == 1 & length(datasettype) == 1, return; % continuous data
            elseif datasettype(1) == 1,                        return; % continuous and epoch data
            end;
            
            allpnts = unique( [ EEG.pnts ] );
            allxmin = unique( [ EEG.xmin ] );
            if length(allpnts) == 1 & length(allxmin) == 1, res = 'yes'; end;
            return;

            case 'chanconsist'  % test channel number and name consistency
                                % ----------------------------------------
             res = 'yes';
             chanlen    = unique( [ EEG.nbchan ] );             
             anyempty    = unique( cellfun( 'isempty', { EEG.chanlocs }) );
             if length(chanlen) == 1 & all(anyempty == 0)
                 channame1 = { EEG(1).chanlocs.labels };
                 for i = 2:length(EEG)
                     channame2 = { EEG(i).chanlocs.labels };
                     if length(intersect(channame1, channame2)) ~= length(channame1), res = 'no'; end;
                 end;
             else res = 'no';
             end;
             return;
             
         case 'icaconsist'  % test ICA decomposition consistency
                            % ----------------------------------
          res = 'yes';
          anyempty    = unique( cellfun( 'isempty', { EEG.icaweights }) );        
          if length(anyempty) == 1 & anyempty(1) == 0
              ica1 = EEG(1).icawinv;
              for i = 2:length(EEG)
                  if ~isequal(EEG(1).icawinv, EEG(i).icawinv)
                      res = 'no';
                  end;
              end;
          else res = 'no';
          end;
          return;
             
        end; 
    end;
    
end;

% reading these option take time because
% of disk access
% --------------
if exist('dipfitdefs'), dipfitdefs; end;
eeglab_options;

% standard checking
% -----------------
ALLEEG = EEG;
for inddataset = 1:length(ALLEEG)
    if ~isempty(ALLEEG(inddataset).data)

        EEG = ALLEEG(inddataset);
        
        if isempty(EEG.data)
            errordlg2(strvcat('Error: no data'), 'Error');
            error('eeg_checkset error: no data'); return;
        end;              

        if ~isempty( varargin)
            if isempty(EEG.data)
                errordlg2('Empty dataset -> File / Import data or File / Load existing dataset', 'Error');
                error('eeg_checkset error: empty dataset'); return;
            end;    
        end;

        % additional checks
        % -----------------
        res = -1; % error code
        if ~isempty( varargin)
            for index = 1:length( varargin )
                switch varargin{ index }
                 case 'data',; % already done at the top 
                 case 'contdata',;
                  if EEG.trials > 1
                      errordlg2(strvcat('Error: function only works on continuous data'), 'Error');
                      return;
                  end;
                 case 'ica', 
                  if isempty(EEG.icaweights)
                      errordlg2(strvcat('Error: no ICA decomposition. use menu "Tools > Run ICA" first.'), 'Error');
                      return;
                  end;
                 case 'epoch', 
                  if EEG.trials == 1
                      errordlg2(strvcat('Extract epochs before running that function', 'Use Tools > Extract epochs'), 'Error');
                      return
                  end;
                 case 'besa', 
                  if ~isfield(EEG, 'sources')
                      errordlg2(strvcat('No dipole information', '1) Export component maps: Tools > Localize ... BESA > Export ...' ...
                                        , '2) Run BESA to localize the equivalent dipoles', ...
                                        '3) Import the BESA dipoles: Tools > Localize ... BESA > Import ...'), 'Error');
                      return
                  end;
                 case 'event', 
                  if isempty(EEG.event)
                      errordlg2(strvcat('Cannot process if no events. Add events.', ...
                                        'Use "Edit > event fields" to create event fields.', ...
                                        'Or use "File > Import event info" or "File > Import epoch info"'), 'Error');
                      return;
                  end;
                 case 'chanloc', 
                  if isempty(EEG.chanlocs) | ~isfield(EEG.chanlocs, 'theta') | ...
                          all(cellfun('isempty', { EEG.chanlocs.theta }))
                      errordlg2( strvcat('Cannot process dataset without channel location information.', ...
                                         'Enter the filename via "Edit > Edit dataset info".', ...
                                         'For file format, enter ''>> help readlocs'' from the command line.'), 'Error');
                      return;
                  end;
                 case 'chanlocs_homogeneous', 
                  if isempty(EEG.chanlocs) | ~isfield(EEG.chanlocs, 'theta') | ...
                          all(cellfun('isempty', { EEG.chanlocs.theta }))
                      errordlg2( strvcat('Cannot process without a channel location information.', ...
                                         'Enter the filename via "Edit > Edit dataset info".', ...
                                         'For file format, enter ''>> help readlocs'' from the command line.'), 'Error');
                      return;
                  end;
                  if ~isfield(EEG.chanlocs, 'X') | isempty(EEG.chanlocs(1).X)
                      EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
                      res = [ inputname(1) ' = eeg_checkset('  inputname(1) ', ''chanlocs_homogeneous'' ); ' ];
                  end;
                 case 'chanlocsize', 
                  if ~isempty(EEG.chanlocs)
                      if length(EEG.chanlocs) > EEG.nbchan
                          questdlg2(strvcat('Warning: there is one more electrode location than', ...
                                            'data channels. EEGLAB will consider the last electrode to be the', ...
                                            'common reference channel. If this is not the case, remove the', ...
                                            'extra channel'), 'Warning', 'Ok', 'Ok');
                      end;    
                  end;
                 case 'makeur', 
                  if ~isempty(EEG.event)
                      if isfield(EEG.event, 'urevent'), 
                          EEG.event = rmfield(EEG.event, 'urevent');
                          disp('eeg_checkset note: re-creating the original event table (EEG.urevent)');
                      else
                          disp('eeg_checkset note: creating the original event table (EEG.urevent)');
                      end;
                      EEG.urevent = EEG.event;
                      for index = 1:length(EEG.event)
                          EEG.event(index).urevent = index;
                      end;
                  end;
                 case 'checkur', 
                  if ~isempty(EEG.event)
                      if isfield(EEG.event, 'urevent') & ~isempty(EEG.urevent)
                          urlatencies = [ EEG.urevent.latency ];
                          [newlat tmpind] = sort(urlatencies);
                          if ~isequal(newlat, urlatencies)
                              EEG.urevent   = EEG.urevent(tmpind);
                              [tmp tmpind2] = sort(tmpind);
                              for index = 1:length(EEG.event)
                                  EEG.event(index).urevent = tmpind2(EEG.event(index).urevent);
                              end;
                          end;
                      end;
                  end;
                 case 'eventconsistency',          
                  [EEG res] = eeg_checkset(EEG);
                  if isempty(EEG.event), return; end;
                  
                  % remove the events which latency are out of boundary
                  % ---------------------------------------------------
                  if isfield(EEG.event, 'latency')
                      try, alllatencies = [ EEG.event.latency ];
                      catch, error('Checkset: error empty latency entry for new events added by user');
                      end;
                      I1 = find(alllatencies < 0.5);
                      I2 = find(alllatencies > EEG.pnts*EEG.trials);
                      if (length(I1) + length(I2)) > 0 
                          fprintf('eeg_checkset warning: %d/%d events had out-of-bounds latencies and were removed\n', ...
                                  length(I1) + length(I2), length(EEG.event));
                          EEG.event(union(I1, I2)) = [];
                      end;
                  end;
                  
                  % save information for non latency fields updates
                  % -----------------------------------------------
                  difffield = [];
                  if ~isempty(EEG.event) & isfield(EEG.event, 'epoch')
                      % remove fields with empty epochs
                      % -------------------------------
                      removeevent = [];
                      try, allepochs = [ EEG.event.epoch ];
                          removeevent = find( allepochs < 1 | allepochs > EEG.trials);
                          if ~isempty(removeevent)
                              disp([ 'eeg_checkset warning: ' int2str(length(removeevent)) ' event had invalid epoch numbers and were removed']);
                          end;
                      catch, 
                          for indexevent = 1:length(EEG.event)
                              if isempty( EEG.event(indexevent).epoch ) | ~isnumeric(EEG.event(indexevent).epoch) ...
                                      | EEG.event(indexevent).epoch < 1 | EEG.event(indexevent).epoch > EEG.trials
                                  removeevent = [removeevent indexevent];
                                  disp([ 'eeg_checkset warning: event ' int2str(indexevent) ' has an invalid epoch number: removed']);
                              end;
                          end;
                      end;
                      EEG.event(removeevent) = [];
                      allepochs = [ EEG.event.epoch ];
                      
                      % uniformize fields content for the different epochs
                      % --------------------------------------------------
                      % THIS WAS REMOVED SINCE SOME FIELDS ARE ASSOCIATED WITH THE EVENT AND NOT WITH THE EPOCH
                      % I PUT IT BACK, BUT IT DOES NOT ERASE NON-EMPTY VALUES
                      difffield = setdiff( fieldnames(EEG.event), { 'latency' 'epoch' 'type' });
                      for index = 1:length(difffield)
                          eval(['allvalues = { EEG.event.' difffield{index} ' };']);
                          try,   eval(['valempt = cellfun(''isempty'', allvalues);']);
                          catch, valempt = mycellfun('isempty', allvalues);
                          end;
                          arraytmpinfo = cell(1,EEG.trials);
                          
                          % spetial case of duration
                          % ------------------------
                          if strcmp( difffield{index}, 'duration')
                              if any(valempt)
                                  fprintf(['eeg_checkset: found empty values for field ''' difffield{index} ...
                                           ''' (filling with 0)\n']);
                              end;
                              for indexevent = find(valempt)
                                  EEG.event(indexevent).duration = 0;
                              end;
                          else
                              
                              % get the field content
                              % ---------------------
                              for indexevent = 1:length(EEG.event)
                                  if ~valempt(indexevent)
                                      arraytmpinfo{allepochs(indexevent)} = allvalues{indexevent};
                                  end;
                              end;
                              % uniformize content for all epochs
                              % ---------------------------------
                              for indexevent = 1:length(EEG.event)
                                  if valempt(indexevent)
                                      EEG.event = setfield( EEG.event, { indexevent }, difffield{index}, ...
                                                                       arraytmpinfo{allepochs(indexevent)});
                                  end;
                              end;
                              if any(valempt)
                                  fprintf(['eeg_checkset: found empty values for field ''' difffield{index} '''\n']);
                                  fprintf(['              filling with values of other events in the same epochs\n']);
                              end;
                          end;
                      end;
                  end;
                  
                  % uniformize fields (str or int) if necessary
                  % -------------------------------------------
                  allfields = fieldnames(EEG.event);
                  for index = 1:length(allfields)
                      eval(['allvalues = { EEG.event.' allfields{index} ' };']);
                      try,   eval(['valreal = ~cellfun(''isclass'', allvalues, ''char'');']);
                      catch, valreal = mycellfun('isclass', allvalues, 'double');
                      end;
                      
                      format = 'ok';
                      if ~all(valreal) % all valreal ok
                          format = 'str';
                          if all(valreal == 0) % all valreal=0 ok
                              format = 'ok';
                          end;
                      end;
                      if strcmp(format, 'str')
                          fprintf('eeg_checkset note: value format of event field ''%s'' made uniform\n', allfields{index});
                          % get the field content
                          % ---------------------
                          for indexevent = 1:length(EEG.event)
                              if valreal(indexevent)
                                  EEG.event = setfield(EEG.event, { indexevent }, allfields{index}, num2str(allvalues{indexevent}) );
                              end;
                          end;
                      end;
                  end;

                  % check that numeric format is double (Matlab 7)
                  % -----------------------------------
                  allf = fieldnames(EEG.event);
                  if ~isempty(EEG.event)
                      for index = 1:length(allfields)
                          clear tmpval; tmpval = getfield(EEG.event,{ 1 },allf{index});
                          if isnumeric(tmpval) & ~isa(tmpval, 'double')
                              for indexevent = 1:length(EEG.event)
                                  tmpval  =   getfield(EEG.event, { indexevent }, allf{index} );
                                  EEG.event = setfield(EEG.event, { indexevent }, allf{index}, double(tmpval));
                              end;
                          end;
                      end;
                  end;
                  
                  % check duration field, replace empty by 0
                  % ----------------------------------------
                  if isfield(EEG.event, 'duration')
                      try,   valempt = cellfun('isempty'  , { EEG.event.duration });
                      catch, valempt = mycellfun('isempty', { EEG.event.duration });
                      end;
                      if any(valempt),
                          for index = find(valempt)
                              EEG.event(index).duration = 0;
                          end;
                      end;
                  end;

                  % resort events
                  % -------------
                  if isfield(EEG.event, 'latency')
                      try, 
                          if isfield(EEG.event, 'epoch')
                              TMPEEG = pop_editeventvals(EEG, 'sort', { 'epoch' 0 'latency' 0 });
                          else
                              TMPEEG = pop_editeventvals(EEG, 'sort', { 'latency' 0 });
                          end;
                          if ~isequal(TMPEEG.event, EEG.event)
                              EEG = TMPEEG;
                              disp('Event resorted by increasing latencies. Some event indices have changed.');
                          end;
                      catch,
                          disp('eeg_checkset: problem when attempting to resort event latencies.');     
                      end;
                  end;
                  
                  % build epoch structure
                  % ---------------------
                  try,
                      if EEG.trials > 1 & ~isempty(EEG.event)
                          maxlen = 0;
                          EEG.epoch = [];
                          EEG.epoch(1).event = [];    
                          EEG.epoch(EEG.trials).event = [];    
                          for index = 1:length(EEG.event)
                              currentepoch = EEG.event(index).epoch;
                              if currentepoch <= length(EEG.epoch)
                                  EEG.epoch(currentepoch).event = [ EEG.epoch(currentepoch).event index ];
                              else
                                  EEG.epoch(currentepoch).event = [ index ];
                              end;
                              maxlen = max(length(EEG.epoch(currentepoch).event), maxlen);
                          end;
                          
                          % copy event information into the epoch array
                          % -------------------------------------------
                          eventfields = fieldnames(EEG.event);
                          eventfields = setdiff(eventfields, 'epoch');
                          for fieldnum = 1:length(eventfields)
                              eval( ['allfieldvals = { EEG.event.' eventfields{fieldnum} '};'] );
                              for trial = 1:EEG.trials
                                  valfield = allfieldvals( EEG.epoch(trial).event );
                                  if ~isempty(valfield) & strcmp(eventfields{fieldnum}, 'latency')
                                      valfield = eeg_point2lat([ valfield{:} ] ,trial,EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
                                      valfield = round( valfield * 10^8 )/10^8;
                                      valfield = mattocell(valfield);
                                  end;
                                  if ~isempty(valfield) & strcmp(eventfields{fieldnum}, 'duration')
                                      valfield = [ valfield{:} ]/EEG.srate*1000;
                                      valfield = mattocell(valfield);
                                  end;
                                  if ~isempty(valfield)
                                      if maxlen == 1, EEG.epoch = setfield(EEG.epoch, { trial }, ['event' eventfields{fieldnum}], valfield{1});
                                      else            EEG.epoch = setfield(EEG.epoch, { trial }, ['event' eventfields{fieldnum}], valfield);
                                      end;
                                  end;
                              end;
                          end;    
                      end;
                  catch, errordlg2(['Warning: minor problem encountered when generating' 10 ...
                                    'the EEG.epoch structure (used only in user scripts)']); return;
                  end;
                 case { 'loaddata', 'savedata'},;
                 otherwise, error('eeg_checkset: unknown option');
                end;        
            end;
        end;            

        res = [];

        % check name consistency
        % ----------------------
        if ~isempty(EEG.setname)
            if ~isstr(EEG.setname)
                EEG.setname = '';
            else
                if size(EEG.setname,1) > 1
                    disp('eeg_checkset warning: invalid dataset name, removed');
                    EEG.setname = '';
                end;
            end;
        else
            EEG.setname = '';
        end;    

        % checking history and convert if necessary
        % -----------------------------------------
        if isfield(EEG, 'history') & size(EEG.history,1) > 1
            allcoms = cellstr(EEG.history);
            EEG.history = deblank(allcoms{1});
            for index = 2:length(allcoms)
                EEG.history = [ EEG.history 10 deblank(allcoms{index}) ];
            end;
        end;

        % read data if necessary
        % ----------------------
        if isstr(EEG.data) & nargin > 1
            if strcmpi(varargin{1}, 'loaddata')

                if strcmpi(EEG.data, 'in set file')
                    filename = fullfile(EEG.filepath, EEG.filename);
                    EEG = pop_loadset(filename);
                else
                    % opening data file
                    % -----------------
                    filename = fullfile(EEG.filepath, EEG.data);
                    fid = fopen( filename, 'r', 'ieee-le'); %little endian (see also pop_saveset)
                    if fid == -1
                        error( ['file ' filename ' not found. If you have renamed/moved' 10 ...
                                'the .set file, you must also rename/move the associated data file.' ]);
                    else 
                        fprintf('Reading float file ''%s''...\n', filename);
                    end;
                    
                    % old format = .fdt; new format = .dat (transposed)
                    % -------------------------------------------------
                    datformat = 0;
                    if length(filename) > 3
                        if strcmpi(filename(end-2:end), 'dat')
                            datformat = 1;
                        end;
                    end;
                    EEG.datfile = EEG.data;
                    if datformat
                        EEG.data = fread(fid, [EEG.trials*EEG.pnts EEG.nbchan], 'float32')';
                    else
                        EEG.data = fread(fid, [EEG.nbchan Inf], 'float32');
                    end;
                    fclose(fid);
                end;
            end;
        end;

        % save data if necessary
        % ----------------------
        if nargin > 1
            
            % datfile available?
            % ------------------
            datfile = 0;
            if isfield(EEG, 'datfile')
                if ~isempty(EEG.datfile)
                    datfile = 1;
                end;
            end;
            
            % save data
            % ---------
            if strcmpi(varargin{1}, 'savedata') & option_storedisk
                error('eeg_checkset: cannot call savedata any more');
                
                if ~isstr(EEG.data) % not already saved
                    disp('Writing previous dataset to disk...');
                    
                    if datfile
                        tmpdata = reshape(EEG.data, EEG.nbchan,  EEG.pnts*EEG.trials);
                        floatwrite( tmpdata', fullfile(EEG.filepath, EEG.datfile), 'ieee-le');
                        EEG.data   = EEG.datfile;
                    end;            
                    EEG.icaact = [];
                    
                    % saving dataset
                    % --------------
                    filename = fullfile(EEG(1).filepath, EEG(1).filename);
                    if ~isstr(EEG.data), EEG.data = single(EEG.data); end;
                    v = version;
                    if str2num(v(1)) >= 7, save( filename, '-v6', '-mat', 'EEG'); % Matlab 7
                    else                   save( filename, '-mat', 'EEG');
                    end;
                    if ~isstr(EEG.data), EEG.data = 'in set file'; end;
                    
                    res = sprintf('%s = eeg_checkset( %s, ''savedata'');', inputname(1), inputname(1));
                end;
            end;
        end;

        % numerical format
        % ----------------
        if isnumeric(EEG.data)
            v = version;
            EEG.icawinv    = double(EEG.icawinv); % required for dipole fitting, otherwise it crashes
            EEG.icaweights = double(EEG.icaweights);
            EEG.icasphere  = double(EEG.icasphere);
            if ~isempty(findstr(v, 'R11')) | ~isempty(findstr(v, 'R12')) | ~isempty(findstr(v, 'R13'))
                EEG.data       = double(EEG.data);
                EEG.icaact     = double(EEG.icaact);
            else
                try,
                    EEG.data       = single(EEG.data);
                    EEG.icaact     = single(EEG.icaact);
                catch,
                    disp('WARNING: EEGLAB ran out of memory while converting dataset to single precision.');
                    disp('         Save dataset (preferably saving data to a separate file; see File > Memory options).'); 
                    disp('         Then reload it.');
                end;
            end;
        end;

        % verify the type of the variables
        % --------------------------------
        % data dimensions -------------------------
        if isnumeric(EEG.data)
            if size(EEG.data,1) ~= EEG.nbchan
                disp( [ 'eeg_checkset warning: number of columns in data (' int2str(size(EEG.data,1)) ...
                        ') does not match the number of channels (' int2str(EEG.nbchan) '): corrected' ]); 
                res = com;
                EEG.nbchan = size(EEG.data,1);
            end;    

            if (ndims(EEG.data)) < 3 & (EEG.pnts > 1)
                if mod(size(EEG.data,2), EEG.pnts) ~= 0
                    if popask( [ 'eeg_checkset error: the number of frames does not divide the number of columns in the data.'  10 ...
                                 'Should EEGLAB attempt to abort operation ?' 10 '(press Cancel to fix the problem from the command line)']) 
                        error('eeg_checkset error: user abort');
                        %res = com;
                        %EEG.pnts = size(EEG.data,2);
                        %EEG = eeg_checkset(EEG);
                        %return;
                    else
                        res = com;
                        return;
                        %error( 'eeg_checkset error: number of points does not divide the number of columns in data');
                    end;        
                else
                    if EEG.trials > 1
                        disp( 'eeg_checkset note: data array made 3-D'); 
                        res = com;
                    end;    
                    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, size(EEG.data,2)/EEG.pnts);         
                end;    
            end;

            % size of data -----------
            if size(EEG.data,3) ~= EEG.trials 
                disp( ['eeg_checkset warning: 3rd dimension size of data (' int2str(size(EEG.data,3)) ...
                       ') does not match the number of epochs (' int2str(EEG.trials) '), corrected' ]); 
                res = com;
                EEG.trials = size(EEG.data,3);
            end;    
            if size(EEG.data,2) ~= EEG.pnts 
                disp( [ 'eeg_checkset warning: number of columns in data (' int2str(size(EEG.data,2)) ...
                        ') does not match the number of points (' int2str(EEG.pnts) '): corrected' ]); 
                res = com;
                EEG.pnts = size(EEG.data,2);
            end;    
        end;

        % parameters consistency -------------------------
        if     round(EEG.srate*(EEG.xmax-EEG.xmin)+1) ~= EEG.pnts          
            fprintf( 'eeg_checkset note: upper time limit (xmax) adjusted so (xmax-xmin)*srate+1 = number of frames\n'); 
            if EEG.srate == 0
                EEG.srate = 1;
            end;
            EEG.xmax = (EEG.pnts-1)/EEG.srate+EEG.xmin;
            res = com;
        end;
        
        % deal with event arrays
        % ----------------------
        if ~isfield(EEG, 'event'), EEG.event = []; res = com; end;
        if ~isempty(EEG.event)
            if EEG.trials > 1 & ~isfield(EEG.event, 'epoch')
                if popask( [ 'eeg_checkset error: the event info structure does not contain an ''epoch'' field.'  ...
                             'Should EEGLAB attempt to abort operation ?' 10 '(press Cancel to fix the problem from the commandline)']) 
                    error('eeg_checkset error(): user abort');
                    %res = com;
                    %EEG.event = [];
                    %EEG = eeg_checkset(EEG);
                    %return;
                else 
                    res = com;
                    return;
                    %error('eeg_checkset error: no epoch field in event structure');
                end;
            end;
        else
            EEG.event = [];
        end;
        if isempty(EEG.event)
            EEG.eventdescription = {};
        end;
        if ~isfield(EEG, 'eventdescription') | ~iscell(EEG.eventdescription)
            EEG.eventdescription = cell(1, length(fieldnames(EEG.event)));
            res = com; 
        else 
            if ~isempty(EEG.event)
                if length(EEG.eventdescription) > length( fieldnames(EEG.event))
                    EEG.eventdescription = EEG.eventdescription(1:length( fieldnames(EEG.event)));
                elseif length(EEG.eventdescription) < length( fieldnames(EEG.event))
                    EEG.eventdescription(end+1:length( fieldnames(EEG.event))) = {''};
                end;
            end;
        end;
        % create urevent if continuous data
        % ---------------------------------
        %if ~isempty(EEG.event) & ~isfield(EEG, 'urevent')
        %    EEG.urevent = EEG.event;
        %   disp('eeg_checkset note: creating the original event table (EEG.urevent)');
        %    for index = 1:length(EEG.event)
        %        EEG.event(index).urevent = index;
        %    end;
        %end;
        if isfield(EEG, 'urevent') & isfield(EEG.urevent, 'urevent')
            EEG.urevent = rmfield(EEG.urevent, 'urevent');
        end;
        
        % deal with epoch arrays
        % ----------------------
        if ~isfield(EEG, 'epoch'), EEG.epoch = []; res = com; end;
        if ~isfield(EEG, 'epochdescription'), EEG.epochdescription = {}; res = com; end;
        if ~isempty(EEG.epoch)
            if isstruct(EEG.epoch),  l = length( EEG.epoch);
            else                     l = size( EEG.epoch, 2); 
            end;   
            if l ~= EEG.trials
                if popask( [ 'eeg_checkset error: the number of epoch indices in the epoch array/struct (' ...
                             int2str(l) ') is different from the number of epochs in the data (' int2str(EEG.trials) ').' 10 ...
                             'Should EEGLAB attempt to abort operation ?' 10 '(press Cancel to fix the problem from the commandline)']) 
                    error('eeg_checkset error: user abort');
                    %res = com;
                    %EEG.epoch = [];
                    %EEG = eeg_checkset(EEG);
                    %return;
                else
                    res = com;
                    return;
                    %error('eeg_checkset error: epoch structure size invalid');
                end;
            end;
        else
            EEG.epoch = [];
        end;

        % check ica
        % ---------
        if ~isfield(EEG, 'icachansind') 
            if isempty(EEG.icaweights)
                EEG.icachansind = []; res = com; 
            else
                EEG.icachansind = [1:EEG.nbchan]; res = com; 
            end;
        elseif isempty(EEG.icachansind) 
            if isempty(EEG.icaweights)
                EEG.icachansind = []; res = com; 
            else
                EEG.icachansind = [1:EEG.nbchan]; res = com; 
            end;
        end;
        if ~isempty(EEG.icasphere)
            if ~isempty(EEG.icaweights)
                if size(EEG.icaweights,2) ~= size(EEG.icasphere,1)
                    if popask( [ 'eeg_checkset error: number of columns in weights array (' int2str(size(EEG.icaweights,2)) ')' 10 ...
                                 'does not match the number of rows in the sphere array (' int2str(size(EEG.icasphere,1)) ')' 10 ...
                                 'Should EEGLAB remove ICA information ?' 10 '(press Cancel to fix the problem from the commandline)']) 
                        res = com;
                        EEG.icasphere = [];
                        EEG.icaweights = [];
                        EEG = eeg_checkset(EEG);
                        return;
                    else
                        error('eeg_checkset error: user abort');
                        res = com;
                        return;
                        %error('eeg_checkset error: invalid weight and sphere array sizes');
                    end;    
                end;
                if isnumeric(EEG.data)
                    if length(EEG.icachansind) ~= size(EEG.icasphere,2)
                        if popask( [ 'eeg_checkset error: number of elements in ''icachansind'' (' int2str(length(EEG.icachansind)) ')' 10 ...
                                     'does not match the number of columns in the sphere array (' int2str(size(EEG.icasphere,2)) ')' 10 ...
                                     'Should EEGLAB remove ICA information ?' 10 '(press Cancel to fix the problem from the commandline)']) 
                            res = com;
                            EEG.icasphere = [];
                            EEG.icaweights = [];
                            EEG = eeg_checkset(EEG);
                            return;
                        else
                            error('eeg_checkset error: user abort');
                            res = com;
                            return;
                            %error('eeg_checkset error: invalid weight and sphere array sizes');
                        end;    
                    end;
                    if isempty(EEG.icaact) | (size(EEG.icaact,1) ~= size(EEG.icaweights,1)) | (size(EEG.icaact,2) ~= size(EEG.data,2))
                        EEG.icaweights = double(EEG.icaweights);
                        EEG.icawinv = double(EEG.icawinv);
                        
                        % scale ICA components to RMS microvolt
                        if option_scaleicarms
                            if mean(mean(abs(pinv(EEG.icaweights * EEG.icasphere)-EEG.icawinv))) < 0.0001 
                                disp('Scaling components to RMS microvolt');
                                scaling = repmat(sqrt(mean(EEG(1).icawinv(:,:).^2))', [1 size(EEG.icaweights,2)]);
                                EEG.etc.icaweights_beforerms = EEG.icaweights;
                                EEG.etc.icasphere_beforerms = EEG.icasphere;
                                
                                EEG.icaweights = EEG.icaweights .* scaling;
                                EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);
                            end;
                        end;

                        if option_computeica
                            fprintf('eeg_checkset: recomputing the ICA activation matrix ...\n'); 
                            res = com;
                            % Make compatible with Matlab 7
                            if any(isnan(EEG.data(:)))
                                fprintf('eeg_checkset: recomputing using NaN indices [in first channel] ...\n');
                                tmpindices = find(~sum(isnan(EEG.data))); % was: tmpindices = find(~isnan(EEG.data(1,:)));
                                EEG.icaact = zeros(size(EEG.icaweights,1), size(EEG.data,2)); EEG.icaact(:) = NaN;
                                EEG.icaact(:,tmpindices) = (EEG.icaweights*EEG.icasphere)*EEG.data(:,tmpindices);
                            else
                                v = version;
                                if ~isempty(findstr(v, 'R11')) | ~isempty(findstr(v, 'R12')) | ~isempty(findstr(v, 'R13'))
                                    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                                else
                                    % recomputing as double
                                    try
                                        EEG.icaact = single((EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:));
                                    catch,
                                        disp('Weight matrix computed using single precision because of memory limitations');
                                        EEG.icaact = (single(tmpmat))*EEG.data(EEG.icachansind,:);
                                    end;
                                end;
                            end;
                            EEG.icaact    = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
                        end;
                    end;
                end;
                if isempty(EEG.icawinv)
                    EEG.icawinv = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv
                    res         = com;
                end;     
            else
                disp( [ 'eeg_checkset warning: weights matrix cannot be empty if sphere matrix is not, correcting ...' ]); 
                res = com;
                EEG.icasphere = [];
            end;
            if option_computeica
                if ~isempty(EEG.icaact) & ndims(EEG.icaact) < 3 & (EEG.trials > 1)
                    disp( [ 'eeg_checkset note: independent component made 3-D' ]); 
                    res = com;
                    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);        
                end;
            else 
                if ~isempty(EEG.icaact)
                    fprintf('eeg_checkset: removing ICA activation matrix (as per edit options) ...\n'); 
                end;
                EEG.icaact     = [];
            end;
        else
            if ~isempty( EEG.icaweights ), EEG.icaweights = []; res = com; end;
            if ~isempty( EEG.icawinv ),    EEG.icawinv = []; res = com; end;
            if ~isempty( EEG.icaact ),     EEG.icaact = []; res = com; end;
        end;
        if isempty(EEG.icaact)
            EEG.icaact = [];
        end;
        
        % check chanlocs
        % -------------
        if ~isfield(EEG, 'chaninfo')
            EEG.chaninfo = [];
        end;
        if ~isempty( EEG.chanlocs )
            
            if ~isstruct( EEG.chanlocs)
                if exist( EEG.chanlocs ) ~= 2
                    disp( [ 'eeg_checkset warning: channel file does not exist or is not in Matlab path: filename removed from EEG struct' ]); 
                    EEG.chanlocs = [];
                    res = com;
                else
                    res = com;
                    try, EEG.chanlocs = readlocs( EEG.chanlocs );
                        disp( [ 'eeg_checkset: channel file read' ]); 
                    catch, EEG.chanlocs = []; end;
                end;     
            end;
            if isstruct( EEG.chanlocs)
                if length( EEG.chanlocs) ~= EEG.nbchan & length( EEG.chanlocs) ~= EEG.nbchan+1
                    disp( [ 'eeg_checkset warning: number of channels different in data and channel file/struct: channel file/struct removed' ]); 
                    EEG.chanlocs = [];
                    res = com;
                end;
            end;
            if isstruct( EEG.chanlocs)
                if ~isstr(EEG.chanlocs(1).labels)
                    for index = 1:length(EEG.chanlocs)
                        if ~isstr(EEG.chanlocs(index).labels)
                            EEG.chanlocs(index).labels = [ 'E' int2str(EEG.chanlocs(index).labels) ];
                        end;
                    end;
                end;
                if isfield(EEG.chanlocs, 'shrink') & isempty(EEG.chanlocs(end).shrink)
                    for index = 1:length(EEG.chanlocs)
                        EEG.chanlocs(index).shrink = EEG.chanlocs(1).shrink;
                    end;
                end;
            end;
            if isfield( EEG.chanlocs, 'plotrad')
                EEG.chaninfo.plotrad = EEG.chanlocs(1).plotrad;
                EEG.chanlocs = rmfield( EEG.chanlocs, 'plotrad');
            end;
            if isfield( EEG.chanlocs, 'shrink')
                EEG.chaninfo.shrink = EEG.chanlocs(1).shrink;
                EEG.chanlocs = rmfield( EEG.chanlocs, 'shrink');
            end;

            % check if duplicate channel label
            % --------------------------------
            if isfield(EEG.chanlocs, 'labels')
                if length( { EEG.chanlocs.labels } ) > length( unique({ EEG.chanlocs.labels } ) )
                    disp('Warning: some channels have the same label');
                end;
            end;
            
            % force Nosedir to +X
            % -------------------
            if isfield(EEG.chaninfo, 'nosedir')
                if strcmpi(EEG.chaninfo.nosedir, '+x')
                    rotate = 0;
                else
                    disp('EEG checkset note for expert users: Noze direction now set to default +X in EEG.chanlocs and EEG.dipfit.');
                    if strcmpi(EEG.chaninfo.nosedir, '+y')
                        rotate = 270;
                    elseif strcmpi(EEG.chaninfo.nosedir, '-x')
                        rotate = 180;
                    else rotate = 90;
                    end;
                    for index = 1:length(EEG.chanlocs)
                        if ~isempty(EEG.chanlocs(index).theta)
                            rotategrad = rotate/180*pi;
                            coord = (EEG.chanlocs(index).Y + EEG.chanlocs(index).X*sqrt(-1))*exp(sqrt(-1)*-rotategrad);
                            EEG.chanlocs(index).Y = real(coord);
                            EEG.chanlocs(index).X = imag(coord);
                            
                            EEG.chanlocs(index).theta     = EEG.chanlocs(index).theta    -rotate;
                            EEG.chanlocs(index).sph_theta = EEG.chanlocs(index).sph_theta+rotate;
                            if EEG.chanlocs(index).theta    <-180, EEG.chanlocs(index).theta    =EEG.chanlocs(index).theta    +360; end;
                            if EEG.chanlocs(index).sph_theta>180 , EEG.chanlocs(index).sph_theta=EEG.chanlocs(index).sph_theta-360; end;
                        end;
                    end;
                    
                    if isfield(EEG, 'dipfit')
                        if isfield(EEG.dipfit, 'coord_transform')
                            if isempty(EEG.dipfit.coord_transform)
                                EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];
                            end;
                            EEG.dipfit.coord_transform(6) = EEG.dipfit.coord_transform(6)+rotategrad;
                        end;
                    end;

                end;
                EEG.chaninfo.nosedir = '+X';
            end;
            
        end;
        EEG.chaninfo.icachansind = EEG.icachansind; % just a copy for programming convinience

        %if ~isfield(EEG, 'urchanlocs')
        %    EEG.urchanlocs = EEG.chanlocs;
        %    for index = 1:length(EEG.chanlocs)
        %        EEG.chanlocs(index).urchan = index;
        %    end;
        %    disp('eeg_checkset note: creating backup chanlocs structure (urchanlocs)');
        %end;

        % check reference
        % ---------------
        if ~isfield(EEG, 'ref')
            EEG.ref = 'common';
        end;
        if isstr(EEG.ref) & strcmpi(EEG.ref, 'common')
            if length(EEG.chanlocs) > EEG.nbchan
                disp('Extra common reference electrode location detected');
                EEG.ref = EEG.nbchan+1;
            end;
        end;

        % DIPFIT structure
        % ----------------
        if ~isfield(EEG, 'dipfit') 
            EEG.dipfit = []; res = com; 
        elseif exist('dipfitdefs')
            if isfield(EEG.dipfit, 'vol') & ~isfield(EEG.dipfit, 'hdmfile')
                if exist('pop_dipfit_settings')
                    disp('Old DIPFIT structure detected: converting to DIPFIT 2 format');
                    EEG.dipfit.hdmfile     = template_models(1).hdmfile;
                    EEG.dipfit.coordformat = template_models(1).coordformat;
                    EEG.dipfit.mrifile     = template_models(1).mrifile;
                    EEG.dipfit.chanfile    = template_models(1).chanfile;
                    EEG.dipfit.coord_transform = [];
                    res = com;
                end;
            end;
            if isfield(EEG.dipfit, 'hdmfile')
                if length(EEG.dipfit.hdmfile) > 8
                    if strcmpi(EEG.dipfit.hdmfile(end-8), template_models(1).hdmfile(end-8)), EEG.dipfit.hdmfile = template_models(1).hdmfile; end;
                    if strcmpi(EEG.dipfit.hdmfile(end-8), template_models(2).hdmfile(end-8)), EEG.dipfit.hdmfile = template_models(2).hdmfile; end;
                end;
                if length(EEG.dipfit.mrifile) > 8
                    if strcmpi(EEG.dipfit.mrifile(end-8), template_models(1).mrifile(end-8)), EEG.dipfit.mrifile = template_models(1).mrifile; end;
                    if strcmpi(EEG.dipfit.mrifile(end-8), template_models(2).mrifile(end-8)), EEG.dipfit.mrifile = template_models(2).mrifile; end;
                end;
                if length(EEG.dipfit.chanfile) > 8
                    if strcmpi(EEG.dipfit.chanfile(end-8), template_models(1).chanfile(end-8)), EEG.dipfit.chanfile = template_models(1).chanfile; end;
                    if strcmpi(EEG.dipfit.chanfile(end-8), template_models(2).chanfile(end-8)), EEG.dipfit.chanfile = template_models(2).chanfile; end;
                end;
            end;

            if isfield(EEG.dipfit, 'coord_transform')
                if isempty(EEG.dipfit.coord_transform)
                    EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];
                end;
            elseif ~isempty(EEG.dipfit)
                EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];
            end;

        end;

        % EEG.times (only for epoched datasets)
        % ---------
        if (EEG.trials > 1)
            EEG.times = linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);
        else
            EEG.times = [];
        end;

        if ~isfield(EEG, 'history')    EEG.history    = ''; res = com; end;
        if ~isfield(EEG, 'splinefile') EEG.splinefile = ''; res = com; end;
        if ~isfield(EEG, 'icasplinefile') EEG.icasplinefile = ''; res = com; end;
        if ~isfield(EEG, 'saved')      EEG.saved      = 'no'; res = com; end;
        if ~isfield(EEG, 'subject')    EEG.subject    = ''; res = com; end;
        if ~isfield(EEG, 'condition')  EEG.condition  = ''; res = com; end;
        if ~isfield(EEG, 'group')      EEG.group      = ''; res = com; end;
        if ~isfield(EEG, 'session')    EEG.session    = []; res = com; end;
        if ~isfield(EEG, 'urchanlocs') EEG.urchanlocs = []; res = com; end;
        if ~isfield(EEG, 'specdata')   EEG.specdata   = []; res = com; end;
        if ~isfield(EEG, 'specicaact') EEG.specicaact = []; res = com; end;
        if ~isfield(EEG, 'comments')   EEG.comments   = ''; res = com; end;
        if ~isfield(EEG, 'etc'     )   EEG.etc        = []; res = com; end;
        if ~isfield(EEG, 'urevent' )   EEG.urevent    = []; res = com; end;
        if ~isfield(EEG, 'ref') | isempty(EEG.ref) EEG.ref = 'common'; res = com; end;

        % create fields if absent
        % -----------------------
        if ~isfield(EEG, 'reject')                    EEG.reject.rejjp = []; res = com; end;

        listf = { 'rejjp' 'rejkurt' 'rejmanual' 'rejthresh' 'rejconst', 'rejfreq' ...
                  'icarejjp' 'icarejkurt' 'icarejmanual' 'icarejthresh' 'icarejconst', 'icarejfreq'};
        for index = 1:length(listf)    
            if ~isfield(EEG.reject, listf{index}),    EEG.reject = setfield(EEG.reject, listf{index}, []); res = com; end;
            elecfield = [listf{index} 'E'];
            if ~isfield(EEG.reject, elecfield),     EEG.reject = setfield(EEG.reject, elecfield, []); res = com; end;
            % check if electrode array is empty with rejection array is not
            if ~isempty(getfield(EEG.reject, listf{index})) & isempty(getfield(EEG.reject, elecfield))
                nbchan = fastif( strcmp(listf{index}, 'ica'), size(EEG.icaweights,1), EEG.nbchan);
                EEG.reject = setfield(EEG.reject, elecfield, zeros(nbchan, length(getfield(EEG.reject, listf{index})))); res = com;
            end;
        end;
        if ~isfield(EEG.reject, 'rejglobal')        EEG.reject.rejglobal = []; res = com; end;
        if ~isfield(EEG.reject, 'rejglobalE')        EEG.reject.rejglobalE = []; res = com; end;

        % default colors for rejection
        % ----------------------------
        if ~isfield(EEG.reject, 'rejmanualcol')   EEG.reject.rejmanualcol = [1.0000    1     0.783]; res = com; end;
        if ~isfield(EEG.reject, 'rejthreshcol')   EEG.reject.rejthreshcol = [0.8487    1.0000    0.5008]; res = com; end;
        if ~isfield(EEG.reject, 'rejconstcol')    EEG.reject.rejconstcol  = [0.6940    1.0000    0.7008]; res = com; end;
        if ~isfield(EEG.reject, 'rejjpcol')       EEG.reject.rejjpcol     = [1.0000    0.6991    0.7537]; res = com; end;
        if ~isfield(EEG.reject, 'rejkurtcol')     EEG.reject.rejkurtcol   = [0.6880    0.7042    1.0000]; res = com; end;
        if ~isfield(EEG.reject, 'rejfreqcol')     EEG.reject.rejfreqcol   = [0.9596    0.7193    1.0000]; res = com; end;
        if ~isfield(EEG.reject, 'disprej')        EEG.reject.disprej      = { }; end;
        
        if ~isfield(EEG, 'stats')           EEG.stats.jp = []; res = com; end;
        if ~isfield(EEG.stats, 'jp')        EEG.stats.jp = []; res = com; end;
        if ~isfield(EEG.stats, 'jpE')       EEG.stats.jpE = []; res = com; end;
        if ~isfield(EEG.stats, 'icajp')     EEG.stats.icajp = []; res = com; end;
        if ~isfield(EEG.stats, 'icajpE')    EEG.stats.icajpE = []; res = com; end;
        if ~isfield(EEG.stats, 'kurt')      EEG.stats.kurt = []; res = com; end;
        if ~isfield(EEG.stats, 'kurtE')     EEG.stats.kurtE = []; res = com; end;
        if ~isfield(EEG.stats, 'icakurt')   EEG.stats.icakurt = []; res = com; end;
        if ~isfield(EEG.stats, 'icakurtE')  EEG.stats.icakurtE = []; res = com; end;

        % component rejection
        % -------------------
        if ~isfield(EEG.stats, 'compenta')        EEG.stats.compenta = []; res = com; end;
        if ~isfield(EEG.stats, 'compentr')        EEG.stats.compentr = []; res = com; end;
        if ~isfield(EEG.stats, 'compkurta')       EEG.stats.compkurta = []; res = com; end;
        if ~isfield(EEG.stats, 'compkurtr')       EEG.stats.compkurtr = []; res = com; end;
        if ~isfield(EEG.stats, 'compkurtdist')    EEG.stats.compkurtdist = []; res = com; end;
        if ~isfield(EEG.reject, 'threshold')      EEG.reject.threshold = [0.8 0.8 0.8]; res = com; end;
        if ~isfield(EEG.reject, 'threshentropy')  EEG.reject.threshentropy = 600; res = com; end;
        if ~isfield(EEG.reject, 'threshkurtact')  EEG.reject.threshkurtact = 600; res = com; end;
        if ~isfield(EEG.reject, 'threshkurtdist') EEG.reject.threshkurtdist = 600; res = com; end;
        if ~isfield(EEG.reject, 'gcompreject')    EEG.reject.gcompreject = []; res = com; end;
        if length(EEG.reject.gcompreject) ~= size(EEG.icaweights,1)
            EEG.reject.gcompreject = zeros(1, size(EEG.icaweights,1));
        end;

        % remove old fields
        % -----------------
        if isfield(EEG, 'averef'), EEG = rmfield(EEG, 'averef'); end;
        if isfield(EEG, 'rt'    ), EEG = rmfield(EEG, 'rt');     end;

        % store in new structure
        % ----------------------
        ALLEEGNEW(inddataset) = EEG;
    end;
end;

% recorder fields
% ---------------
fieldorder = { 'setname' ...
               'filename' ...
               'filepath' ...
               'subject' ...
               'group' ...
               'condition' ...
               'session' ...
               'comments' ...
               'nbchan' ...
               'trials' ...
               'pnts' ...
               'srate' ...
               'xmin' ...
               'xmax' ...
               'times' ...
               'data' ...
               'icaact' ...
               'icawinv' ...
               'icasphere' ...
               'icaweights' ...
               'icachansind' ...
               'chanlocs' ...
               'urchanlocs' ...
               'chaninfo' ...
               'ref' ...
               'event' ...
               'urevent' ...
               'eventdescription' ...
               'epoch' ...
               'epochdescription' ...
               'reject' ...
               'stats' ...
               'specdata' ...
               'specicaact' ...
               'splinefile' ...
               'icasplinefile' ...
               'dipfit' ...
               'history' ...
               'saved' ... 
               'etc' };

fieldorder = [ fieldorder, setdiff(fieldnames(EEG), fieldorder) ];
try 
    ALLEEGNEW = orderfields(ALLEEGNEW, fieldorder);
    EEG = ALLEEGNEW;
catch, disp('error trying to order fields');
end;

if exist('ALLEEGNEW')
    EEG = ALLEEGNEW;
end;

return;

function num = popask( text )
     ButtonName=questdlg2( text, ...
            'Confirmation', 'Cancel', 'Yes','Yes');
     switch lower(ButtonName),
          case 'cancel', num = 0;
          case 'yes',    num = 1;
     end;

function res = mycellfun(com, vals, classtype);
    res = zeros(1, length(vals));
    switch com
     case 'isempty', 
      for index = 1:length(vals), res(index) = isempty(vals{index}); end;
     case 'isclass'
      if strcmp(classtype, 'double')
          for index = 1:length(vals), res(index) = isnumeric(vals{index}); end;
      else 
          error('unknown cellfun command');
      end;
     otherwise error('unknown cellfun command');
    end;
        
    
