% eeg_checkset() - check the consistency of fields of an EEG dataset 
%
% Structure of an EEG dataset under EEGLAB:
%    EEG.data         - two-dimensional continuous data array (chans, frames)
%                       OR three-dim. epoched data array (chans, frames, epochs)
%    EEG.setname      - name of the dataset
%    EEG.filename     - filename of the dataset
%    EEG.filepath     - filepath of the dataset
%    EEG.chanlocs     - structure array containing names and positions 
%                       of the channels on the scalp
%    EEG.pnts         - number of time points (data frames) per epoch (trial).
%                       OR if data is continuous, total number of time points
%    EEG.nbchan       - number of channels
%    EEG.trials       - number of epochs (trials) in the dataset. If data
%                       is continuous, automatically set to 1.
%    EEG.srate        - channel sampling rate (in Hz)
%    EEG.xmin         - epoch start time (in seconds)
%    EEG.xmax         - epoch end time (in seconds)
%    EEG.times        - time vector (one value per time point)
%    EEG.ref          - ['common'|'averef'|integer] reference index or type
%    EEG.comments     - comments about the dataset
%
% ICA variables:
%    EEG.icaact       - ICA activations (components, frames, epochs)
%                       [] means compute_ica option is set to 0 under
%                       EEGLAB options -> activations are computed on the fly.
%    EEG.icasphere    - sphere array returned by linear (ICA) decomposition
%    EEG.icaweights   - weight array returned by linear (ICA) decomposition
%    EEG.icawinv      - inverse (ICA) weight matrix giving the projected
%                       activity of the components to the electrodes.
%                       NOTE: Any linear unmixing matrix may be used. 
%
% Event and epoch structures:    
%       EEG.event     - event structure (any number of events per epoch)
%       EEG.epoch     - epoch structure (one structure per epoch)
%       EEG.eventdescription - cell array of strings describing event fields.
%       EEG.epochdescription - cell array of strings describing epoch fields.
%       --> See the http://sccn.ucsd.edu/eeglab/maintut/eeglabscript.html 
%           for details
%
% Variables used for manual and semi-automatic data rejection:
%      EEG.specdata          - data spectrum for every single trial
%      EEG.specica           - data spectrum for every single trial
%      EEG.stats.kurtc       - component kurtosis values
%      EEG.stats.kurtg       - global kurtosis of components      
%      EEG.stats.kurta       - kurtosis of accepted epochs      
%      EEG.stats.kurtr       - kurtosis of rejected epochs      
%      EEG.stats.kurtd       - kurtosis of spatial distribution      
%      EEG.reject.entropy    - entropy of epochs  
%      EEG.reject.entropyc   - entropy of components
%      EEG.reject.threshold  - rejection thresholds 
%      EEG.reject.icareject  - epochs rejected by ICA criteria
%      EEG.reject.gcompreject - rejected ICA components
%      EEG.reject.sigreject  - epochs rejected by single-channel criteria
%      EEG.reject.elecreject - epochs rejected by raw data criteria
%
%      EEG.reject.compreject - deprecated
%      EEG.reject.comptrial  - deprecated
%      EEG.reject.eegentropy - deprecated
%      EEG.reject.eegkurt    - deprecated
%      EEG.reject.eegkurtg   - deprecated
%
% Usage:
%       >> [EEGOUT, res] = eeg_checkset( EEG ); % check consistency of EEG
%
% Inputs:
%       EEG        - EEGLAB dataset structure
%
% Outputs:
%       EEGOUT     - output EEGLAB dataset
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

if nargin < 1
    help eeg_checkset;
    return;
end;

% checking multiple datasets
if isempty(EEG), return; end;
if ~isfield(EEG, 'data'), return; end;
if length(EEG) > 1
    for index = 1:length(EEG)
        if ~isempty(EEG(index))
            if ~isempty( varargin)
                [TMP, res] = eeg_checkset(EEG(index), varargin{:});
            else
                [TMP, res] = eeg_checkset(EEG(index));
            end;
            [EEG TMP] = eeg_store(EEG, TMP, index);
        end;
    end;
    return;
end;
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

com = sprintf('%s = eeg_checkset( %s );', inputname(1), inputname(1));
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

% read data if necessary
% ----------------------
if isa(EEG.data, 'single')
    EEG.data = double(EEG.data);
end;
if isstr(EEG.data)
    fid = fopen([EEG.filepath EEG.data], 'r', 'ieee-le'); %little endian (see also pop_saveset)
    if fid == -1
        disp(['file ' [EEG.filepath EEG.data] ' not found, trying local folder']);
        fid = fopen(EEG.data, 'r', 'ieee-le'); %little endian (see also pop_saveset)
        if fid == -1
            errordlg2(['Cannot open data file ''' [EEG.data] ''''], 'error');
            error('File not found');
        end;
        fprintf('Reading float file ''%s''...\n', [EEG.data]);
    else 
        fprintf('Reading float file ''%s''...\n', [EEG.filepath EEG.data]);
    end;
    EEG.data = fread(fid, [EEG.nbchan Inf], 'float32');
end;

% verify the type of the variables
% --------------------------------
    % data dimensions -------------------------
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

    % parameters coherence -------------------------
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
    if ~isempty(EEG.event) & ~isfield(EEG, 'urevent')
        EEG.urevent = EEG.event;
        disp('eeg_checkset note: creating backup event table (urevent)');
        for index = 1:length(EEG.event)
            EEG.event(index).urevent = index;
        end;
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
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
    if ~isempty(EEG.icasphere)
        if ~isempty(EEG.icaweights)
            if size(EEG.icaweights,2) ~= size(EEG.icasphere,1)
                  if popask( [ 'eeg_checkset error: number of columns in weights array (' int2str(size(EEG.icaweights,2)) 10 ')' ...
                   'does not match the number of rows in the sphere array (' int2str(size(EEG.icasphere,1)) ')' 10 ...
                   'Should EEGLAB attempt to abort operation ?' 10 '(press Cancel to fix the problem from the commandline)']) 
                    error('eeg_checkset error: user abort');
                    %res = com;
                    %EEG.icasphere = [];
                    %EEG.icaweights = [];
                    %EEG = eeg_checkset(EEG);
                    %return;
                else
                    res = com;
                    return;
                    %error('eeg_checkset error: invalid weight and sphere array sizes');
                end;    
            end;
            if size(EEG.icasphere,2) ~= size(EEG.data,1)
                   disp( [ 'eeg_checkset warning: number of columns in ica matrix (' int2str(size(EEG.icasphere,2)) ...
                   ') does not match the number of rows in data (' int2str(size(EEG.data,1)) ')' ]); 
                res = com;
            end;
            if isempty(EEG.icaact) | (size(EEG.icaact,1) ~= size(EEG.icaweights,1)) | (size(EEG.icaact,2) ~= size(EEG.data,2))
                if size(EEG.data,1) ~= size(EEG.icasphere,2)
                       if popask( [ 'eeg_checkset error: number of columns in sphere array (' int2str(size(EEG.icasphere,2)) 10 ')' ...
                       'does not match the number of rows in data(' int2str(size(EEG.data,1)) ')' 10 ...
                       'Do you want to want to abort operation ?' 10 '(press Cancel to fix the problem from the commandline)']) 
                        error('eeg_checkset error: user abort');
                        %res = com;
                        %EEG.icasphere = [];
                        %EEG.icaweights = [];
                        %EEG = eeg_checkset(EEG);
                        %return;
                    else
                        res = com;
                        return;
                        %error('eeg_checkset error: invalid weight and sphere array size');
                    end;    
                end;
                if option_computeica
                     fprintf('eeg_checkset: recomputing the ICA activation matrix ...\n'); 
                    res = com;
                    if any(isnan(EEG.data(:)))
                        fprintf('eeg_checkset: recomputing using NaN indices in first channel ...\n'); 
                        tmpindices = find(~isnan(EEG.data(1,:)));
                        EEG.icaact = zeros(size(EEG.icaweights,1), size(EEG.data,2)); EEG.icaact(:) = NaN;
                        EEG.icaact(:,tmpindices) = (EEG.icaweights*EEG.icasphere)*EEG.data(:,tmpindices);
                    else
                        EEG.icaact    = (EEG.icaweights*EEG.icasphere)*EEG.data(:,:);
                    end;
                    EEG.icaact    = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
                end;
             end;
            if isempty(EEG.icawinv)
                EEG.icawinv    = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv
                res = com;
            end;     
        else
               disp( [ 'eeg_checkset warning: weights matrix cannot be empty if sphere matrix is not, correcting ...' ]); 
            res = com;
               EEG.icasphere = [];
        end;
        if option_computeica
            if (ndims(EEG.icaact)) < 3 & (EEG.trials > 1)
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
                    EEG.chanlocs(index).labels = int2str(EEG.chanlocs(index).labels);
                end;
            end;
        end;
        if isfield(EEG.chanlocs, 'shrink') & isempty(EEG.chanlocs(end).shrink)
            for index = 1:length(EEG.chanlocs)
                EEG.chanlocs(index).shrink = EEG.chanlocs(1).shrink;
            end;
        end;
    end;
end;

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

% EEG.times (only for epoched datasets)
% ---------
if (EEG.trials > 1)
    EEG.times = linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);
else
    if isfield(EEG, 'times')
        EEG = rmfield(EEG, 'times');
    end;
end;
if ~isfield(EEG, 'specdata') EEG.specdata = []; res = com; end;
if ~isfield(EEG, 'specicaact') EEG.specicaact = []; res = com; end;
if ~isfield(EEG, 'comments') EEG.comments = ''; res = com; end;
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
if ~isfield(EEG.reject, 'rejmanualcol')        EEG.reject.rejmanualcol = [1.0000    1     0.783]; res = com; end;
if ~isfield(EEG.reject, 'rejthreshcol')        EEG.reject.rejthreshcol = [0.8487    1.0000    0.5008]; res = com; end;
if ~isfield(EEG.reject, 'rejconstcol')        EEG.reject.rejconstcol  = [0.6940    1.0000    0.7008]; res = com; end;
if ~isfield(EEG.reject, 'rejjpcol')            EEG.reject.rejjpcol     = [1.0000    0.6991    0.7537]; res = com; end;
if ~isfield(EEG.reject, 'rejkurtcol')        EEG.reject.rejkurtcol   = [0.6880    0.7042    1.0000]; res = com; end;
if ~isfield(EEG.reject, 'rejfreqcol')        EEG.reject.rejfreqcol   = [0.9596    0.7193    1.0000]; res = com; end;
if ~isfield(EEG.reject, 'disprej')            EEG.reject.disprej      = { }; end;
    
if ~isfield(EEG, 'stats')            EEG.stats.jp = []; res = com; end;
if ~isfield(EEG.stats, 'jp')        EEG.stats.jp = []; res = com; end;
if ~isfield(EEG.stats, 'jpE')        EEG.stats.jpE = []; res = com; end;
if ~isfield(EEG.stats, 'icajp')        EEG.stats.icajp = []; res = com; end;
if ~isfield(EEG.stats, 'icajpE')    EEG.stats.icajpE = []; res = com; end;
if ~isfield(EEG.stats, 'kurt')        EEG.stats.kurt = []; res = com; end;
if ~isfield(EEG.stats, 'kurtE')        EEG.stats.kurtE = []; res = com; end;
if ~isfield(EEG.stats, 'icakurt')    EEG.stats.icakurt = []; res = com; end;
if ~isfield(EEG.stats, 'icakurtE')    EEG.stats.icakurtE = []; res = com; end;

% component rejection
% -------------------
if ~isfield(EEG.stats, 'compenta')        EEG.stats.compenta = []; res = com; end;
if ~isfield(EEG.stats, 'compentr')        EEG.stats.compentr = []; res = com; end;
if ~isfield(EEG.stats, 'compkurta')        EEG.stats.compkurta = []; res = com; end;
if ~isfield(EEG.stats, 'compkurtr')        EEG.stats.compkurtr = []; res = com; end;
if ~isfield(EEG.stats, 'compkurtdist')    EEG.stats.compkurtdist = []; res = com; end;
if ~isfield(EEG.reject, 'threshold')        EEG.reject.threshold = [0.8 0.8 0.8]; res = com; end;
if ~isfield(EEG.reject, 'threshentropy')    EEG.reject.threshentropy = 600; res = com; end;
if ~isfield(EEG.reject, 'threshkurtact')    EEG.reject.threshkurtact = 600; res = com; end;
if ~isfield(EEG.reject, 'threshkurtdist')    EEG.reject.threshkurtdist = 600; res = com; end;
if ~isfield(EEG.reject, 'gcompreject')        EEG.reject.gcompreject = []; res = com; end;
if length(EEG.reject.gcompreject) ~= size(EEG.icaweights,1)
    EEG.reject.gcompreject = zeros(1, size(EEG.icaweights,1));
end;

% component rejection
% -------------------
% additional checks
% -----------------
if ~isempty( varargin)
    for index = 1:length( varargin )
        switch varargin{ index }
         case 'data',; % already done at the top 
         case 'contdata',;
          if EEG.trials > 1
              errordlg2(strvcat('Error: function only works on continuous data'), 'Error');
              error('eeg_checkset error: data is not continuous'); return;
          end;
         case 'ica', 
          if isempty(EEG.icaweights)
              if ~popask(strvcat('No ICA weights. Compute now?', '(then go back to the function you just called)'))
                  errordlg2('eeg_checkset: ICA components must be derived before running that function'); 
                  error('no ICA components'); return; 
              end;
              [EEG res] = pop_runica(EEG);
              res = [ inputname(1) ' = eeg_checkset('  inputname(1) '); ' res ];
          end;
         case 'epoch', 
          if EEG.trials == 1
              errordlg2(strvcat('Epochs must be extracted before running that function', 'Use Tools > Extract epochs'), 'Error');
              error('eeg_checkset error: epochs must be extracted before running that function'); return
          end;
         case 'besa', 
          if ~isfield(EEG, 'sources')
              errordlg2(strvcat('No dipole information', '1) Component maps must be exported: Tools > Localize ... BESA > Export ...' ...
                                , '2) BESA must be run to localize the equivalent dipoles', ...
                                '3) BESA dipoles must be imported: Tools > Localize ... BESA > Import ...'), 'Error');
              error('eeg_checkset error: no BESA dipole information'); return
          end;
         case 'event', 
          if isempty(EEG.event)
              errordlg2(strvcat('Cannot process if no events. First add events.', 'Use File > Import event info or > Import epoch info'), 'Error');
              error('eeg_checkset: no events'); return;
          end;
         case 'chanloc', 
          if isempty(EEG.chanlocs)
              errordlg2( strvcat('Cannot process without channel location file.', ...
                         'Enter the name of the file via "Edit > Edit dataset info".', ...
                         'For file format, enter ''>> help readlocs'' from the command line.'), 'Error');
              error('eeg_checkset: cannot process dataset without channel location file.'); return;
          end;
         case 'chanlocs_homogenous', 
          if isempty(EEG.chanlocs)
              errordlg2( strvcat('Cannot process without channel location file.', ...
                         'Enter the name of the file via "Edit > Edit dataset info".', ...
                         'For file format, enter ''>> help readlocs'' from the command line.'), 'Error');
              error('eeg_checkset: cannot process dataset without channel location file.'); return;
          end;
          if ~isfield(EEG.chanlocs, 'X') | isempty(EEG.chanlocs(1).X)
              EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
              res = [ inputname(1) ' = eeg_checkset('  inputname(1) ', ''chanlocs_homogenous'' ); ' ];
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
              disp('eeg_checkset note: creating backup event table (urevent)');
              EEG.urevent = EEG.event;
              for index = 1:length(EEG.event)
                  EEG.event(index).urevent = index;
              end;
          end;
         case 'eventconsistency',          
          if isempty(EEG.event), return; end;
          
          % remove the events which latency are out of boundary
          % ---------------------------------------------------
          if isfield(EEG.event, 'latency')
              try, alllatencies = cell2mat( { EEG.event.latency } );
              catch, error('Checkset: error empty latency entry for new events added by user');
              end;
              I1 = find(alllatencies < 1);
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
              try, allepochs = cell2mat( { EEG.event.epoch } );
                  removeevent = find( allepochs < 1 | allepochs > EEG.trials);
                  if ~isempty(removeevent)
                      disp([ 'eeg_checkset warning: ' int2str(length(removeevent)) ' event had invalid epoch number and were removed']);
                  end;
              catch, 
                  for indexevent = 1:length(EEG.event)
                      if isempty( EEG.event(indexevent).epoch ) | ~isnumeric(EEG.event(indexevent).epoch) ...
                              | EEG.event(indexevent).epoch < 1 | EEG.event(indexevent).epoch > EEG.trials
                          removeevent = [removeevent indexevent];
                          disp([ 'eeg_checkset warning: event ' int2str(indexevent) ' has invalid epoch number: removed']);
                      end;
                  end;
              end;
              EEG.event(removeevent) = [];
              allepochs = cell2mat( { EEG.event.epoch } );
              
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
                      fprintf(['eeg_checkset: filling with values of other events in the same epochs\n']);
                  end;
               end;
          end;
          
          % uniformize fields (str or int) if necessary
          % -------------------------------------------
          allfields = fieldnames(EEG.event);
          for index = 1:length(allfields)
              eval(['allvalues = { EEG.event.' allfields{index} ' };']);
              try,   eval(['valreal = cellfun(''isclass'', allvalues, ''double'');']);
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
                  fprintf('eeg_checkset: value format of event field ''%s'' made uniform\n', allfields{index});
                  % get the field content
                  % ---------------------
                  for indexevent = 1:length(EEG.event)
                      if valreal(indexevent)
                          EEG.event = setfield(EEG.event, { indexevent }, allfields{index}, num2str(allvalues{indexevent}) );
                      end;
                  end;
              end;
          end;
          
          % build epoch structure
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
                          valfield = eeg_point2lat(cell2mat(valfield),trial,EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
                          valfield = mat2cell(valfield);
                      end;
                      if ~isempty(valfield)
                          if maxlen == 1, EEG.epoch = setfield(EEG.epoch, { trial }, ['event' eventfields{fieldnum}], valfield{1});
                          else            EEG.epoch = setfield(EEG.epoch, { trial }, ['event' eventfields{fieldnum}], valfield);
                          end;
                      end;
                  end;
              end;    
          end;
          catch, errordlg2(['warning: minor problem encountered when generating' 10 ...
                        'epoch information (only useful for users using command line scripts)']); return;
          end;
         otherwise, error('eeg_checkset: unknown option');
        end;        
    end;
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
    
