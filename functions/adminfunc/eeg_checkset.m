% eeg_checkset() - check consistency of EEGLAB dataset fields.
%
% Structure of an EEG dataset under EEGLAB:
%	EEG.data     	- three-dimensional data array (chans, frames, epochs)
%	EEG.setname  	- name of the dataset
%	EEG.filename 	- filename of the dataset
%	EEG.filepath    - filepath of the dataset
%	EEG.namechan 	- channel labels (will be deprecated)
%	EEG.chanlocs  	- name of file containing names and positions 
%                         of the channels on the scalp
%	EEG.pnts     	- number of frames (time points) per epoch (trial)
%	EEG.nbchan     	- number of channels in each epoch
%	EEG.trials     	- number of epochs (trials) in the dataset
%	EEG.srate      	- sampling rate (in Hz)
%	EEG.xmin      	- epoch start time (in seconds)
%	EEG.xmax      	- epoch end time (in seconds)
%
% ICA variables:
%	EEG.icaact      - ICA activations (components x frames x epochs)  
%	EEG.icasphere   - sphere array returned by ICA
%	EEG.icaweights  - weight array returned by ICA
%	EEG.icawinv     - inverse ICA weight matrix giving the projected
%                         activity of the components at the electrodes.
%                         NOTE: Any linear unmixing matrix may be used.
% Event and epoch structures:	
%       EEG.event       - event structure (any number of events per epoch)
%       EEG.epoch       - epoch structure (one structure per epoch)
%   --> See the web page http://sccn.ucsd.edu/eeglab/xxx.html for details
%
% Variables used for manual and semi-automatic data rejection:
%	EEG.stats.kurtc         - component kurtosis values
%	EEG.stats.kurtg         - global kurtosis of components      
%	EEG.stats.kurta         - kurtosis of accepted epochs      
%	EEG.stats.kurtr         - kurtosis of rejected epochs      
%	EEG.stats.kurtd         - kurtosis of spatial distribution      
%	EEG.reject.entropy  	- entropy of epochs  
%	EEG.reject.entropyc   	- entropy of components
%	EEG.reject.threshold    - rejection thresholds 
%	EEG.reject.icareject    - epochs rejected by ICA criteria
%	EEG.reject.gcompreject  - rejected ICA components
%	EEG.reject.sigreject    - epochs rejected by single-channel criteria
%	EEG.reject.elecreject   - epochs rejected by raw data criteria
%	EEG.reject.compreject   - deprecated
%	EEG.reject.comptrial    - deprecated
%	EEG.reject.eegentropy   - deprecated
%	EEG.reject.eegkurt      - deprecated
%	EEG.reject.eegkurtg     - deprecated
%
% Usage:
%       >> [EEGOUT, res] = eeg_checkset( EEG );
%
% Inputs:
%       EEG        - dataset structure
%
% Outputs:
%       EEGOUT     - output dataset
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
% EEG.averef = 'No' by default -sm
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
% add furhter check for EEG.averef
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

if ~isempty( varargin)
    if isempty(EEG.data)
        helpdlg('Empty dataset -> File / Import data or File / Load existing dataset', 'Error');
        error('eeg_checkset: empty dataset');
    end;    
end;

com = sprintf('%s = eeg_checkset( %s );', inputname(1), inputname(1));
res = [];

% read data if necessary
% ----------------------
if isstr(EEG.data)
	fid = fopen(EEG.data, 'r', 'ieee-le');
	if fid == -1
		error(['Can not open data file ''' EEG.data ''', check directory']);
	end;
	fprintf('Reading float file ''%s''...', EEG.data);
	EEG.data = fread(fid, [EEG.nchan Inf], 'float32');
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
                          'Should EEGLAB attempt to adjust?' 10 '(press Cancel to fix the problem from the command line)']) 
                res = com;
                EEG.pnts = size(EEG.data,2);
                EEG = eeg_checkset(EEG);
                return;
            else
              	 error( 'eeg_checkset error: number of points does not divide the number of columns in data');
            end;  	  
      else
        if EEG.trials > 1
       		disp( 'eeg_checkset note: number of dimensions in data increased to 3'); 
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
    if 	round(EEG.srate*(EEG.xmax-EEG.xmin)+1) ~= EEG.pnts	  	
       fprintf( 'eeg_checkset warning: inconsistency (xmax-xmin)*rate+1 (=%f) must be equal to the number of frames (=%d); xmax corrected\n', ...
          EEG.srate*(EEG.xmax-EEG.xmin)+1, EEG.pnts); 
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
                          'Should EEGLAB remove all events ?' 10 '(press Cancel to fix the problem from the command line)']) 
                res = com;
                EEG.event = [];
                EEG = eeg_checkset(EEG);
                return;
            else
                error('eeg_checkset error: no epoch field in event structure');
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
                   'Do you want to remove them ?' 10 '(press Cancel to fix the problem from the command line)']) 
                res = com;
                EEG.epoch = [];
                EEG = eeg_checkset(EEG);
                return;
            else
                error('eeg_checkset error: epoch structure size invalid');
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
 	  			if popask( [ 'eeg_checkset error: number of columns in weights array (' int2str(size(EEG.icaweights,2)) 10
 	  			') does not match the number of rows in the sphere array (' int2str(size(EEG.icasphere,1)) ')' 10 ...
 	  			'Should EEGLAB clear these matrices?' 10 '(press Cancel to fix the problem from the command line)']) 
                    res = com;
                    EEG.icasphere = [];
                    EEG.icaweights = [];
                    EEG = eeg_checkset(EEG);
                    return;
                else
                    error('eeg_checkset error: invalid weight and sphere array sizes');
                end;    
			end;
			if size(EEG.icasphere,2) ~= size(EEG.data,1)
 	  			disp( [ 'eeg_checkset warning: number of columns in ica matrix (' int2str(size(EEG.icasphere,2)) ...
 	  			') does not match the number of rows in data (' int2str(size(EEG.data,1)) ')' ]); 
                res = com;
			end;
			if isempty(EEG.icaact) | (size(EEG.icaact,1) ~= size(EEG.icaweights,1)) | (size(EEG.icaact,2) ~= size(EEG.data,2))
                if size(EEG.data,1) ~= size(EEG.icasphere,2)
	 	  			if popask( [ 'eeg_checkset error: number of columns in sphere array (' int2str(size(EEG.icasphere,2)) 10
	 	  			') does not match the number of rows in data(' int2str(size(EEG.data,1)) ')' 10 ...
	 	  			'Do you want to remove them ?' 10 '(press Cancel to fix the problem from the command line)']) 
	                    res = com;
	                    EEG.icasphere = [];
	                    EEG.icaweights = [];
	                    EEG = eeg_checkset(EEG);
	                    return;
	                else
	                    error('eeg_checkset error: invalid weight and sphere array size');
	                end;    
                end;
                if option_computeica
 	    			fprintf('eeg_checkset: recomputing ica activation matrix ...\n'); 
                    res = com;
                    EEG.icaact     = (EEG.icaweights*EEG.icasphere)*EEG.data(:,:);
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
				disp( [ 'eeg_checkset note: number of dimensions in independent component array increased to 3' ]); 
				res = com;
				EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);		
			end;
		else 
			if ~isempty(EEG.icaact)
				fprintf('eeg_checkset: removing ica activation matrix (see edit options) ...\n'); 
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
        if length( EEG.chanlocs) ~= EEG.nbchan
			disp( [ 'eeg_checkset warning: number of channels different in data and channel file/struct: channel file/struct removed' ]); 
	        EEG.chanlocs = [];
	        res = com;
	    end;
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
if ~isfield(EEG, 'averef') | isempty(EEG.averef) EEG.averef = 'No'; res = com; end;

% create fields if absent
% -----------------------

if ~isfield(EEG, 'reject')					EEG.reject.rejjp = []; res = com; end;
if ~isfield(EEG.reject, 'rejjp')			EEG.reject.rejjp = []; res = com; end;
if ~isfield(EEG.reject, 'rejjpE')			EEG.reject.rejjpE = []; res = com; end;
if ~isfield(EEG.reject, 'rejkurt')			EEG.reject.rejkurt = []; res = com; end;
if ~isfield(EEG.reject, 'rejkurtE')			EEG.reject.rejkurtE = []; res = com; end;
if ~isfield(EEG.reject, 'rejmanual')		EEG.reject.rejmanual = []; res = com; end;
if ~isfield(EEG.reject, 'rejmanualE')		EEG.reject.rejmanualE = []; res = com; end;
if ~isfield(EEG.reject, 'rejthresh')		EEG.reject.rejthresh = []; res = com; end;
if ~isfield(EEG.reject, 'rejthreshE')		EEG.reject.rejthreshE = []; res = com; end;
if ~isfield(EEG.reject, 'rejfreq')			EEG.reject.rejfreq = []; res = com; end;
if ~isfield(EEG.reject, 'rejfreqE')			EEG.reject.rejfreqE = []; res = com; end;
if ~isfield(EEG.reject, 'rejconst')			EEG.reject.rejconst = []; res = com; end;
if ~isfield(EEG.reject, 'rejconstE')		EEG.reject.rejconstE = []; res = com; end;
if ~isfield(EEG.reject, 'icarejjp')			EEG.reject.icarejjp = []; res = com; end;
if ~isfield(EEG.reject, 'icarejjpE')		EEG.reject.icarejjpE = []; res = com; end;
if ~isfield(EEG.reject, 'icarejkurt')		EEG.reject.icarejkurt = []; res = com; end;
if ~isfield(EEG.reject, 'icarejkurtE')		EEG.reject.icarejkurtE = []; res = com; end;
if ~isfield(EEG.reject, 'icarejmanual')		EEG.reject.icarejmanual = []; res = com; end;
if ~isfield(EEG.reject, 'icarejmanualE')	EEG.reject.icarejmanualE = []; res = com; end;
if ~isfield(EEG.reject, 'icarejthresh')		EEG.reject.icarejthresh = []; res = com; end;
if ~isfield(EEG.reject, 'icarejthreshE')	EEG.reject.icarejthreshE = []; res = com; end;
if ~isfield(EEG.reject, 'icarejfreq')		EEG.reject.icarejfreq = []; res = com; end;
if ~isfield(EEG.reject, 'icarejfreqE')		EEG.reject.icarejfreqE = []; res = com; end;
if ~isfield(EEG.reject, 'icarejconst')		EEG.reject.icarejconst = []; res = com; end;
if ~isfield(EEG.reject, 'icarejconstE')		EEG.reject.icarejconstE = []; res = com; end;

if ~isfield(EEG.reject, 'rejglobal')		EEG.reject.rejglobal = []; res = com; end;
if ~isfield(EEG.reject, 'rejglobalE')		EEG.reject.rejglobalE = []; res = com; end;

if ~isfield(EEG, 'stats')			EEG.stats.jp = []; res = com; end;
if ~isfield(EEG.stats, 'jp')		EEG.stats.jp = []; res = com; end;
if ~isfield(EEG.stats, 'jpE')		EEG.stats.jpE = []; res = com; end;
if ~isfield(EEG.stats, 'icajp')		EEG.stats.icajp = []; res = com; end;
if ~isfield(EEG.stats, 'icajpE')	EEG.stats.icajpE = []; res = com; end;
if ~isfield(EEG.stats, 'kurt')		EEG.stats.kurt = []; res = com; end;
if ~isfield(EEG.stats, 'kurtE')		EEG.stats.kurtE = []; res = com; end;
if ~isfield(EEG.stats, 'icakurt')	EEG.stats.icakurt = []; res = com; end;
if ~isfield(EEG.stats, 'icakurtE')	EEG.stats.icakurtE = []; res = com; end;

% component rejection
% -------------------
if ~isfield(EEG.stats, 'compenta')		EEG.stats.compenta = []; res = com; end;
if ~isfield(EEG.stats, 'compentr')		EEG.stats.compentr = []; res = com; end;
if ~isfield(EEG.stats, 'compkurta')		EEG.stats.compkurta = []; res = com; end;
if ~isfield(EEG.stats, 'compkurtr')		EEG.stats.compkurtr = []; res = com; end;
if ~isfield(EEG.stats, 'compkurtdist')	EEG.stats.compkurtdist = []; res = com; end;
if ~isfield(EEG.reject, 'threshold')		EEG.reject.threshold = [0.8 0.8 0.8]; res = com; end;
if ~isfield(EEG.reject, 'threshentropy')	EEG.reject.threshentropy = 600; res = com; end;
if ~isfield(EEG.reject, 'threshkurtact')	EEG.reject.threshkurtact = 600; res = com; end;
if ~isfield(EEG.reject, 'threshkurtdist')	EEG.reject.threshkurtdist = 600; res = com; end;
if ~isfield(EEG.reject, 'gcompreject')		EEG.reject.gcompreject = []; res = com; end;
if length(EEG.reject.gcompreject) ~= size(EEG.icaact,1)
	EEG.reject.gcompreject = zeros(1, size(EEG.icaact,1));
end;
if isempty(EEG.reject.gcompreject)
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
		 case 'ica', 
		  if isempty(EEG.icaweights)
			  ButtonName=questdlg([ 'No ICA weights. Compute now?' 10 '(then go back to the function you just called)'], ...
								  'Confirmation', 'Cancel', 'Yes','Yes');
			  
			  switch lower(ButtonName),
			   case 'cancel', error('eeg_checkset: ICA components must be derived before running that function'); 
			  end;
			  [EEG res] = pop_runica(EEG);
			  res = [ inputnames(1) ' = eeg_checkset('  inputnames(1) '); ' res ];
		  else, return; end;
		 case 'epoch', 
		  if EEG.trials == 1
			  errordlg([ 'Epochs must be extracted before running that function' 10 'Use /Tools/Extract epochs'], 'Error');
			  error('eeg_checkset: epochs must be extracted before running that function');
		  end;
		 case 'event', 
		  if isempty(EEG.event)
			  errordlg([ 'Cannot process if no events. First add events.' 10 'Use /File/Import event info or /Import epoch info'], 'Error');
			  error('eeg_checkset: no events');
		  end;
		 case 'chanloc', 
		  if isempty(EEG.chanlocs)
			  errordlg( ['Cannot process without channel location file.' 10 ...
						 'Enter the name of the file via "/Edit/Edit dataset info".' 10 ...
						 'For the file format, enter ''>> help totoplot'' from the command line.' ], 'Error');
			  error('eeg_checkset: cannot process without channel location file.');
		  end;
		 case 'eventconsistency',	
		  if isempty(EEG.event), return; end;
		  
		  % remove the events which latency are out of boundary
		  % ---------------------------------------------------
		  if isfield(EEG.event, 'latency')
			  try, alllatencies = cell2mat( { EEG.event.latency } );
			  catch, error('Checkset: error empty latency entry for new events added by user');
			  end;
			  I1 = find(alllatencies < 0);
			  I2 = find(alllatencies > EEG.pnts*EEG.trials);
			  if (length(I1) + length(I2)) > 0 
				  fprintf('Checkset warning: %d/%d events had out-of-bounds latencies and were removed\n', ...
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
 			  	  allvalues = eval(['{ EEG.event.' difffield{index} ' };']);
 				  valempt = cellfun('isempty', allvalues);
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
					  fprintf(['Eeg_checkset: found empty values for field ''' difffield{index} '''\n']);
					  fprintf(['Eeg_checkset: filling with value of other events in the same epochs\n']);
				  end;
 			  end;
		  end;
		  
		  % uniformize fields (str or int) if necessary
		  % -------------------------------------------
		  allfields = fieldnames(EEG.event);
		  for index = 1:length(allfields)
			  eval(['allvalues = { EEG.event.' allfields{index} ' };']);
			  valreal = cellfun('isclass', allvalues, 'double');
			  
			  format = 'ok';
			  if ~all(valreal) % all valreal ok
				  format = 'str';
				  if all(valreal == 0) % all valreal=0 ok
					  format = 'ok';
				  end;
			  end;
			  if strcmp(format, 'str')
				  fprintf('eeg_checkset: uniformize value type of event field ''%s''\n', allfields{index});
				  % get the field content
				  % ---------------------
				  for indexevent = 1:length(EEG.event)
					  if valreal(indexevent)
						  EEG.event = setfield(EEG.event, { indexevent }, allfields{index}, num2str(allvalues{indexevent}) );
					  end;
				  end;
			  end;
		  end;
		 otherwise, error('eeg_checkset: unknown option');
        end;        
    end;
end;            

return;	

function num = popask( text )
	 ButtonName=questdlg( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;
