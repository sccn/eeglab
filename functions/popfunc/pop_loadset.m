% pop_loadset() - load an EEG dataset. If no arguments, pop up an input window.
%
% Usage:
%   >> EEGOUT = pop_loadset; % pop up window to input arguments
%   >> EEGOUT = pop_loadset( 'key1', 'val1', 'key2', 'val2', ...);
%   >> EEGOUT = pop_loadset( filename, filepath); % old calling format
%
% Optional inputs:
%   'filename'  - [string] dataset filename. Default pops up a graphical
%                 interface to browse for a data file.
%   'filepath'  - [string] dataset filepath. Default is current folder. 
%   'loadmode'  - ['all', 'info', integer] 'all' -> load the data and
%                 the dataset structure. 'info' -> load only the dataset 
%                 structure but not the actual data. [integer] ->  load only 
%                 a specific channel. This is efficient when data is stored 
%                 in a separate '.dat' file in which individual channels 
%                 may be loaded independently of each other. {default: 'all'}
%   'eeg'       - [EEG structure] reload current dataset
% Note:
%       Multiple filenames and filepaths may be specified. If more than one,
%       the output EEG variable will be an array of EEG structures.
% Output
%   EEGOUT - EEG dataset structure or array of structures
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001; SCCN/INC/UCSD, 2002-
%
% See also: eeglab(), pop_saveset()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 

function [EEG, command] = pop_loadset( inputname, inputpath, varargin)

command = '';
EEG  = [];

if nargin < 1
    % pop up window
    % -------------
	[inputname, inputpath] = uigetfile2('*.SET*;*.set', 'Load dataset(s) -- pop_loadset()', 'multiselect', 'on');
    drawnow;
	if isequal(inputname, 0) return; end
    options = { 'filename' inputname 'filepath' inputpath };
else
    % account for old calling format
    % ------------------------------
    if ~strcmpi(inputname, 'filename') && ~strcmpi(inputname, 'filepath') && ~strcmpi(inputname, 'eeg') && ~strcmpi(inputname, 'loadmode') && ~strcmpi(inputname, 'check')
        options = { 'filename' inputname }; 
        if nargin > 1
            options = { options{:} 'filepath' inputpath }; 
        end
        if nargin > 2
            options = { options{:} 'loadmode' varargin{1} }; 
        end
    else
        options = { inputname inputpath varargin{:} };
    end
end

% decode input parameters
% -----------------------
g = finputcheck( options, ...
                 { 'filename'   { 'string';'cell' }    []   '';
                   'filepath'   'string'               []   '';
                   'check'      'string'               { 'on';'off' }   'on';
                   'verbose'    'string'               { 'on';'off' }   'on';
                   'loadmode'   { 'string';'integer' } { { 'info' 'all' } [] }  'all';
                   'eeg'        'struct'               []   struct('data',{}) }, 'pop_loadset');
if ischar(g), error(g); end
if ischar(g.filename), g.filename = { g.filename }; end

% reloading EEG structure from disk
% ---------------------------------
if ~isempty(g.eeg)

    EEG = pop_loadset( 'filepath', g.eeg.filepath, 'filename', g.eeg.filename);

else
    eeglab_options;
    ALLEEGLOC = [];
    for ifile = 1:length(g.filename)
        
         if ifile > 1 && option_storedisk
              g.loadmode = 'last';
%             warndlg2(strvcat('You may only load a single dataset','when selecting the "Store at most one', 'dataset in memory" option'));
%             break;
         end
        
        % read file
        % ---------
        filename = fullfile(g.filepath, g.filename{ifile});
        if strcmpi(g.verbose, 'on')
            fprintf('pop_loadset(): loading file %s ...\n', filename);
        end
        if strcmpi(g.loadmode, 'info')
            if ismatlab
                TMPVAR = load('-mat', filename, '-regexp', '^((?!data).)*$');
            else
                TMPVAR = load('-mat', filename);
                if isfield(TMPVAR, 'data')
                    TMPVAR = rmfield(TMPVAR, 'data');
                end
            end
            if isfield(TMPVAR, 'setname')
                TMPVAR.data = 'in set file';
            end
        else
            TMPVAR = load('-mat', filename);
        end

        % variable not found
        % ------------------
        if isempty(TMPVAR)
            error('No dataset info is associated with this file');
        end

        if isfield(TMPVAR, 'EEG')

            % load individual dataset
            % -----------------------
            EEG = checkoldformat(TMPVAR.EEG);
            [ EEG.filepath EEG.filename ext ] = fileparts( filename );
            EEG.filename = [ EEG.filename ext ];

            % account for name changes etc...
            % -------------------------------
            if ischar(EEG.data) && ~strcmpi(EEG.data, 'EEGDATA')

                [tmp EEG.data ext] = fileparts( EEG.data ); EEG.data = [ EEG.data ext];
                if ~isempty(tmp) && ~strcmpi(tmp, EEG.filepath)
                    disp('Warning: updating folder name for .dat|.fdt file');
                end
                if ~strcmp(EEG.filename(1:end-3), EEG.data(1:end-3))
                    disp('Warning: the name of the dataset has changed on disk, updating EEG structure accordingly');
                    EEG.data    = [ EEG.filename(1:end-3) EEG.data(end-2:end) ];
                    EEG.datfile = [ EEG.filename(1:end-3) EEG.data(end-2:end) ];
                    EEG.saved = 'no';
                end

            end

            % copy data to output variable if necessary (deprecated)
            % -----------------------------------------
            if ~strcmpi(g.loadmode, 'info') && isfield(TMPVAR, 'EEGDATA')
                if ~option_storedisk || ifile == length(g.filename)
                    EEG.data = TMPVAR.EEGDATA;
                end
            end

        elseif isfield(TMPVAR, 'ALLEEG') % old format

            eeglab_options;
            if option_storedisk
                error('Cannot load multiple dataset file. Change memory option to allow multiple datasets in memory, then try again. Remember that this file type is OBSOLETE.');
            end

            % this part is deprecated as of EEGLAB 5.00
            % since all dataset data have to be saved in separate files
            % -----------------------------------------------------
            disp('pop_loadset(): appending datasets');
            EEG = TMPVAR.ALLEEG;
            for index=1:length(EEG)
                EEG(index).filename = '';
                EEG(index).filepath = '';        
                if ischar(EEG(index).data), 
                    EEG(index).filepath = g.filepath; 
                    if length(g.filename{ifile}) > 4 && ~strcmp(g.filename{ifile}(1:end-4), EEG(index).data(1:end-4)) && strcmpi(g.filename{ifile}(end-3:end), 'sets')
                        disp('Warning: the name of the dataset has changed on disk, updating .dat data file to the new name');
                        EEG(index).data = [ g.filename{ifile}(1:end-4) 'fdt' int2str(index) ];
                    end
                end
            end
        else
            EEG = checkoldformat(TMPVAR);
            if ~isfield( EEG, 'data')
                error('pop_loadset(): not an EEG dataset file');
            end
            if ischar(EEG.data), EEG.filepath = g.filepath; end
        end
        
        %ALLEEGLOC = pop_newset(ALLEEGLOC, EEG, 1);
        ALLEEGLOC = eeg_store(ALLEEGLOC, EEG, 0, 'verbose', 'off');
                
    end
    EEG = ALLEEGLOC;
end

% load all data or specific data channel
% --------------------------------------
if strcmpi(g.check, 'on')
    EEG = eeg_checkset(EEG);
end
if ischar(g.loadmode)
    if strcmpi(g.loadmode, 'all')
        EEG = eeg_checkset(EEG, 'loaddata');
    elseif strcmpi(g.loadmode, 'last')
        EEG(end) = eeg_checkset(EEG(end), 'loaddata');
    end
else
    % load/select specific channel
    % ----------------------------
    EEG.datachannel = g.loadmode;
    EEG.data   = eeg_getdatact(EEG, 'channel', g.loadmode);
    EEG.nbchan = length(g.loadmode);
    if ~isempty(EEG.chanlocs)
        EEG.chanlocs = EEG.chanlocs(g.loadmode);
    end
    EEG.icachansind = [];
    EEG.icaact = [];
    EEG.icaweights = [];
    EEG.icasphere = [];
    EEG.icawinv = [];
    %if ischar(EEG.data)
    %    EEG.datfile = EEG.data;
    %    fid = fopen(fullfile(EEG.filepath, EEG.data), 'r', 'ieee-le');
    %    fseek(fid, EEG.pnts*EEG.trials*( g.loadmode - 1), 0 );
    %    EEG.data    = fread(fid, EEG.pnts*EEG.trials, 'float32');
    %    fclose(fid);
    %else
    %    EEG.data        = EEG.data(g.loadmode,:,:);
    %end
end

% set file name and path
% ----------------------
if length(EEG) == 1
    tmpfilename = g.filename{1};
    if isempty(g.filepath)
        [g.filepath tmpfilename ext] = fileparts(tmpfilename);
        tmpfilename = [ tmpfilename ext ];
    end
    EEG.filename = tmpfilename;
    EEG.filepath = g.filepath;
end

% set field indicating that the data has not been modified
% --------------------------------------------------------
if isfield(EEG, 'changes_not_saved')
    EEG = rmfield(EEG, 'changes_not_saved');
end
for index=1:length(EEG)
    EEG(index).saved = 'justloaded';
end

command = sprintf('EEG = pop_loadset(%s);', vararg2str(options));
return;

function EEG = checkoldformat(EEG)
	if ~isfield( EEG, 'data')
		fprintf('pop_loadset(): Incompatible with new format, trying old format and converting...\n');
		eegset = EEG.cellArray;
		
		off_setname             = 1;  %= filename
		off_filename            = 2;  %= filename
		off_filepath            = 3;  %= fielpath
		off_type                    = 4;  %= type EEG AVG CNT
		off_chan_names          = 5;  %= chan_names
		off_chanlocs            = 21;  %= filename
		off_pnts                    = 6;  %= pnts
		off_sweeps                  = 7; %= sweeps
		off_rate                    = 8;  %= rate
		off_xmin                    = 9;  %= xmin
		off_xmax                    = 10;  %= xmax
		off_accept                  = 11; %= accept
		off_typeeeg                 = 12; %= typeeeg
		off_rt                      = 13; %= rt
		off_response            = 14; %= response
		off_signal                  = 15; %= signal
		off_variance            = 16; %= variance
		off_winv                    = 17; %= variance
		off_weights             = 18; %= variance
		off_sphere                  = 19; %= variance
		off_activations         = 20; %= variance
		off_entropytrial        = 22; %= variance
		off_entropycompo        = 23; %= variance
		off_threshold       = 24; %= variance
		off_comporeject     = 25; %= variance
		off_sigreject       = 26;
		off_kurtA                       = 29;
		off_kurtR                       = 30;
		off_kurtDST                 = 31;
		off_nbchan                  = 32;
		off_elecreject      = 33;
		off_comptrial       = 34;
		off_kurttrial       = 35; %= variance
		off_kurttrialglob   = 36; %= variance
		off_icareject       = 37; %= variance
		off_gcomporeject    = 38; %= variance
		off_eegentropy          = 27;
		off_eegkurt             = 28;
		off_eegkurtg            = 39;

		off_tmp1                        = 40;
		off_tmp2                        = 40;
		
		% must convert here into new format
		EEG.setname    = eegset{off_setname   };
		EEG.filename   = eegset{off_filename  };
		EEG.filepath   = eegset{off_filepath  };
		EEG.namechan   = eegset{off_chan_names};
		EEG.chanlocs    = eegset{off_chanlocs   };
		EEG.pnts       = eegset{off_pnts      };
		EEG.nbchan     = eegset{off_nbchan    };
		EEG.trials     = eegset{off_sweeps    };
		EEG.srate       = eegset{off_rate      };
		EEG.xmin       = eegset{off_xmin      };
		EEG.xmax       = eegset{off_xmax      };
		EEG.accept     = eegset{off_accept    };
		EEG.eegtype    = eegset{off_typeeeg   };
		EEG.rt         = eegset{off_rt        };
		EEG.eegresp    = eegset{off_response  };
		EEG.data     = eegset{off_signal    };
		EEG.icasphere  = eegset{off_sphere    };
		EEG.icaweights = eegset{off_weights   };
		EEG.icawinv       = eegset{off_winv      };
		EEG.icaact        = eegset{off_activations  };
		EEG.stats.entropy    = eegset{off_entropytrial };
		EEG.stats.kurtc      = eegset{off_kurttrial    };
		EEG.stats.kurtg      = eegset{off_kurttrialglob};
		EEG.stats.entropyc   = eegset{off_entropycompo };
		EEG.reject.threshold  = eegset{off_threshold    };
		EEG.reject.icareject  = eegset{off_icareject    };
		EEG.reject.compreject = eegset{off_comporeject  };
		EEG.reject.gcompreject= eegset{off_gcomporeject };
		EEG.reject.comptrial  = eegset{off_comptrial    };
		EEG.reject.sigreject  = eegset{off_sigreject    };
		EEG.reject.elecreject = eegset{off_elecreject   };
		EEG.stats.kurta      = eegset{off_kurtA        };
		EEG.stats.kurtr      = eegset{off_kurtR        };
		EEG.stats.kurtd      = eegset{off_kurtDST      };
		EEG.stats.eegentropy = eegset{off_eegentropy   };
		EEG.stats.eegkurt    = eegset{off_eegkurt      };
		EEG.stats.eegkurtg   = eegset{off_eegkurtg     };
		%catch
		%	disp('Warning: some variables may not have been assigned');
		%end
		
		% modify the eegtype to match the new one
		
		try
			if EEG.trials > 1
				EEG.events  = [ EEG.rt(:) EEG.eegtype(:) EEG.eegresp(:) ];
			end
		catch, end
	end
	% check modified fields
	% ---------------------
	if isfield(EEG,'icadata'), EEG.icaact = EEG.icadata; end;  
	if isfield(EEG,'poschan'), EEG.chanlocs = EEG.poschan; end;  
	if ~isfield(EEG, 'icaact'), EEG.icaact = []; end
	if ~isfield(EEG, 'chanlocs'), EEG.chanlocs = []; end
	
	if isfield(EEG, 'events') && ~isfield(EEG, 'event')
		try
			if EEG.trials > 1
				EEG.events  = [ EEG.rt(:) ];
				
				EEG = eeg_checkset(EEG);
				EEG = pop_importepoch(EEG, EEG.events, { 'rt'}, {'rt'}, 1E-3);
			end
			if isfield(EEG, 'trialsval')
				EEG = pop_importepoch(EEG, EEG.trialsval(:,2:3), { 'eegtype' 'response' }, {},1,0,0);
			end
			EEG = eeg_checkset(EEG, 'eventconsistency');
		catch, disp('Warning: could not import events'); end;			
	end
	rmfields = {'icadata' 'events' 'accept' 'eegtype' 'eegresp' 'trialsval' 'poschan' 'icadata' 'namechan' };
	for index = 1:length(rmfields)
		if isfield(EEG, rmfields{index}), 
			disp(['Warning: field ' rmfields{index} ' is deprecated']);
		end
	end
