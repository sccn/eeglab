% pop_loadset() - load a dataset (pop out window if no arguments)
%
% Usage:
%   >> VAROUT = pop_loadset( filename, filepath);
%
% Inputs:
%   filename  - file name
%   filepath  - file path
%
% Output
%   VAROUT - variable read
%
% Note: output may be a dataset structure of an array of
%       dataset structures.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), pop_saveset()

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
% Revision 1.14  2002/10/14 23:04:02  arno
% default output
%
% Revision 1.13  2002/10/08 15:58:51  arno
% debugging .fdt files
%
% Revision 1.12  2002/08/22 01:41:31  arno
% further checks
%
% Revision 1.11  2002/08/22 00:44:59  arno
% compatibility with previous eeglab version
%
% Revision 1.10  2002/08/22 00:04:44  arno
% old format compatibility
%
% Revision 1.9  2002/08/11 19:20:36  arno
% removing eeg_consts
%
% Revision 1.8  2002/04/23 18:00:17  arno
% modifying help message
%
% Revision 1.7  2002/04/23 17:58:40  arno
% standalone function
%
% Revision 1.6  2002/04/12 00:52:53  arno
% smart load for non-identical structures
%
% Revision 1.5  2002/04/11 23:26:27  arno
% changing structure check
%
% Revision 1.4  2002/04/11 03:34:26  arno
% fully synthesis between pop_loadset and pop_loadwks, pop_loadwks removed
%
% Revision 1.3  2002/04/09 20:58:47  arno
% adding check of event consistency
%
% Revision 1.2  2002/04/08 23:33:26  arno
% format conversion for events
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [VAROUT, command] = pop_loadset( inputname, inputpath);

command = '';
VAROUT  = [];
%eeg_emptyset;
if nargin < 2  
	[inputname, inputpath] = uigetfile('*.set*', 'Load dataset(s) -- pop_loadset()');
    drawnow;
	if inputname == 0 return; end;
end;

fprintf('Pop_loadset: loading file %s ...\n', inputname);
TMPVAR = load([ inputpath inputname ], '-mat');

if isfield(TMPVAR, 'EEG') %individual dataset
	% load individual dataset
	% -----------------------
	VAROUT = checkoldformat(TMPVAR.EEG);
    if isstr(VAROUT.data), VAROUT.filepath = inputpath; end;
elseif isfield(TMPVAR, 'ALLEEG') %multiple dataset
	disp('Pop_loadset: appending datasets');
	VAROUT = TMPVAR.ALLEEG;
    for index=1:length(VAROUT)
        if isstr(VAROUT(index).data), VAROUT(index).filepath = inputpath; end;
    end;
else
	VAROUT = checkoldformat(TMPVAR);
	if ~isfield( VAROUT, 'data')
		error('Pop_loadset: non-EEGLAB dataset file');
	end;
    if isstr(VAROUT.data), VAROUT.filepath = inputpath; end;
end;
command = sprintf('EEG = pop_loadset( ''%s'', ''%s'');', inputname, inputpath);
return;

function EEG = checkoldformat(EEG)
	if ~isfield( EEG, 'data')
		fprintf('Incompatible with new format, trying old format and converting...\n');
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
		%	disp('Warning: some variable may not have been assigned');
		%end;
		
		% modify the eegtype to match the new one
		
		try, 
			if EEG.trials > 1
				EEG.events  = [ EEG.rt(:) EEG.eegtype(:) EEG.eegresp(:) ];
			end;
		catch, end;
	end;
	% check modified fields
	% ---------------------
	if isfield(EEG,'icadata'), EEG.icaact = EEG.icadata; end;  
	if isfield(EEG,'poschan'), EEG.chanlocs = EEG.poschan; end;  
	if ~isfield(EEG, 'icaact'), EEG.icaact = []; end;
	if ~isfield(EEG, 'chanlocs'), EEG.chanlocs = []; end;
	
	if isfield(EEG, 'events') & ~isfield(EEG, 'event')
		try, 
			if EEG.trials > 1
				EEG.events  = [ EEG.rt(:) ];
				
				EEG = eeg_checkset(EEG);;
				EEG = pop_importepoch(EEG, EEG.events, { 'rt'}, {'rt'}, 1E-3);
			end;
			if isfield(EEG, 'trialsval')
				EEG = pop_importepoch(EEG, EEG.trialsval(:,2:3), { 'eegtype' 'response' }, {},1,0,0);
			end;
			EEG = eeg_checkset(EEG, 'eventconsistency');
		catch, disp('warning: could not import events'); end;			
	end;
	rmfields = {'icadata' 'events' 'accept' 'eegtype' 'eegresp' 'trialsval' 'poschan' 'icadata' 'namechan' };
	for index = 1:length(rmfields)
		if isfield(EEG, rmfields{index}), 
			disp(['Warning: field ' rmfields{index} ' is deprecated']);
		end;
	end;
	
	