% pop_loadset() - load a dataset (pop out window if no arguments)
%
% Usage:
%   >> [EEG] = pop_loadset( filename, filepath);
%
% Inputs:
%   filename  - file name
%   filepath  - file path
%
% Outputs:
%   EEG       - data structure
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
% Revision 1.2  2002/04/08 23:33:26  arno
% format conversion for events
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, command] = pop_loadset( inputname, inputpath);
eeg_consts;

command = '';
%eeg_emptyset;
if nargin < 2
   
	[filename, filepath] = uigetfile('*.set', 'Load a dataset -- pop_loadset()');
	if filename == 0 return; end;
	inputname = filename;
	inputpath = filepath;

end;

fprintf('loading set %s ...\n', inputname);
eval( [ 'load -mat ' inputpath inputname ] );

if ~isfield( EEG, 'data')
	fprintf('Incompatible with new format, trying old format...\n');
   eegset = cellArray;

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
	EEG.epoch  = [ EEG.rt(:) EEG.eegtype(:) EEG.eegresp(:) EEG.accept(:) ];
	I = find( EEG.rt < 1000 ); 	
	EEG.event      = [ ones(1, length(I)); EEG.rt(I) ]'; % type is 1, meaing reaction time
	EEG = rmfield( EEG, {'accept', 'eegtype', 'eegresp', 'accept' });
	
	if EEG.trials > 1
		if size(EEG.event,1) == EEG.trials
			EEG.event = eeg_eventformat([  EEG.event [1:EEG.trials]'], 'struct', {'type', 'latency' 'epoch' });
		else
			EEG.event = eeg_eventformat( EEG.event, 'struct', {'type', 'latency'});
		end;
	else
		EEG.event = eeg_eventformat([ ones(length(EEG.trials)) EEG.event], 'struct', {'type', 'latency'});
	end;
	
	EEG.epoch = eeg_epochformat(EEG.epoch, 'struct', { 'rt' 'type' 'resp' 'accept' });
end;

% check modified fields
% ---------------------
if isfield(EEG,'icadata')
    EEG.icaact = EEG.icadata;
    EEG = rmfield(EEG, 'icadata');
end;  
if isfield(EEG,'poschan')
    EEG.chanlocs = EEG.poschan;
    EEG = rmfield(EEG, 'poschan');
end;  
EEG = eeg_checkset(EEG, 'eventconsistency');
command = sprintf('EEG = pop_loadset( ''%s'', ''%s'');', inputname, inputpath);
return;
