% eeg_emptyset() - Initialize an EEG dataset structure with default values.
%
% Usage:
%   >> EEG = eeg_emptyset();
%
% Outputs:
%   EEG    - empty dataset structure with default values.
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
% Revision 1.7  2004/03/02 16:37:49  arno
% reference unspecified
%
% Revision 1.6  2003/07/28 15:22:03  arno
% averef -> ref
%
% Revision 1.5  2002/09/25 23:52:41  arno
% average referece
%
% Revision 1.4  2002/08/08 18:08:57  arno
% removing rt
%
% Revision 1.3  2002/06/25 14:37:26  arno
% header problem fix
%
% Revision 1.2  2002/04/08 20:47:54  arno
% adding comments field
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function EEG = eeg_emptyset();

EEG.setname    = '';
EEG.filename   = '';
EEG.filepath   = '';
EEG.pnts       = 0;
EEG.nbchan     = 0;
EEG.trials     = 0;
EEG.srate      = 1;
EEG.xmin       = 0;
EEG.xmax       = 0;
EEG.data       = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact    = [];
EEG.event     = [];
EEG.epoch  = [];
EEG.chanlocs    = '';
EEG.chaninfo    = '';
EEG.comments    = '';
EEG.ref         = []; % unspecified

%EEG.reject.threshold  = [1 0.8 0.85];
%EEG.reject.icareject  = [];
%EEG.reject.compreject = [];
%EEG.reject.gcompreject= [];
%EEG.reject.comptrial  = [];
%EEG.reject.sigreject  = [];
%EEG.reject.elecreject = [];

%EEG.stats.kurta      = [];
%EEG.stats.kurtr      = [];
%EEG.stats.kurtd      = [];		
%EEG.stats.eegentropy = [];
%EEG.stats.eegkurt    = [];
%EEG.stats.eegkurtg   = [];
%EEG.stats.entropy    = [];
%EEG.stats.kurtc      = [];
%EEG.stats.kurtt      = [];
%EEG.stats.entropyc   = [];

return;
