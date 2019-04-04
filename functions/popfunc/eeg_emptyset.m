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

function EEG = eeg_emptyset();

EEG.setname     = '';
EEG.filename    = '';
EEG.filepath    = '';
EEG.subject     = '';
EEG.group       = '';
EEG.condition   = '';
EEG.session     = [];
EEG.comments    = '';
EEG.nbchan      = 0;
EEG.trials      = 0;
EEG.pnts        = 0;
EEG.srate       = 1;
EEG.xmin        = 0;
EEG.xmax        = 0;
EEG.times       = [];
EEG.data        = [];
EEG.icaact      = [];
EEG.icawinv     = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.chanlocs    = [];
EEG.urchanlocs  = [];
EEG.chaninfo    = [];
EEG.ref         = [];
EEG.event       = [];
EEG.urevent     = [];
EEG.eventdescription = {};
EEG.epoch       = [];
EEG.epochdescription = {};
EEG.reject      = [];
EEG.stats       = [];
EEG.specdata    = [];
EEG.specicaact  = [];
EEG.splinefile  = '';
EEG.icasplinefile = '';
EEG.dipfit      = [];
EEG.history     = '';
EEG.saved       = 'no';
EEG.etc         = [];

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
