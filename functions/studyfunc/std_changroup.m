% std_changroup() - Create channel groups for plotting.
%
% Usage:    
%                >> STUDY = std_changroup(STUDY, ALLEEG);   
%                >> STUDY = std_changroup(STUDY, ALLEEG, chanlocs, 'interp');   
% Inputs:
%   ALLEEG     - Top-level EEGLAB vector of loaded EEG structures for the dataset(s) 
%                in the STUDY. ALLEEG for a STUDY set is typically loaded using 
%                pop_loadstudy(), or in creating a new STUDY, using pop_createstudy().  
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   chanlocs   - EEGLAB channel structure. Only construct the STUDY.changrp
%                structure for a subset of channels.
%   'interp'   - optional input in case channel locations are interpolated
%
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user 
%                edits, if any. The STUDY.changrp structure is created. It contains as
%                many elements as there are channels. For example, STUDY.changrp(1)
%                is the first channel. Fields of the changrp structure created at this
%                point are 
%                    STUDY.changrp.name      : name of the channel group
%                    STUDY.changrp.channels  : cell array containing channel labels
%                                              for the group.
%                    STUDY.changrp.setinds   : indices of datasets containing the
%                                              selected channels.
%                    STUDY.changrp.allinds   : indices of channels within the datasets 
%                                              above.
%
% Authors: Arnaud Delorme, CERCO, 2006

% Copyright (C) Arnaud Delorme, CERCO, arno@salk.edu
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

function STUDY = std_changroup(STUDY, ALLEEG, alllocs, interp);

if nargin < 4
    interp = 'off';
end

% union of all channel structures
% -------------------------------
inputloc = 0;
if nargin >= 3
    if ~isempty(alllocs)
        inputloc = 1;
    end
end
if ~inputloc
    alllocs = eeg_mergelocs(ALLEEG.chanlocs);
end

% create group for each electrode
% -------------------------------
if isstruct(alllocs)
    alllocs = { alllocs.labels };
end
STUDY.changrp = [];
for indc = 1:length(alllocs)
    STUDY.changrp(indc).name     = alllocs{indc};
    STUDY.changrp(indc).channels = { alllocs{indc} };
end
