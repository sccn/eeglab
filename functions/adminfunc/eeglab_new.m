  % EEGLAB_NEW - script called when a new dataset is created and require 
%              input from user. Use workspace variables EEG, ALLEEG,
%              CURRENTSET, LASTCOM
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2022
%
% See also: EEGLAB

% Copyright (C) 2022 Arnaud Delorme
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

% functions below modify the current dataset, but should not create new datasets
functionsEEGLAB = { 'pop_chanedit' 'pop_editset' 'pop_comment' 'pop_headplot' 'pop_selectcomps' ...
              'pop_runica'   'pop_editeventfield' 'pop_editeventvals' 'pop_adjustevents' ...
              'pop_saveset'  'pop_importepoch'    'pop_importevent'   'pop_chanevent' ...
              'pop_importpres' 'pop_importerplab' 'pop_iclabel' 'pop_icflag' 'pop_dipfit_headmodel' ...
              'pop_dipfit_settings' 'pop_dipfit_gridsearch' 'pop_dipfit_nonlinear' 'pop_multifit' 'pop_leadfield' 'topoplot' ';;;' ...
              'pop_roi_connect' 'pop_importmff' 'pop_roi_activity' 'pop_eegstats' };

posEqual = find('=' == LASTCOM);
if ~isempty(LASTCOM) && ~isempty(EEG) ...
        && ~isempty(posEqual) ...
        && contains(LASTCOM(1:posEqual(1)), 'EEG') ... % when no equal sign found, do not create new datasets
        && ~contains(LASTCOM, functionsEEGLAB) ... % when using one of the functions in the list, do not create new datasets
        && ~exist('DEBUG_EEGLAB_MENUS', 'var') % when in debug menu mode, do not create new datasets (see eeglab_execmenu)

    % Here the dataset was modified 
    EEG = eegh(LASTCOM, EEG); 
    [ALLEEG, EEG, CURRENTSET, LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, 'study', ~isempty(STUDY)+0);
    eegh(LASTCOM);
    disp('Done');
    eeglab('redraw');

elseif ~isempty(LASTCOM)

    % Here the dataset was not modified (plotting only) 
    EEG = eegh(LASTCOM, EEG); 
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    if ~isempty(posEqual) && contains(LASTCOM(1:posEqual(1)), 'EEG')
        % dataset modified, store in history
        eegh('[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);');
    end
    eeglab('redraw');
    % even though EEG is modified by eegh call, do not keep eeg_store in history
    % so the history does not reflect modifications to EEG.history for
    % better redeability

end

clear functionsEEGLAB posEqual;