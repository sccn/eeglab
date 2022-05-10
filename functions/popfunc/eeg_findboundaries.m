% eeg_findboundaries() - return indices of boundary events
%
% Usage:
%   >> boundaryIndices = eeg_boundaryevent(INEEG);
%   >> boundaryIndices = eeg_boundaryevent(INEEG.event);
%
% Inputs:
%   INEEG         - input EEG dataset structure. Used to check if event
%                   types are string or numerical.
%
% Outputs:
%   boundaryIndices - indices of boundary events
%
% Author: Arnaud Delorme, 2022
% 
% see also: eeglab()

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

function boundaries = eeg_findboundaries(EEG)

    if nargin < 1
        help eeg_findboundaries;
        return
    end

    if isempty(EEG)
        boundaries = [];
        return;
    end

    if isfield(EEG, 'event') && isfield(EEG, 'setname')
        tmpevent = EEG.event;
    else
        tmpevent = EEG;
    end

    % type of boundary event
    eeglab_options;
    if ischar(tmpevent(1).type)
        boundaries = strmatch('boundary', { tmpevent.type });
    elseif option_boundary99
        boundaries = find([ tmpevent.type ] == -99);
    else
        boundaries = [];
    end
