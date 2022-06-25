% eeg_boundarytype() - return boundary event. It is usually "boundary". 
%                Based on EEGLAB options, it can be -99.
% Usage:
%   >> boundaryType = eeg_boundarytype(INEEG);
%
% Inputs:
%   INEEG         - input EEG dataset structure. Used to check if event
%                   types are string or numerical.
%
% Outputs:
%   boundaryType - boundary event type. Usually "boundary" but based on EEGLAB
%                   preferences can be -99 if event types are numerical.
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

function boundaryType = eeg_boundarytype(INEEG1, INEEG2)

    if nargin < 1
        help eeg_boundarytype;
        return
    end

    if nargin < 2 
        INEEG2 = INEEG1;
    end

    if isfield(INEEG1, 'event') && isfield(INEEG1, 'setname') tmpevent1 = INEEG1.event; else tmpevent1 = INEEG1; end
    if isfield(INEEG2, 'event') && isfield(INEEG2, 'setname') tmpevent2 = INEEG2.event; else tmpevent2 = INEEG2; end

    % type of boundary event
    eeglab_options;
    boundaryType = 'boundary';
    cond1 = isempty(tmpevent1) || (isfield(tmpevent1, 'type') && isnumeric(tmpevent1(1).type));
    cond2 = isempty(tmpevent2) || (isfield(tmpevent2, 'type') && isnumeric(tmpevent2(1).type));
    if option_boundary99 && ~(isempty(tmpevent1) && isempty(tmpevent2))
        if cond1 && cond2
            boundaryType = -99;
        end
    end