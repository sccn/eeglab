% eeg_isboundary() - check if an event is a boundary event
%
% Usage:
%   >> logical = eeg_boundaryevent(INEEG);
%   >> logical = eeg_isboundary(event);
%
% Inputs:
%   event         - EEGLAB event
%
% Outputs:
%   logical - true or false
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

function logi = eeg_isboundary(tmpevent)

    if nargin < 1
        help eeg_isboundary;
        return
    end

    logi = false;
    if isempty(tmpevent) || ~isfield(tmpevent, 'type')
        return;
    else
        % type of boundary event
        eeglab_options;
        if ischar(tmpevent(1).type)
            if strcmpi('boundary', tmpevent(1).type)
                logi = true;
            end
        elseif option_boundary99
            if tmpevent(1).type == -99
                logi = true;
            end
        end
    end