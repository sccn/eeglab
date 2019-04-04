% getchanlist() - Obtain indices of specified channel types. 
%
% Usage:
%       >> chanlist = getchanlist(chanlocs, type)
%
% Inputs:
%     chanlocs - EEGLAB channel location structure
%     type     - [string] select channel of specified type
%                can enter a cell array to select several channel types
%
% Output:
%     chanlist - list of channel indices for the selected type(s)
%                sorted in increasing order
%
% Note: this function is not case sensitive
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2004
%
% See also: pop_chanedit()

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2004
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

function chanlist = getchanlist(chanlocs, type)
    
    if nargin < 1
        help getchanlist;
        return;
    end
    
    if nargin < 2 || ~isfield(chanlocs, 'type')
        chanlist = [1:length(chanlocs)];
        return;
    end
    
    % search channel types
    % --------------------
    if ischar(type), type = { type }; end
    type = lower(type);
    chanlist = [];
    alltypes = lower({ chanlocs.type });
    for index = 1:length(type)
        tmplist = strmatch(type{index}, alltypes, 'exact');
        if isempty(tmplist)
            fprintf('Warning: no channel of type ''%s'' found\n', type{index});
        end
        chanlist = [ chanlist tmplist' ];
    end
    chanlist = sort(chanlist);
    
    
