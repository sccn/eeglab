% eeg_cache() - Store data in cache with hashcode
%
% >> [cache data] = eeg_cache(cache, hashcode, data);  
%
% Inputs:
%   cache      - Cell array containing cached data
%   hashcode   - MD5 hashcode
%
% Optional input:
%   data       - This can be anything
%
% Outputs:
%   cache      - Updated cache information
%   data       - Data corresponding to the hash code. If the hash code
%                is not found, this output is empty.
%
% Note: this script checks if the cache variable is larger than option_cachesize
% If it is the case, elements are removed until an acceptable size is
% reached.
%
% Examples
% cache = eeg_cache([], 'dyufisf8da0df', 3); % store 3 with hashcode dyufisf8da0df
% [cache data] = eeg_cache(cache, 'dyufisf8da0df'); % retreive 3
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, 2015

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, June 07, 2007, arno@sccn.ucsd.edu
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

function [cache data] = eeg_cache(cache, hashcode, data)

eeglab_options;

ind = [];
if ~isempty(cache)
    allHashCodes = { cache.hashcode };
    maxlen       = max(max(cellfun(@length, allHashCodes)), length(hashcode));
    ind = find(strncmp(hashcode, allHashCodes, maxlen));
end
    
if isempty(ind)
    if nargin > 2
        cache(end+1).hashcode = hashcode;
        cache(end  ).data     = data;

        % make sure the cache is not too large
        while getfield(whos('cache'), 'bytes') > option_cachesize*1000000
            cache(1) = [];
        end
    else
        data = [];
    end
else
    data = cache(ind).data;
end
