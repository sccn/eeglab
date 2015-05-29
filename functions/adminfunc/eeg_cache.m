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

function [cache data] = eeg_cache(cache, hashcode, data)

eeglab_options;

ind = [];
if ~isempty(cache)
    allHashCodes = { cache.hashcode };
    maxlen       = max(max(cellfun(@length, allHashCodes)), length(hashcode));
    ind = find(strncmp(hashcode, allHashCodes, maxlen));
end;
    
if isempty(ind)
    if nargin > 2
        cache(end+1).hashcode = hashcode;
        cache(end  ).data     = data;

        % make sure the cache is not too large
        while getfield(whos('cache'), 'bytes') > option_cachesize
            cache(1) = [];
        end;
    else
        data = [];
    end;
else
    data = cache(ind).data;
end;
