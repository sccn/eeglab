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

function chanlist = getchanlist(chanlocs, type)
    
    if nargin < 1
        help getchanlist;
        return;
    end;
    
    if nargin < 2 || ~isfield(chanlocs, 'type')
        chanlist = [1:length(chanlocs)];
        return;
    end;
    
    % search channel types
    % --------------------
    if isstr(type), type = { type }; end;
    type = lower(type);
    chanlist = [];
    alltypes = lower({ chanlocs.type });
    for index = 1:length(type)
        tmplist = strmatch(type{index}, alltypes, 'exact');
        if isempty(tmplist)
            fprintf('Warning: no channel of type ''%s'' found\n', type{index});
        end;
        chanlist = [ chanlist tmplist' ];
    end;
    chanlist = sort(chanlist);
    
    
