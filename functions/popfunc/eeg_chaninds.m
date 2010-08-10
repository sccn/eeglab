% std_chaninds() - look up channel indices in a EEG structure
%
% Usage:
%         >> inds = std_chaninds(EEG, channames);
% Inputs:
%         EEG   - EEG structure containing a chanlocs substructure.
%                 the chanlocs structure may also be used as input.
%     channames - [cell] channel names. May also be a string containing
%                 one or several channel names.
%
% Outputs:
%       inds - [integer array] channel indices
%
% Author: Arnaud Delorme, CERCO, 2009-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function finalinds = eeg_chaninds(EEG, channames, errorifnotfound);

    if nargin < 2
        help eeg_chaninds;
        return;
    end;
    if nargin < 3
        errorifnotfound = 1;
    end;
    
    if isfield(EEG, 'chanlocs')
         chanlocs = EEG.chanlocs;
    else chanlocs = EEG;
    end;
    
    % decode string if necessary
    % --------------------------
    if isstr(channames)
        channames = parsetxt( channames );
    end;

    finalinds   = [];
    if isempty(chanlocs)
         tmpallchans = [];
    else tmpallchans = lower({ chanlocs.labels });
    end;
    if isempty(channames), finalinds = [1:length(chanlocs)]; return; end;
    for c = 1:length(channames)
        chanind = strmatch( lower(channames{c}), tmpallchans, 'exact');
        if isempty(chanind), 
            chanind = str2double(channames{c});
            if isnan(chanind), chanind = []; end;
            if errorifnotfound && isempty(chanind)
                error(sprintf('Channel %s not found', channames{c})); 
            end;
        end;
        finalinds   = [ finalinds chanind ];
    end;
    finalinds = sort(finalinds);
