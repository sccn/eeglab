% std_chaninds() - look up channel indices in a STUDY
%
% Usage:
%         >> inds = std_chaninds(STUDY,  channames);
%         >> inds = std_chaninds(EEG, channames);
%         >> inds = std_chaninds(chanlocs, channames);
% Inputs:
%         STUDY    - studyset structure containing a changrp substructure.
%         EEG      - EEG structure containing channel location structure
%         chanlocs - channel location structure
%     channames - [cell] channel names
%
% Outputs:
%       inds - [integer array] channel indices
%
% Author: Arnaud Delorme, CERCO, 2006-

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

function finalinds = std_chaninds(instruct, channames);

    finalinds   = [];
    if isfield(instruct, 'chanlocs')
        EEG = instruct;
        tmpchanlocs = EEG.chanlocs;
        tmpallchans = lower({ tmpchanlocs.labels });
    elseif isfield(instruct, 'filename')
        tmpallchans = lower({ instruct.changrp.name });
    else
        tmpallchans = instruct;
    end;
    if ~iscell(channames), channames = { channames }; end;
    
    if isempty(channames), finalinds = [1:length(tmpallchans)]; return; end;
    for c = 1:length(channames)
        if isnumeric(channames)
            chanind = channames(c);
        else
            chanind = strmatch( lower(channames{c}), tmpallchans, 'exact');
            if isempty(chanind), warning(sprintf([ 'Channel "%s" and maybe others was not' 10 'found in pre-computed data file' ], channames{c})); end;
        end;
        finalinds   = [ finalinds chanind ];
    end;
