% std_changroup() - Create channel groups for plotting.
%
% Usage:    
%                >> STUDY = std_changroup(STUDY, ALLEEG);   
%                >> STUDY = std_changroup(STUDY, ALLEEG, chanlocs, 'interp');   
% Inputs:
%   ALLEEG     - Top-level EEGLAB vector of loaded EEG structures for the dataset(s) 
%                in the STUDY. ALLEEG for a STUDY set is typically loaded using 
%                pop_loadstudy(), or in creating a new STUDY, using pop_createstudy().  
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   chanlocs   - EEGLAB channel structure. Only construct the STUDY.changrp
%                structure for a subset of channels.
%   'interp'   - optional input in case channel locations are interpolated
%
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user 
%                edits, if any. The STUDY.changrp structure is created. It contains as
%                many elements as there are channels. For example, STUDY.changrp(1)
%                is the first channel. Fields of the changrp structure created at this
%                point are 
%                    STUDY.changrp.name      : name of the channel group
%                    STUDY.changrp.channels  : cell array containing channel labels
%                                              for the group.
%                    STUDY.changrp.setinds   : indices of datasets containing the
%                                              selected channels.
%                    STUDY.changrp.allinds   : indices of channels within the datasets 
%                                              above.
%
% Authors: Arnaud Delorme, CERCO, 2006

% Copyright (C) Arnaud Delorme, CERCO, arno@salk.edu
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

function STUDY = std_changroup(STUDY, ALLEEG, alllocs, interp);

if nargin < 4
    interp = 'off';
end;

% union of all channel structures
% -------------------------------
inputloc = 0;
if nargin >= 3
    if ~isempty(alllocs)
        inputloc = 1;
    end;
end;
if ~inputloc
    alllocs = eeg_mergelocs(ALLEEG.chanlocs);
end;

% create group for each electrode
% -------------------------------
if isstruct(alllocs)
    alllocs = { alllocs.labels };
end;
STUDY.changrp = [];
for indc = 1:length(alllocs)
    STUDY.changrp(indc).name     = alllocs{indc};
    STUDY.changrp(indc).channels = { alllocs{indc} };
    tmp = std_chanlookupnew( STUDY, ALLEEG, STUDY.changrp(indc), interp);
    STUDY.changrp(indc).setinds = tmp.setinds;
    STUDY.changrp(indc).allinds = tmp.allinds;
    STUDY.changrp(indc).centroid = [];
end;

% if strcmpi(interp, 'off')
%     if length(unique( cellfun(@length, { ALLEEG.chanlocs }))) ~= 1
%          STUDY.changrpstatus = 'some channels missing in some datasets';
%     else STUDY.changrpstatus = 'all channels present in all datasets';
%     end;
% else STUDY.changrpstatus = 'all channels present in all datasets - interpolated';
% end;

%STUDY.changrp(indc).name = [ 'full montage' ];
%STUDY.changrp(indc).channels = { alllocs.labels };
%tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(indc));
%STUDY.changrp(indc).chaninds = tmp.chaninds;
return; 
    
% find datasets and channel indices
% ---------------------------------
function changrp = std_chanlookupnew( STUDY, ALLEEG, changrp, interp);

    setinfo       = STUDY.design(STUDY.currentdesign).cell;
    allconditions = STUDY.design(STUDY.currentdesign).variable(1).value;
    allgroups     = STUDY.design(STUDY.currentdesign).variable(2).value;
    nc = max(length(allconditions),1);
    ng = max(length(allgroups),    1);
    changrp.allinds = cell( nc, ng );
    changrp.setinds = cell( nc, ng );
    for index = 1:length(setinfo)
        % get index of independent variables
        % ----------------------------------
        condind = std_indvarmatch( setinfo(index).value{1}, allconditions);
        grpind  = std_indvarmatch( setinfo(index).value{2}, allgroups    );
        if isempty(allconditions), condind = 1; end;
        if isempty(allgroups),     grpind  = 1; end;

        % scan all channel labels
        % -----------------------
        if strcmpi(interp, 'off')
            datind  = setinfo(index).dataset;
            tmpchanlocs = ALLEEG(datind(1)).chanlocs;
            tmplocs = { tmpchanlocs.labels };
            for indc = 1:length(changrp.channels) % usually just one channel
                ind = strmatch( changrp.channels{indc}, tmplocs, 'exact');
                if length(ind) > 1, error([ 'Duplicate channel label ''' tmplocs{ind(1)} ''' for dataset ' int2str(datind) ]); end;
                if ~isempty(ind)
                    changrp.allinds{ condind, grpind } = [ changrp.allinds{ condind, grpind } ind ];
                    changrp.setinds{ condind, grpind } = [ changrp.setinds{ condind, grpind } index ];
                end;
            end;
        else % interpolation is "on", all channels for all datasets
            alllocs = { STUDY.changrp.name };
            ind = strmatch( changrp.name, alllocs, 'exact');
            changrp.allinds{ condind, grpind } = [ changrp.allinds{ condind, grpind } ind   ];
            changrp.setinds{ condind, grpind } = [ changrp.setinds{ condind, grpind } index ];
        end;
    end;
