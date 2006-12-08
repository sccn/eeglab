% std_changroup() - Create channel groups for plotting.
%
% Usage:    
%                >> STUDY = std_changroup(STUDY, ALLEEG);   
% Inputs:
%   ALLEEG     - Top-level EEGLAB vector of loaded EEG structures for the dataset(s) 
%                in the STUDY. ALLEEG for a STUDY set is typically loaded using 
%                pop_loadstudy(), or in creating a new STUDY, using pop_createstudy().  
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user edits,
%                if any. Plotted channel measure means (maps, ERSPs, etc.) are added to 
%                the STUDY structure after they are first plotted to allow quick replotting.  
%
% Authors: Arnaud Delorme, CERCO, 2006

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.3  2006/12/08 19:40:49  arno
% error for duplicate channel label
%
% Revision 1.2  2006/11/15 22:55:46  arno
% last channel problem
%
% Revision 1.1  2006/09/12 18:45:16  arno
% Initial revision
%

function STUDY = std_changroup(STUDY, ALLEEG);

% union of all channel structures
% -------------------------------
alllocs = ALLEEG(STUDY.datasetinfo(1).index).chanlocs;
alllabs = { alllocs.labels };
for index = 2:length(STUDY.datasetinfo)
   tmplocs = ALLEEG(STUDY.datasetinfo(index).index).chanlocs;
   alllocs = eeg_mergechan(alllocs, tmplocs);
end;

% create group for each electrode
% -------------------------------
for indc = 1:length(alllocs)
    STUDY.changrp(indc).name = [ alllocs(indc).labels ];
    STUDY.changrp(indc).channels = { alllocs(indc).labels };
    tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(indc));
    STUDY.changrp(indc).chaninds = tmp.chaninds;
    STUDY.changrp(indc).centroid = [];
end;
%STUDY.changrp(indc).name = [ 'full montage' ];
%STUDY.changrp(indc).channels = { alllocs.labels };
%tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(indc));
%STUDY.changrp(indc).chaninds = tmp.chaninds;

% ---------------
% channel look-up
% ---------------
function changrp = std_chanlookup( STUDY, ALLEEG, changrp);

    changrp.chaninds = [];
    changrp.chaninds = zeros(size(STUDY.setind));
    for ir = 1:size(STUDY.setind,1)
        for ic = 1:size(STUDY.setind,2)
            datind  = STUDY.setind(ir,ic);
            tmplocs = { ALLEEG(STUDY.datasetinfo(datind).index).chanlocs.labels };

            for indc = 1:length(changrp.channels)
                ind = strmatch( changrp.channels(indc), tmplocs, 'exact');
                if length(ind) > 1, error([ 'Duplicate channel label ''' tmplocs{ind(1)} ''' for dataset ' int2str(datind) ]); end;
                if ~isempty(ind)
                    changrp.chaninds(ir,ic) = ind;
                end;
            end;
        end;    
    end;
