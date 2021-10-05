% std_mergeclust()  - Commandline function, to merge several clusters. 
% Usage:    
%                   >> [STUDY] = std_mergeclust(STUDY, ALLEEG, mrg_cls, name);   
% Inputs:
%   STUDY         - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG        - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                       ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   mrg_cls        - clusters indexes to merge.  
% Optional inputs:
%   name           - [string] a mnemonic cluster name for the merged cluster. 
%                       {default: 'Cls #', where '#' is the next available cluster number}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with the merged cluster. 
%
%   Example:
%                         >> mrg_cls = [3 7 9]; name = 'eyes';  
%                         >> [STUDY] = std_mergecluster(STUDY,ALLEEG, mrg_cls, name);
%                    Merge clusters 3, 7 and 9 to a new cluster named 'eyes'. 
%
%  See also  pop_clustedit         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, July, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, July 11, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_mergeclust(STUDY, ALLEEG, mrg_cls, varargin)

if isempty(varargin) || strcmpi(varargin,'')
    name = 'Cls';
else
    name = varargin{1}; 
end        

% Cannot merge clusters if any of the clusters is a 'Notclust' or 'Outlier'
% cluster, or if has children
comps = [];
sets = [];
for k = 1:length(mrg_cls)
    if strncmpi('Notclust',STUDY.cluster(mrg_cls(k)).name,8) || strncmpi('Outliers',STUDY.cluster(mrg_cls(k)).name,8) || ...
        ~isempty(STUDY.cluster(mrg_cls(k)).child)
        warndlg2([ 'std_mergeclust: cannot merge clusters if one of the clusters '...
        'is a ''Notclust'' or ''Outliers'' cluster, or if it has children clusters.']);
    end
    parent{k} = STUDY.cluster(mrg_cls(k)).name;
    comps = [comps STUDY.cluster(mrg_cls(k)).comps];
    sets = [sets STUDY.cluster(mrg_cls(k)).sets];
end

%sort by sets
[tmp,sind] = sort(sets(1,:));
sets = sets(:,sind);
comps = comps(sind);
% sort component indexes within a set
diffsets = unique_bc(sets(1,:));
for k = 1:length(diffsets)
    ci = find(sets(1,:) == diffsets(k)); % find the components belonging to each set
    [tmp,cind] = sort(comps(ci));
    comps(ci) = comps(ci(cind));
end

% Create a new empty cluster
[STUDY] = std_createclust(STUDY, ALLEEG, name);  
% Update merge cluster with parent clusters
STUDY.cluster(end).parent = parent;
STUDY.cluster(end).sets = sets; % Update merge cluster with merged component sets
STUDY.cluster(end).comps = comps; % Update merge cluster with the merged components
for k = 1:length(mrg_cls) % update the merge cluster as a child for the parent clusters
    STUDY.cluster(mrg_cls(k)).child{end + 1} = STUDY.cluster(end).name;
end
STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign);


