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

function STUDY = std_mergeclust(STUDY, ALLEEG, mrg_cls, varargin)

if isempty(varargin) | strcmpi(varargin,'')
    name = 'Cls';
else
    name = varargin{1}; 
end        

% Cannot merge clusters if any of the clusters is a 'Notclust' or 'Outlier'
% cluster, or if has children
comps = [];
sets = [];
for k = 1:length(mrg_cls)
    if strncmpi('Notclust',STUDY.cluster(mrg_cls(k)).name,8) | strncmpi('Outliers',STUDY.cluster(mrg_cls(k)).name,8) | ...
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
diffsets = unique(sets(1,:));
for k = 1:length(diffsets)
    ci = find(sets(1,:) == diffsets(k)); % find the compnents belonging to each set
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


