% std_createclust()  - dreate a new empty cluster.  After creation, components 
%                      may be (re)assigned to it using std_movecomp().
% Usage:
%                    >> [STUDY] = std_createclust(STUDY, ALLEEG, 'key', val);
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - vector of EEG datasets included in the STUDY, typically created 
%              using load_ALLEEG().
%
% Optional inputs:
%   'name'     - ['string'] name of the new cluster {default: 'Cls #', where 
%                '#' is the next available cluster number}
%   'clusterind' - [integer] cluster for each of the component. Ex: 61 components 
%                and 2 clusters: 'clusterind' will be a 61x1 vector of 1 and
%                2 (and 0=outlisers)
%   'centroid' - centroid for clusters. If 2 clusters, size will be 2 x
%                width of the preclustering matrix. This is a deprecated 
%                functionality.
%   'algorithm' - [cell] algorithm parameters used to obtain the clusters
%   'parentcluster' - ['on'|'off'] use the parent cluster (cluster 1) to
%                 perform clustering (this cluster contains all the selected
%                 components by default). Otherwise, the cluster defined in
%                 STUDY.etc.preclust.clustlevel is used as parent.
%   
% Outputs:
%   STUDY    - the input STUDY set structure modified with the new cluster.
%
% Example: 
%   >> [STUDY] = std_createclust(STUDY, ALLEEG, 'name', 'eye_movements', ...
%                   'clusterind', [0 1 0 1 0 1], 'parentcluster', 'on');
%   % Create a new cluster named 'eye_movements' with components 2, 4, and
%   % of 6 the default parent cluster defined in 
%
%  See also  pop_clustedit(), std_movecomp()
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

function [STUDY] = std_createclust(STUDY, ALLEEG, varargin)

if nargin< 2
    help std_createclust;
    return;
end

% decoding options for backward compatibility
% -------------------------------------------
options = {};
if length(varargin) > 0 && ~ischar(varargin{1})
    % STUDY, IDX, algorithm, parentClusterNumber
    if isnumeric(ALLEEG)
        options = { options{:} 'clusterind' ALLEEG };
        if nargin > 3, options = { options{:} 'centroid'      varargin{1} }; end
        if nargin > 4, options = { options{:} 'algorithm'     varargin{2} }; end
        ALLEEG = [];
    end
elseif length(varargin) < 2
    options = { options{:} 'name' varargin{1} };
else
    options =  varargin;
end
opt = finputcheck(options, { 'name'             'string'   []  'Cls';
                             'clusterind'       'integer'  []  length(STUDY.cluster)+1;
                             'parentcluster'    'string'   { 'on','off' }  'off';
                             'algorithm'        'cell'     []  {};
                             'ignore0'          'string'   { 'on','off' }  'off';
                             'centroid'         'real'     []  [] }, 'std_createclust');
if ischar(opt), error(opt); end

% opt.clusterind - index of cluster for each component. Ex: 63 components and 2
% clusters: opt.clusterind will be a 61x1 vector of 1 and 2 (and 0=outlisers)
% C - centroid for clusters. If 2 clusters, size will be 2 x
%     width of the preclustering matrix

if strcmpi(opt.parentcluster, 'on')
    firstind = 1;
    cls      = 1;
    sameica  = std_findsameica(ALLEEG);
    sets     = [];
    comps    = [];
    STUDY.cluster = [];
    for index = 1:length(sameica)
        newcomps = STUDY.datasetinfo(sameica{index}(1)).comps;
        if isempty(newcomps), newcomps = [1:size(ALLEEG(sameica{index}(1)).icaweights,1)]; end
        comps = [ comps newcomps ];
        sets(length(sameica{index}):-1:1,end+1:end+length(newcomps)) = repmat( sameica{index}', [1 length(newcomps) ] );
    end
    sets(find(sets == 0))   = NaN;
    STUDY.cluster(1).name   = 'Parentcluster 1';
    STUDY.cluster(1).sets   = sets;
    STUDY.cluster(1).comps  = comps;
    STUDY.cluster(1).parent = {};
    STUDY.cluster(1).child  = {};
    STUDY.cluster.preclust.preclustparams = [];    
    STUDY.cluster.preclust.preclustdata   = []; 
else
    % Find the next available cluster index
    % -------------------------------------
    cls = min(max(opt.clusterind), length(unique(opt.clusterind)));
    nc  = 0; % index of last cluster
    for k =  1:length(STUDY.cluster)
        ti = strfind(STUDY.cluster(k).name, ' ');
        tmp = STUDY.cluster(k).name(ti(end) + 1:end);
        nc = max(nc,str2num(tmp));
        % check if there is a cluster of Notclust components
        if isfield(STUDY.etc, 'preclust') && isfield(STUDY.etc.preclust, 'preclustparams')
            if ~isempty(STUDY.cluster(k).parent)
                %strcmp(STUDY.cluster(k).parent,STUDY.cluster(STUDY.etc.preclust.clustlevel).name) 
                STUDY.cluster(k).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
            end
        end
    end
    len = length(STUDY.cluster);

    if ~isempty(find(opt.clusterind==0)) && strcmpi(opt.ignore0, 'off') %outliers exist
        firstind = 0;
        nc  = nc + 1;
        len = len + 1;
    else
        firstind = 1;
    end

    % create clustlevel if it does not exist
    % --------------------------------------
    if ~isfield(STUDY.etc, 'preclust')
        STUDY.etc.preclust.clustlevel = 1;
        STUDY.etc.preclust.preclustdata = [];
    elseif ~isfield(STUDY.etc.preclust, 'clustlevel')
        STUDY.etc.preclust.clustlevel = 1;
        STUDY.etc.preclust.preclustdata = [];
    end
    
    % create all clusters
    % -------------------
    for k = firstind:cls 

        % cluster name
        % ------------
        if k == 0
             STUDY.cluster(len).name   = [ 'outlier ' num2str(k+nc)];
        else STUDY.cluster(k+len).name = [ opt.name ' ' num2str(k+nc)];
        end

        % find indices
        % ------------
        tmp = find(opt.clusterind==k); % opt.clust.erind contains the cluster index for each component
        STUDY.cluster(k+len).sets  = STUDY.cluster(STUDY.etc.preclust.clustlevel).sets(:,tmp);
        STUDY.cluster(k+len).comps = STUDY.cluster(STUDY.etc.preclust.clustlevel).comps(tmp);
        STUDY.cluster(k+len).algorithm = opt.algorithm;
        STUDY.cluster(k+len).parent{end+1} = STUDY.cluster(STUDY.etc.preclust.clustlevel).name;
        STUDY.cluster(k+len).child = [];
        if ~isempty(STUDY.etc.preclust.preclustdata) && all(tmp <= size(STUDY.etc.preclust.preclustdata,1))
             STUDY.cluster(k+len).preclust.preclustdata   = STUDY.etc.preclust.preclustdata(tmp,:);
             STUDY.cluster(k+len).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
        else STUDY.cluster(k+len).preclust.preclustdata   = [];
        end

        %update parents clusters with cluster child indices
        % -------------------------------------------------
        STUDY.cluster(STUDY.etc.preclust.clustlevel).child{end+1} = STUDY.cluster(k+nc).name;
    end
end


% Find out the highst cluster id number (in cluster name), to find
% next available cluster index


% % find max cluster ID
% 
% max_id = 0;
% if ~isfield(STUDY, 'cluster'), STUDY.cluster = []; end
% for k = 1:length(STUDY.cluster)
%     ti = strfind(STUDY.cluster(k).name, ' ');
%     clus_id = STUDY.cluster(k).name(ti(end) + 1:end);
%     max_id = max(max_id, str2num(clus_id));
% end
% max_id = max_id + 1;
% opt.name = sprintf('%s %d', opt.name, max_id);
% clustind = length(STUDY.cluster)+1;
% % Initialize the new cluster fields.
% if length(STUDY.cluster) > 0
%     STUDY.cluster(clustind).parent{1} = STUDY.cluster(1).name;
%     if ~iscell(STUDY.cluster(1).child)
%          STUDY.cluster(1).child = { opt.name };
%     else STUDY.cluster(1).child = { STUDY.cluster(1).child{:} opt.name };
%     end
% else  
%     STUDY.cluster(clustind).parent{1} = 'manual'; % update parent cluster if exists.
% end
% STUDY.cluster(clustind).name = opt.name;
% STUDY.cluster(clustind).child = [];
% STUDY.cluster(clustind).comps = [];
% STUDY.cluster(clustind).sets = [];
% STUDY.cluster(clustind).algorithm = [];
% STUDY.cluster(clustind).centroid = [];
% STUDY.cluster(clustind).preclust.preclustparams = [];
% STUDY.cluster(clustind).preclust.preclustdata = [];
% 
% if (~isempty(opt.datasets) || ~isempty(opt.subjects)) && ~isempty(opt.components)
%     
%     % convert subjects to dataset indices
%     % -----------------------------------
%     if ~isempty(opt.subjects)
%         if length(opt.subjects) ~= length(opt.components)
%             error('If subjects are specified, the length of the cell array must be the same as for the components');
%         end
%         alls = { ALLEEG.subject };
%         for index = 1:length(opt.subjects)
%             tmpinds = strmatch(opt.subjects{index}, alls, 'exact');
%             if isempty(tmpinds)
%                 error('Cannot find subject');
%             end
%             opt.datasets(1:length(tmpinds),index) = tmpinds;
%         end
%         opt.datasets(opt.datasets(:) == 0) = NaN;
%     end
%     
%     % deal with cell array inputs
%     % ---------------------------
%     if iscell(opt.components)
%         newcomps = [];
%         newdats  = [];
%         for ind1 = 1:length(opt.components)
%             for ind2 = 1:length(opt.components{ind1})
%                 if iscell(opt.datasets)
%                      newdats  = [ newdats  opt.datasets{ind1}' ];
%                 else newdats  = [ newdats  opt.datasets(:,ind1) ];
%                 end
%                 newcomps = [ newcomps opt.components{ind1}(ind2) ];
%             end
%         end
%         opt.datasets   = newdats;
%         opt.components = newcomps;
%     end
%     
%     % create .sets, .comps, .setinds, .allinds fields
%     % -----------------------------------------------
%     [tmp setinds allinds] = std_setcomps2cell( STUDY, opt.datasets, opt.components);
%     STUDY.cluster(clustind).setinds = setinds;
%     STUDY.cluster(clustind).allinds = allinds;
%     STUDY.cluster(clustind) = std_cell2setcomps( STUDY, ALLEEG, clustind); 
%     STUDY.cluster(clustind) = std_setcomps2cell( STUDY, clustind);
%     %[ STUDY.cluster(finalinds(ind)) setinds allinds ] =
%         %std_setcomps2cell(STUDY, finalinds(ind));
% end
