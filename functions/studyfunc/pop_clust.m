% pop_clust() - select and run a clustering algorithm on components from an EEGLAB STUDY 
%               structure of EEG datasets. Clustering data should be prepared beforehand using 
%               pop_preclust() and/or std_preclust(). The number of clusters must be
%               specified in advance. If called in gui mode, the pop_clustedit() window
%               appears when the clustering is complete to display clustering results
%               and allow the user to review and edit them.
% Usage: 
%               >> STUDY = pop_clust( STUDY, ALLEEG); % pop up a graphic interface
%               >> STUDY = pop_clust( STUDY, ALLEEG, 'key1', 'val1', ...); % no pop-up
% Inputs:
%   STUDY       - an EEGLAB STUDY set containing some or all of the EEG sets in ALLEEG. 
%   ALLEEG      - a vector of loaded EEG dataset structures of all sets in the STUDY set.
%
% Optional Inputs:
%   'algorithm' - ['kmeans'|'kmeanscluster'|'Neural Network'] algorithm to be used for 
%                 clustering. The 'kmeans' options requires the statistical toolbox. The
%                 'kmeanscluster' option is included in EEGLAB. The 'Neural Network' 
%                  option requires the Matlab Neural Net toolbox {default: 'kmeans'} 
%   'clus_num'  - [integer] the number of desired clusters (must be > 1)
%                 {default: 20}. Not necessary when using Affinity Propagation algorithm 
%   'maxiter'   - maximum number of iterations when using Affinity Propagation algorithm 
%   'outliers'  - [integer] identify outliers further than the given number of standard
%                 deviations from any cluster centroid. Inf --> identify no such outliers.
%                 {default: Inf from the command line; 3 for 'kmeans' from the pop window}
%   'save'      - ['on' | 'off'] save the updated STUDY to disk {default: 'off'} 
%   'filename'  - [string] if save option is 'on', save the STUDY under this file name
%                    {default: current STUDY filename}
%   'filepath'  - [string] if save option is 'on', will save the STUDY in this directory 
%                    {default: current STUDY filepath}
% Outputs:
%   STUDY       - as input, but modified adding the clustering results.
%
% Graphic interface buttons:
%  "Clustering algorithm" - [list box] display/choose among the available clustering 
%                           algorithms. 
%  "Number of clusters to compute" - [edit box] the number of desired clusters (>2)
%  "Identify outliers"  - [check box] check to detect outliers. 
%  "Save STUDY"         - [check box] check to save the updated STUDY after clustering 
%                         is performed. If no file entered, overwrites the current STUDY. 
%
%  See also  pop_clustedit(), pop_preclust(), std_preclust(), pop_clust()
%
% Authors:  Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

function [STUDY, ALLEEG, command] = pop_clust(STUDY, ALLEEG, varargin)

command = '';
if nargin < 2
    help pop_clust;
    return;
end

if isempty(STUDY.etc)
    error('No pre-clustering information, pre-cluster first!');
end
if ~isfield(STUDY.etc, 'preclust')
    error('No pre-clustering information, pre-cluster first!');
end
if isempty(STUDY.etc.preclust)
    error('No pre-clustering information, pre-cluster first!');
end

% Check that path to the stats toolbox comes first (conflict with Fieldtrip)
flagstats = strcmp(regexp(which('kmeans'), '(?<=[\\/]toolbox[\\/])[^\\/]+', 'match', 'once'),'stats');
if ~flagstats
    kmeansPath = fileparts(which('kmeans'));
    rmpath(kmeansPath);
    addpath(kmeansPath);
end
        
if isempty(varargin) %GUI call
    
    % remove clusters below clustering level (done also after GUI)
    % --------------------------------------
    rmindex    = [];
    clustlevel = STUDY.etc.preclust.clustlevel;
    nameclustbase = STUDY.cluster(clustlevel).name;
    if clustlevel == 1
        rmindex = [2:length(STUDY.cluster)];
    else
        for index = 2:length(STUDY.cluster)
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) && ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [ rmindex index ];
            end
        end;        
    end
    
    if length(STUDY.cluster) > 2 && ~isempty(rmindex)
        resp = questdlg2('Clustering again will delete the last clustering results', 'Warning', 'Cancel', 'Ok', 'Ok');
        if strcmpi(resp, 'cancel'), return; end
    end
    
	alg_options = {'Kmeans (stat. toolbox)' 'Neural Network (stat. toolbox)' 'Kmeanscluster (no toolbox)' 'Affinity Propagation' }; %'Hierarchical tree' 
	set_outliers = ['set(findobj(''parent'', gcbf, ''tag'', ''outliers_std''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));'...
                            'set(findobj(''parent'', gcbf, ''tag'', ''std_txt''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));']; 
                        
	algoptions = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''kmeans''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''), ''visible'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ...
                  'if  get(findobj(''parent'', gcbf, ''tag'', ''clust_algorithm''),''value'') == 4 ;'...
                  'set(findobj(''parent'', gcbf, ''userdata'', ''clust_num''), ''enable'', ''off'', ''visible'', ''off''); else; set(findobj(''parent'', gcbf, ''userdata'', ''clust_num''), ''enable'', ''on'', ''visible'', ''on''); end;'];
	
    saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    if ~exist('kmeans'), valalg = 3; else valalg = 1; end
                  
    strclust = '';
    if STUDY.etc.preclust.clustlevel > length(STUDY.cluster)
        STUDY.etc.preclust.clustlevel = 1;
    end
    if STUDY.etc.preclust.clustlevel == 1
        strclust = [ 'Performing clustering on cluster ''' STUDY.cluster(STUDY.etc.preclust.clustlevel).name '''' ];
    else
        strclust = [ 'Performing sub-clustering on cluster ''' STUDY.cluster(STUDY.etc.preclust.clustlevel).name '''' ];
    end
    
    numClust = ceil(mean(cellfun(@length, { STUDY.datasetinfo.comps })));
    if numClust > 2, numClustStr = num2str(numClust);
    else             numClustStr = '10';
    end
    
	clust_param = inputgui( { [1] [1] [1 1] [1 0.5 0.5 ] [ 1 0.5 0.5 ] }, ...
	{ {'style' 'text'       'string' strclust 'fontweight' 'bold'  } {} ...
      {'style' 'text'       'string' 'Clustering algorithm:' } ...
      {'style' 'popupmenu'  'string' alg_options  'value' valalg 'tag' 'clust_algorithm'  'Callback' algoptions } ...
      {'style' 'text'       'string' 'Number of clusters to compute:' 'userdata' 'clust_num' } ...
      {'style' 'edit'       'string' numClustStr 'tag' 'clust_num' 'userdata' 'clust_num' } {} ...
      {'style' 'checkbox'   'string' 'Separate outliers (enter std.)' 'tag' 'outliers_on' 'value' 0 'Callback' set_outliers 'userdata' 'kmeans' 'enable' 'on' } ...
      {'style' 'edit'       'string' '3' 'tag' 'outliers_std' 'userdata' 'kmeans' 'enable' 'off' } {} },...
                            'pophelp(''pop_clust'')', 'Set clustering algorithm -- pop_clust()' , [] , 'normal', [ 1 .5 1 1 1]);
	
	if ~isempty(clust_param)
        
        % removing previous cluster information
        % -------------------------------------
        if ~isempty(rmindex)
            fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
            STUDY.cluster(rmindex)          = [];
            STUDY.cluster(clustlevel).child = [];
            if clustlevel == 1 && length(STUDY.cluster) > 1
                STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
            end
        end
        
        clus_alg = alg_options{clust_param{1}};
        clus_num = str2num(clust_param{2});
        outliers_on = clust_param{3};
        stdval = clust_param{4};
                
        outliers = [];
        try
            clustdata = STUDY.etc.preclust.preclustdata;
        catch
            error('Error accessing preclustering data. Perform pre-clustering.');
        end
        command = '[STUDY] = pop_clust(STUDY, ALLEEG,';
        
        if ~isempty(findstr(clus_alg, 'Kmeanscluster')), clus_alg = 'kmeanscluster'; end
        if ~isempty(findstr(clus_alg, 'Kmeans ')), clus_alg = 'kmeans'; end
        if ~isempty(findstr(clus_alg, 'Neural ')), clus_alg = 'neural network'; end
        
        % Cleaning cache
        STUDY.cache = [];
        
        disp('Clustering ...');
        
        switch clus_alg
            case { 'kmeans' 'kmeanscluster' }
                command = sprintf('%s %s%s%s %d %s', command, '''algorithm'',''', clus_alg, ''',''clus_num'', ', clus_num, ',');
                if outliers_on
                    command = sprintf('%s %s %s %s', command, '''outliers'', ', stdval, ',');
                    [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,str2num(stdval),5,lower(clus_alg));
                    [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'robust_kmeans', clus_num});
                else
                    if strcmpi(clus_alg, 'kmeans')
                        [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'replicates',10,'emptyaction','drop');
                    else
                        %[IDX,C,sumd,D] = kmeanscluster(clustdata,clus_num);
                        [C,IDX,sumd] =kmeans_st(real(clustdata),clus_num,150);
                    end
                    [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Kmeans', clus_num});
                end    
         case 'Hierarchical tree'
             %[IDX,C] = hierarchical_tree(clustdata,clus_num);
             %[STUDY] = std_createclust(STUDY,IDX,C,  {'Neural Network', clus_num});
         case 'neural network'
             [IDX,C] = neural_net(clustdata,clus_num);
             [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm',  {'Neural Network', clus_num});
             command = sprintf('%s %s %d %s', command, '''algorithm'', ''Neural Network'',''clus_num'', ', clus_num, ',');
             
         case 'Affinity Propagation'
             command = sprintf('%s %s%s%s %d %s', command, '''algorithm'',''Affinity Propagation'',');
             plugin_askinstall('limo_eeg', 'apcluster');
             [IDX,C,sumd] = std_apcluster(clustdata,'maxits',200);
             [STUDY]      = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Affinity Propagation',size(C,1)});
        end
        disp('Done.');
        
        % If save updated STUDY to disk
        save_on = 0; % old option to save STUDY
        if save_on
            command = sprintf('%s %s', command, '''save'', ''on'',');
            if ~isempty(clust_param{6})
                [filepath filename ext] = fileparts(clust_param{6});
                command = sprintf('%s%s%s%s%s%s', command, '''filename'', ''', [filename ext], ', ''filepath'', ''', filepath, ''');' );
                STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [filename ext], 'filepath', filepath);
              else
                command(end:end+1) = ');';
                if (~isempty(STUDY.filename)) && (~isempty(STUDY.filepath))
                    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
                else
                    STUDY = pop_savestudy(STUDY, ALLEEG);
                end
           end
       else
           command(end:end+1) = ');';
        end
           
       % Call menu to plot clusters (use EEGLAB menu which include std_envtopo)
       LASTCOM = command;
       eval([ get(findobj(findobj('tag', 'EEGLAB'), 'Label', 'Edit/plot component clusters'), 'callback') ] );
       %[STUDY com] = pop_clustedit(STUDY, ALLEEG); 
	end
    
else %command line call
    % remove clusters below clustering level (done also after GUI)
    % --------------------------------------
    rmindex    = [];
    clustlevel = STUDY.etc.preclust.clustlevel;
    nameclustbase = STUDY.cluster(clustlevel).name;
    if clustlevel == 1
        rmindex = [2:length(STUDY.cluster)];
    else
        for index = 2:length(STUDY.cluster)
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) && ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [ rmindex index ];
            end
        end;        
    end
    if ~isempty(rmindex)
        fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
        STUDY.cluster(rmindex)          = [];
        STUDY.cluster(clustlevel).child = [];
        if clustlevel == 1 && length(STUDY.cluster) > 1
            STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
        end
    end

    %default values
    algorithm = 'kmeans';
    clus_num  = 20;
    save     = 'off';
    filename = STUDY.filename;
    filepath = STUDY.filepath;
    outliers = Inf; % default std is Inf - no outliers
    maxiter  = 200;
    
    if mod(length(varargin),2) ~= 0
        error('pop_clust(): input variables must be specified in pairs: keywords, values');
    end
    
    for k = 1:2:length(varargin)
        switch(varargin{k})
            case 'algorithm'
                algorithm  = varargin{k+1};
            case 'clus_num'
                clus_num   = varargin{k+1};
            case 'outliers'
                outliers   =  varargin{k+1};                
            case 'save'
                save       = varargin{k+1};
            case 'filename' 
                filename   = varargin{k+1};
            case 'filepath'
                filepath   = varargin{k+1};
            case 'maxiter' 
                maxiter    = varargin{k+1};
        end
    end    
    if clus_num < 2
        clus_num = 2;
    end
    
    clustdata = STUDY.etc.preclust.preclustdata;
    switch lower(algorithm)
        case { 'kmeans' 'kmeanscluster' }
            if outliers == Inf
                if strcmpi(algorithm, 'kmeans')
                    [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'replicates',10,'emptyaction','drop');
                else
                    [IDX,C,sumd,D] = kmeanscluster(clustdata,clus_num);
                end
                [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Kmeans', clus_num});
            else
                [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,outliers,5, algorithm);
                [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'robust_kmeans', clus_num});
            end
        case 'neural network'
            [IDX,C] = neural_net(clustdata,clus_num);
            [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm',  {'Neural Network', clus_num});
        case 'affinity propagation'           
             [IDX,C,sumd] = std_apcluster(clustdata,'maxits',maxiter);
             [STUDY]      = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Affinity Propagation',size(C,1)});
        otherwise
            disp('pop_clust: unknown algorithm return');
            return
    end
            % If save updated STUDY to disk
    if strcmpi(save,'on')
        if (~isempty(STUDY.filename)) && (~isempty(STUDY.filepath))
            STUDY = pop_savestudy(STUDY, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
        else
            STUDY = pop_savestudy(STUDY);
        end
   end       
end
STUDY.saved = 'no';


% IDX - index of cluster for each component. Ex: 63 components and 2
% clusters: IDX will be a 61x1 vector of 1 and 2 (and 0=outlisers)
% C - centroid for clusters. If 2 clusters, size will be 2 x
%     width of the preclustering matrix
function [STUDY] = std_createclust2_old(STUDY,IDX,C, algorithm)

% Find the next available cluster index
% -------------------------------------
clusters = [];
cls = size(C,1); % number of cluster = number of row of centroid matrix
nc  = 0; % index of last cluster 
for k =  1:length(STUDY.cluster)
    ti = strfind(STUDY.cluster(k).name, ' ');
    tmp = STUDY.cluster(k).name(ti(end) + 1:end);
    nc = max(nc,str2num(tmp));
    % check if there is a cluster of Notclust components
    if strcmp(STUDY.cluster(k).parent,STUDY.cluster(STUDY.etc.preclust.clustlevel).name) 
        STUDY.cluster(k).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
        clusters = [clusters k];
    end
end
len = length(STUDY.cluster);

if ~isempty(find(IDX==0)) %outliers exist
    firstind = 0;
    nc  = nc + 1;
    len = len + 1;
else
    firstind = 1;
end

% create all clusters
% -------------------
for k = firstind:cls 
    
    % cluster name
    % ------------
    if k == 0
         STUDY.cluster(len).name   = [ 'outlier ' num2str(k+nc)];
    else STUDY.cluster(k+len).name = [ 'Cls ' num2str(k+nc)];
    end

    % find indices
    % ------------
    tmp = find(IDX==k); % IDX contains the cluster index for each component
    STUDY.cluster(k+len).sets  = STUDY.cluster(STUDY.etc.preclust.clustlevel).sets(:,tmp);
    STUDY.cluster(k+len).comps = STUDY.cluster(STUDY.etc.preclust.clustlevel).comps(tmp);
    STUDY.cluster(k+len).algorithm = algorithm;
    STUDY.cluster(k+len).parent{end+1} = STUDY.cluster(STUDY.etc.preclust.clustlevel).name;
    STUDY.cluster(k+len).child = [];
    STUDY.cluster(k+len).preclust.preclustdata = STUDY.etc.preclust.preclustdata(tmp,:);
    STUDY.cluster(k+len).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
    STUDY.cluster(k+len).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;

    %update parents clusters with cluster child indices
    % -------------------------------------------------
    STUDY.cluster(STUDY.etc.preclust.clustlevel).child{end+1} = STUDY.cluster(k+nc).name;
end

clusters = [ clusters firstind+len:cls+len];%the new created clusters indices.
