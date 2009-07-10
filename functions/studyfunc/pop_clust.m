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
%   'clus_num'  - [integer] the number of desired clusters (must be > 1) {default: 20}
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

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

% $Log: not supported by cvs2svn $
% Revision 1.37  2009/06/24 22:10:30  arno
% fixed history
%
% Revision 1.36  2009/04/17 15:01:20  julie
% Remove old clusters when called from commandline
%
% Revision 1.35  2007/12/09 01:06:59  arno
% docuemntation and auto-disable checkbox
%
% Revision 1.34  2007/10/25 21:33:14  nima
% _
%
% Revision 1.33  2007/09/11 10:33:38  arno
% add new free algorithm for kmeans
%
% Revision 1.32  2007/08/07 18:29:22  arno
% fix help from guy (bug 283)
%
% Revision 1.31  2007/08/07 18:26:58  arno
% added help message
%
% Revision 1.30  2007/08/06 18:25:51  scott
% further help msg edits for accuracy  ('kmean' -> 'kmeans', etc).
%
% Revision 1.29  2007/08/06 18:11:48  arno
% update message
%
% Revision 1.28  2007/08/06 18:00:26  scott
% the default value '3' is entered in line 176, I believe...
%
% Revision 1.27  2007/08/06 17:35:41  arno
% update help message according to bug 282
%
% Revision 1.26  2007/05/01 21:05:14  arno
% better help message
%
% Revision 1.25  2006/03/21 15:45:36  arno
% new .sets format; embeding component creation
%
% Revision 1.23  2006/03/12 03:20:35  arno
% study saving
%
% Revision 1.22  2006/03/12 02:56:42  arno
% gui aspect ratio
%
% Revision 1.21  2006/03/11 06:05:50  arno
% header
%
% Revision 1.20  2006/03/11 00:27:40  arno
% header
%
% Revision 1.19  2006/03/10 22:53:17  arno
% typo
%
% Revision 1.18  2006/03/10 18:40:04  arno
% same
%
% Revision 1.17  2006/03/10 18:39:13  arno
% indices removal
%
% Revision 1.16  2006/03/08 20:07:36  arno
% rename
%
% Revision 1.15  2006/03/03 00:43:32  scott
% editing msgs -sm
%
% Revision 1.14  2006/02/22 22:44:39  arno
% now edit all clusters by default
%
% Revision 1.13  2006/02/22 22:40:45  arno
% smarter warning message
%
% Revision 1.12  2006/02/22 22:37:41  arno
% allowing subclustering
%
% Revision 1.11  2006/02/22 21:34:13  arno
% add cluster description to GUI
%
% Revision 1.10  2006/02/14 00:13:32  arno
% adding log
%

function [STUDY, ALLEEG, command] = pop_clust(STUDY, ALLEEG, varargin)

command = '';
if nargin < 2
    help pop_clust;
    return;
end;

if isempty(STUDY.etc)
    error('No pre-clustering information, pre-cluster first!');
end;
if ~isfield(STUDY.etc, 'preclust')
    error('No pre-clustering information, pre-cluster first!');
end;
if isempty(STUDY.etc.preclust)
    error('No pre-clustering information, pre-cluster first!');
end;

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
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) & ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [ rmindex index ];
            end;
        end;        
    end;
    
    if length(STUDY.cluster) > 2 & ~isempty(rmindex)
        resp = questdlg2('Clustering again will delete the last clustering results', 'Warning', 'Cancel', 'Ok', 'Ok');
        if strcmpi(resp, 'cancel'), return; end;
    end;
    
	alg_options = {'Kmeans (stat. toolbox)' 'Neural Network (stat. toolbox)' 'Kmeanscluster (no toolbox)' }; %'Hierarchical tree' 
	set_outliers = ['set(findobj(''parent'', gcbf, ''tag'', ''outliers_std''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));'...
                            'set(findobj(''parent'', gcbf, ''tag'', ''std_txt''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));']; 
	algoptions = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''kmeans''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    if ~exist('kmeans'), valalg = 3; else valalg = 1; end;
                  
    strclust = '';
    if STUDY.etc.preclust.clustlevel > length(STUDY.cluster)
        STUDY.etc.preclust.clustlevel = 1;
    end;
    if STUDY.etc.preclust.clustlevel == 1
        strclust = [ 'Performing clustering on cluster ''' STUDY.cluster(STUDY.etc.preclust.clustlevel).name '''' ];
    else
        strclust = [ 'Performing sub-clustering on cluster ''' STUDY.cluster(STUDY.etc.preclust.clustlevel).name '''' ];
    end;
    
	clust_param = inputgui( { [1] [1] [3 1] [3 1] [ 0.2 2.73 1 ] [1] [0.29 2 2.5 0.8] }, ...
	{ {'style' 'text'       'string' strclust 'fontweight' 'bold'  } {} ...
      {'style' 'text'       'string' 'Clustering algorithm:' } ...
      {'style' 'popupmenu'  'string' alg_options  'value' valalg 'tag' 'clust_algorithm'  'Callback' algoptions } ...
      {'style' 'text'       'string' 'Number of clusters to compute:' } ...
      {'style' 'edit'       'string' '2' 'tag' 'clust_num' }...
      {'style' 'checkbox'   'string' '' 'tag' 'outliers_on' 'value' 0 'Callback' set_outliers 'userdata' 'kmeans' 'enable' 'on' } ...
      {'style' 'text'       'string' 'Separate outliers %gt; [N] std.dev. from any cluster center' 'userdata' 'kmeans' 'enable' 'on' } ...
      {'style' 'edit'       'string' '3' 'tag' 'outliers_std' 'enable' 'off' } {}...
      {'style' 'checkbox'   'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} ...
      {'style' 'text'       'string' 'Save STUDY set to disk'} ...
      {'style' 'edit'       'string' fullfile(STUDY.filepath, STUDY.filename) 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
      {'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} },...
                            'pophelp(''pop_clust'')', 'Set clustering algorithm -- pop_clust()' , [] , 'normal');
	
	if ~isempty(clust_param)
        
        % removing previous cluster information
        % -------------------------------------
        if ~isempty(rmindex)
            fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
            STUDY.cluster(rmindex)          = [];
            STUDY.cluster(clustlevel).child = [];
            if clustlevel == 1 & length(STUDY.cluster) > 1
                STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
            end;
        end;
        
        clus_alg = alg_options{clust_param{1}};
        clus_num = str2num(clust_param{2});
        outliers_on = clust_param{3};
        stdval = clust_param{4};
        save_on = clust_param{5};
                
        outliers = [];
        try
            clustdata = STUDY.etc.preclust.preclustdata;
        catch
            error('Error accesing preclustering data. Perform pre-clustering.');
        end;
        command = '[STUDY] = pop_clust(STUDY, ALLEEG,';
        
        if ~isempty(findstr(clus_alg, 'Kmeanscluster')), clus_alg = 'kmeanscluster'; end;
        if ~isempty(findstr(clus_alg, 'Kmeans ')), clus_alg = 'kmeans'; end;
        if ~isempty(findstr(clus_alg, 'Neural ')), clus_alg = 'neural network'; end;
        
        disp('Clustering ...');
        switch clus_alg
            case { 'kmeans' 'kmeanscluster' }
                command = sprintf('%s %s%s%s %d %s', command, '''algorithm'',''', clus_alg, ''',''clus_num'', ', clus_num, ',');
                if outliers_on
                    command = sprintf('%s %s %s %s', command, '''outliers'', ', stdval, ',');
                    [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,str2num(stdval),5,lower(clus_alg));
                    [STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'robust_kmeans', clus_num});
                else
                    if strcmpi(clus_alg, 'kmeans')
                        [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'replicates',10,'emptyaction','drop');
                    else
                        %[IDX,C,sumd,D] = kmeanscluster(clustdata,clus_num);
                        [C,IDX,sumd] =kmeans_st(real(clustdata),clus_num,150);
                    end;
                    [STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'Kmeans', clus_num});
                end    
         case 'Hierarchical tree'
             %[IDX,C] = hierarchical_tree(clustdata,clus_num);
             %[STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'Neural Network', clus_num});
         case 'neural network'
             [IDX,C] = neural_net(clustdata,clus_num);
             [STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'Neural Network', clus_num});
             command = sprintf('%s %s %d %s', command, '''algorithm'', ''Neural Network'',''clus_num'', ', clus_num, ',');
        end
        disp('Done.');
        
        % If save updated STUDY to disk
        if save_on
            command = sprintf('%s %s', command, '''save'', ''on'',');
            if ~isempty(clust_param{6})
                [filepath filename ext] = fileparts(clust_param{6});
                command = sprintf('%s%s%s%s%s%s', command, '''filename'', ''', [filename ext], ', ''filepath'', ''', filepath, ''');' );
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);
                STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [filename ext], 'filepath', filepath);
              else
                command(end:end+1) = ');';
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, command); 
                if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
                    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
                else
                    STUDY = pop_savestudy(STUDY, ALLEEG);
                end
           end
       else
           command(end:end+1) = ');';
           STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);            
       end
           
       [STUDY com] = pop_clustedit(STUDY, ALLEEG); 
       command = [ command com];
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
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) & ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [ rmindex index ];
            end;
        end;        
    end;
    if ~isempty(rmindex)
        fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
        STUDY.cluster(rmindex)          = [];
        STUDY.cluster(clustlevel).child = [];
        if clustlevel == 1 & length(STUDY.cluster) > 1
            STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
        end;
    end;

    %default values
    algorithm = 'kmeans';
    clus_num = 20;
    save = 'off';
    filename = STUDY.filename;
    filepath = STUDY.filepath;
    outliers = Inf; % default std is Inf - no outliers
    
    if mod(length(varargin),2) ~= 0
        error('pop_clust(): input variables must be specified in pairs: keywords, values');
    end
    
    for k = 1:2:length(varargin)
        switch(varargin{k})
            case 'algorithm'
                algorithm = varargin{k+1};
            case 'clus_num'
                clus_num = varargin{k+1};
            case 'outliers'
                outliers =  varargin{k+1};                
            case 'save'
                save = varargin{k+1};
            case 'filename' 
                filename = varargin{k+1};
            case 'filepath'
                filepath = varargin{k+1};
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
                end;
                [STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'Kmeans', clus_num});
            else
                [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,outliers,5, algorithm);
                [STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'robust_kmeans', clus_num});
            end
        case 'Neural Network'
            [IDX,C] = neural_net(clustdata,clus_num);
            [STUDY, clusters] = std_createclust2(STUDY,IDX,C,  {'Neural Network', clus_num});
        otherwise
            disp('pop_clust: unknown algorithm return');
            return
    end
            % If save updated STUDY to disk
    if strcmpi(save,'on')
        if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
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
function [STUDY, clusters] = std_createclust2(STUDY,IDX,C, algorithm)

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
