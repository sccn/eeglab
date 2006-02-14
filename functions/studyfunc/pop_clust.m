% pop_clust() - select and run a clustering algorithm on components from an EEGLAB STUDY 
%               set of EEG datasets. Clustering data should be prepared beforehand using 
%               pop_preclust() and/or eeg_preclust(). The number of clusters must be
%               specified in advance. If called in GUI mode, the pop_clustedit() window
%               appears when the clustering is complete to display clustering results
%               and allow the user to review and edit them.
% Usage: 
%               >> [ STUDY] = pop_clust( STUDY, ALLEEG); % pop up a graphic interface
%               >> [ STUDY] = pop_clust( STUDY, ALLEEG, 'key1', 'val1', ...); % no pop-up
% Inputs:
%   STUDY       - an EEGLAB STUDY set containing some or all of the EEG sets in ALLEEG. 
%   ALLEEG      - a vector of loaded EEG dataset structures of all sets in the STUDY set.
%
% Optional Inputs:
%   'algorithm' - ['kmeans' | 'Neural Network'] algorithm to be used for clustering.
%                 the 'Neural Network' option requires the Matlab Neural Net toolbox 
%                    {default: 'kmeans'} 
%   'clust_num' - [num] the number of desired clusters (must be > 1)
%                    {default: 20}
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
%  "Identify outliers'' - [check box] check to detect outliers. 
%  "Save STUDY"         - [check box] check to save the updated STUDY after clustering 
%                         is performed. If no file entered, overwrites the current STUDY. 
%
%  See also  pop_clustedit(), pop_preclust(), eeg_preclust(), pop_clust()
%
% Authors:  Hilit Serby, SCCN, INC, UCSD, October 11, 2004

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

function [STUDY, ALLEEG, com] = pop_clust(STUDY, ALLEEG, varargin)

com = '';
if isempty(varargin) %GUI call
	alg_options = {'Kmeans' 'Neural Network' }; %'Hierarchical tree' 
	set_outliers = ['set(findobj(''parent'', gcbf, ''tag'', ''outliers_std''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));'...
                            'set(findobj(''parent'', gcbf, ''tag'', ''std_txt''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));']; 
	algoptions = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''kmeans''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
	
	clust_param = inputgui( { [3 1] [3 1] [ 0.2 2.7 1 ] [1] [0.29 2 2.5 0.8] }, ...
	{ {'style' 'text'       'string' 'Clustering algorithm:' } ...
      {'style' 'popupmenu'  'string' alg_options  'value' 1 'tag' 'clust_algorithm'  'Callback' algoptions } ...
      {'style' 'text'       'string' 'Number of clusters to compute:' } ...
      {'style' 'edit'       'string' '2' 'tag' 'clust_num' }...
      {'style' 'checkbox'   'string' '' 'tag' 'outliers_on' 'value' 0 'Callback' set_outliers 'userdata' 'kmeans' 'enable' 'on' } ...
      {'style' 'text'       'string' 'Identify outliers more than X std from cluster centers' 'userdata' 'kmeans' 'enable' 'on' } ...
      {'style' 'edit'       'string' '3' 'tag' 'outliers_std' 'enable' 'off' } {}...
      {'style' 'checkbox'   'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} ...
      {'style' 'text'       'string' 'Save STUDY set to disk'} ...
      {'style' 'edit'       'string' fullfile(STUDY.filepath, STUDY.filename) 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
      {'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} },...
                            'pophelp(''eeg_clust'')', 'Set clustering algorithm -- pop_clust()' , [] , 'normal');
	
	if ~isempty(clust_param)
        com = '% no history yet for clustering';
        
        clus_alg = alg_options{clust_param{1}};
        clus_num = str2num(clust_param{2});
        outliers_on = clust_param{3};
        stdval = clust_param{4};
        save_on = clust_param{5};
        
        outliers = [];
        try
            clustdata = STUDY.etc.preclust.preclustdata;
        catch
            error('Error accesing preclustering data. Study must go through pre-clustering first.');
        end;
        command = '[STUDY] = pop_clust(STUDY, ALLEEG,';
        
        disp('Performing clustering...');
        switch clus_alg
            case 'Kmeans'
                command = sprintf('%s %s %d %s', command, '''algorithm'', ''kmeans'',''clus_num'', ', clus_num, ',');
                if outliers_on
                    [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,str2num(stdval),5);
                    [STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'robust_kmeans', clus_num});
                    command = sprintf('%s %s %s %s', command, '''outliers'', ', stdval, ',');
                else
                    [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'replicates',10,'emptyaction','drop');
                    [STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'Kmeans', clus_num});
                end    
         case 'Hierarchical tree'
             %[IDX,C] = hierarchical_tree(clustdata,clus_num);
             %[STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'Neural Network', clus_num});
         case 'Neural Network'
             [IDX,C] = neural_net(clustdata,clus_num);
             [STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'Neural Network', clus_num});
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
                STUDY = pop_savestudy(STUDY, 'filename', [filename ext], 'filepath', filepath);
              else
                command(end:end+1) = ');';
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, command); 
                if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
                    STUDY = pop_savestudy(STUDY, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
                else
                    STUDY = pop_savestudy(STUDY);
                end
           end
       else
           command(end:end+1) = ');';
           STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);            
       end
           
        [STUDY] = pop_clustedit(STUDY, ALLEEG, clusters);           
	end
    
else %command line call
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
    switch algorithm
        case 'kmeans'                
            if outliers == Inf
                [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'replicates',10,'emptyaction','drop');
                [STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'Kmeans', clus_num});
            else
                [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,outliers,5);
                [STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'robust_kmeans', clus_num});
            end
        case 'Neural Network'
            [IDX,C] = neural_net(clustdata,clus_num);
            [STUDY, clusters] = create_cluster(STUDY,IDX,C,  {'Neural Network', clus_num});
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
