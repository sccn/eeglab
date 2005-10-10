% Usage: >> [ STUDY] = pop_clust(ALLEEG, STUDY); %Graphical interface call
%            >> [ STUDY] = pop_clust(ALLEEG, STUDY, 'key1', 'val1', ...); %command line call 
%
% pop_clust() - a function (both GUI and command line modes) to choose and run a
% clustering algorithm on a pre-prepared data. The data is prepared using
% the GUI pop_preclust function or the command line eeg_preclust  function.
% The function allows you to select a clustering algorithm and requires to set 
% the number of clusters in advance. In the GUI mode when the clustering is
% completed the pop_clustedit function is called, to present the clustering results
% and allow the user to edit them.
%
%
% Inputs:
%   ALLEEG     - an EEGLAB data structure, which holds EEG sets (can also be one EEG set). 
%                      ALLEEG must contain all the datasets that were used
%                      in the STUDY. 
%   STUDY      - an eeglab STUDY set, can contain some or all of the EEGs
%                      structures in ALLEEG. 
% Optional Inputs:
%   'algorithm'  - ['kmeans' | 'Neural Network'] what algorithm will be used for clustering.
%                         {default: 'kmeans'} 
%   'clust_num' - [num] the number of desired clusters (must be bigger than one). 
%                         {default: 20}
%   'save'          - ['on' | 'off'] save the updated STUDY. {default: 'off'} 
%   'filename'   - [string] if save option is on, will save the STUDY under this file name. 
%                         {default: current STUDY filename}
%   'filepath'   - [string] if save option is on, will save the STUDY in this directory. 
%                         {default: current STUDY filepath}
% Outputs:
%   STUDY - same as inputs only modified with the clustering algorithm results.
%
% Graphic interface buttons:
%   ''Clustering algorithm'' - [list box] displays an optional clustering algorithms. 
%                         At the moment only the k-means algorithm is implemented.
%   ''Number of clusters to compute'' - [edit box] the number of desired clusters,
%                         minimum is 2.
%   ''Identify outliers'' - [check box] optional algorithm choice to detect outliers. 
%   'Save STUDY'' - [check box] an option to save the STUDY after clustering is performed. 
%                         If  no file is entered will overwrite the current file of STUDY. 
%
%  See also  pop_clustedit, pop_preclust, eeg_preclust, pop_clust         
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

function [STUDY] = pop_clust(ALLEEG, STUDY, varargin)

if isempty(varargin) %GUI call
	alg_options = {'Kmeans' 'Neural Network' }; %'Hierarchical tree' 
	set_outliers = ['set(findobj(''parent'', gcbf, ''tag'', ''outliers_std''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));'...
                            'set(findobj(''parent'', gcbf, ''tag'', ''std_txt''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));']; 
	algoptions = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''kmeans''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
	
	clust_param = inputgui( { [1 1] [1] [2 1] [1] [ 1 2 1 2 ] [1] [0.5 2 2.5 0.8] [1] }, ...
	{ { 'style' 'text' 'string' 'Clustering algorithm:' } {'style' 'listbox' 'string' alg_options  'value' 1 'tag' 'clust_algorithm'  'Callback' algoptions } {}...
	{'style' 'text' 'string' 'Number of clusters to compute:' } { 'style' 'edit' 'string' '2' 'tag' 'clust_num' } {}...
	{ 'style' 'checkbox' 'tag' 'outliers_on' 'string' '' 'value' 0 'Callback' set_outliers 'userdata' 'kmeans' 'enable' 'on' } ...
	{'style' 'text' 'string' 'Identify outliers more than' 'userdata' 'kmeans' 'enable' 'on' } ...
	{ 'style' 'edit' 'string' '3' 'tag' 'outliers_std' 'enable' 'off' } ...
	{'style' 'text' 'string' 'std from cluster centers'  'tag'  'std_txt' 'enable' 'off' } {}...
	{'style' 'checkbox'  'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} {'style' 'text' 'string' 'Save STUDY set to disk'} ...
	{'style' 'edit' 'string' '' 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
	{'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} {} },...
	'pophelp(''eeg_clust'')', 'Set clustering algorithm -- pop_clust()' , [] , 'normal', ...
	[2 1 1 1 1 1 1 1 ]);
	
	if ~isempty(clust_param)
        clus_alg = alg_options{clust_param{1}};
        clus_num = str2num(clust_param{2});
        outliers_on = clust_param{3};
        stdval = clust_param{4};
        save_on = clust_param{5};
        
        outliers = [];
        clustdata = STUDY.etc.preclust.preclustdata;
        command = '[STUDY] = pop_clust(ALLEEG, STUDY,';
        
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
           
        [STUDY] = pop_clustedit(ALLEEG, STUDY,clusters);           
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
        error('pop_clust: input varaibles must come in pairs of keyx and valx');
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

