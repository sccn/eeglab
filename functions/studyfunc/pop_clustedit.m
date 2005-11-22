% pop_clustedit() - graphic user interface (GUI)-based function with editing and plotting 
%                options for visualizing and manipulating an input STUDY structure. 
%                Only component measures (e.g., dipole locations, scalp maps, spectra, 
%                ERPs, ERSPs, ITCs) that have been computed and saved in the study EEG 
%                datasets can be visualized. These can be computed during pre-clustering 
%                using the GUI-based function pop_preclust() or the equivalent command
%                line functions eeg_createdata() and eeg_preclust(). To use dipole 
%                locations for clustering, they must first be stored in the EEG dataset 
%                structures using dipfit(). Supported cluster editing functions include 
%                new cluster creation, cluster merging, outlier rejection, and cluster 
%                renaming. Components can also be moved from one cluster to another 
%                or to the outlier cluster. 
% Usage:    
%                >> STUDY = pop_clustedit(STUDY, ALLEEG, clusters);   
% Inputs:
%   ALLEEG     - Top-level EEGLAB vector of loaded EEG structures for the dataset(s) 
%                in the STUDY. ALLEEG for a STUDY set is typically loaded using 
%                pop_loadstudy(), or in creating a new STUDY, using pop_createstudy().  
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%
% Optional inputs:
%   clusters   - [integer vector] of cluster numbers. These clusters will be visualized 
%                and manipulated in the pop_clustedit() graphic interface. There are 
%                restrictions on which clusters can be loaded together. The clusters must 
%                either originate from the same clustering (same pre_clustering() and 
%                subsequent pop_clust() execution), or they must all be leaf clusters 
%                (i.e., clusters with no child clusters) {default: all leaf clusters}.
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user edits,
%                if any. Plotted cluster measure means (maps, ERSPs, etc.) are added to 
%                the STUDY structure after they are first plotted to allow quick replotting.  
%
% Graphic interface buttons:
%  "Select cluster to plot" - [list box] Displays available clusters to plot (format: 
%                'cluster name (number of components)'). The presented clusters depend 
%                on the optional input variable 'clusters'. Selecting (clicking on) a 
%                cluster from the list will display the selected cluster components in the 
%                "Select component(s) to plot" list box. Use the plotting buttons below 
%                to plot selected measures of the selected cluster. Additional editing 
%                options (renaming the cluster, rejecting outliers, moving components to 
%                another cluster) are also available. The option 'All N cluster centroids' 
%                at the top of the list displays all the clusters in the list except the 
%                'Notcluster', 'Outlier' and 'ParentCluster' clusters. Selecting this option 
%                will plot the cluster centroids (i.e. ERP, ERSP, ...) in a single figure.
%  "Select component(s) to plot" - [list box] Displays the ICA components of the currently 
%                selected cluster (in the "Select cluster to plot" list box). Each component 
%                has the format: 'subject name, component index'. Multiple components can be 
%                selected from the list. Use the plotting buttons below to plot different 
%                measures of the selected components on different figures. Selecting the 
%                "All components" option is  equivalent to using the cluster plotting buttons. 
%                Additional editing options are reassigning the selected components to 
%                another cluster or moving them to the outlier cluster.
%  "Plot Cluster properties" - [button] Displays in one figure all the mean cluster measures
%                (e.g., dipole locations, scalp maps, spectra, etc.) that were calculated
%                and saved in the EEG datsets. If there is more than one condition, the ERP 
%                and the spectrum will have different colors for each condition. The ERSP 
%                and ITC plots will show only the first condition; clicking on the subplot 
%                will open a new figure with the different conditions displayed together. 
%                Uses the command line function cls_plotclust().
%  "Plot scalp maps"  - [button] Displays the scalp maps of cluster components.
%                If applied to a cluster, scalp maps of the cluster components
%                are plotted along with the cluster mean scalp map in one figure. 
%                If "All # cluster centroids" option is selected, all cluster scalp map
%                means are plotted in the same figure. If applied to components, displays
%                the scalp maps of the specified cluster components in separate figures.
%                Uses the command line functions cls_plotclustmap() and cls_plotcompmap().
%  "Plot ERSPs" - [button] Displays the cluster component ERSPs. 
%                If applied to a cluster, component ERSPs are plotted in one figure  
%                (per condition) with the cluster mean ERSP. If "All # cluster centroids" 
%                option is selected, plots all average ERSPs of the clusters in one figure 
%                per condition. If applied to components, display the ERSP images of specified 
%                cluster components in separate figures, using one figure for all conditions.
%                Uses the command line functions cls_plotclustersp() and cls_plotcompersp().
%  "Plot ITCs" - [button] Same as  "Plot ERSPs" but with ITC.
%                Uses the command line functions cls_plotclustitc() and cls_plotcompitc().
%  "Plot dipoles" - [button] Displays the dipoles of the cluster components.
%                If applied to a cluster, plots the cluster component dipoles (in blue) 
%                plus the average cluster dipole (in red). If "All # cluster centroids" option 
%                is selected, all cluster plots are displayed in one figure each cluster in 
%                a separate subplot. If applied to components, displays the ERSP images of the
%                specified cluster. For specific components displays components dipole (in blue) 
%                plus the average cluster dipole (in Red) in separate figures. 
%                Uses the command line functions cls_plotclustdip() and cls_plotcompdip().
%  "Plot spectra" - [button] Displays the cluster component spectra.   
%                If applied to a cluster, displays component spectra plus the average cluster 
%                spectrum in bold. For a specific cluster, displays the cluster component 
%                spectra plus the average cluster spectrum (in bold) in one figure per condition.
%                If the "All # cluster centroids" option is selected, displays the average 
%                spectrum of all clusters in the same figure, with spectrum for different 
%                conditions (if any) plotted in different colors.  
%                If applied to components, displays the spectrum of specified cluster 
%                components in separate figures using one figure for all conditions.  
%                Uses the command line functions cls_plotclustspec() and cls_plotcompspec().
%  "Plot ERPs" - [button] Same as "Plot spectra" but for ERPs.
%                Uses the command line functions cls_plotclusterp() and cls_plotcomperp().
%  "Create new cluster" - [button] Creates a new empty cluster.
%                Opens a popup window in which a name for the new cluster can be entered.
%                If no name is given the default name is 'Cls #', where '#' is the next
%                available cluster number. For changes to take place, press the popup 
%                window 'OK' button, else press the 'Cancel' button. After the empty 
%                cluster is created, components can be moved into it using, 
%                'Reassign selected component(s)' (see below). Uses the command line 
%                function cls_createclust().
%  "Rename selected cluster" - [button] Renames a cluster using the selected (mnemonic) name. 
%                Opens a popup window in which a new name for the selected cluster can be 
%                entered. For changes to take place, press the popup window 'OK' button, 
%                else press the 'Cancel' button. Uses the command line function cls_renameclust().
%  "Reject outlier components" - [button] rejects outlier components to an outlier cluster.
%                Opens a popup window to specify the outlier threshold. Move outlier 
%                components that are more than x standard deviations devs from the 
%                cluster centroid to an outlier cluster. For changes to take place, 
%                press the popup window 'OK' button, else press the 'Cancel' button. 
%                Uses the command line function cls_rejectoutliers().
%  "Merge clusters" - [button] Merges several clusters into one cluster.
%                Opens a popup window in which the clusters to merge may be specified 
%                An optional name can be given to the merged cluster. If no name is given, 
%                the default name is 'Cls #', where '#' is the next available cluster number.   
%                For changes to take place, press the popup window 'OK' button, else press
%                the 'Cancel' button. Uses the command line function cls_mergeclust().
%  "Remove selected outlier component(s)" - [button] Moves selected component(s) to the 
%                outlier cluster. The components that will be moved are the ones selected 
%                in the "Select component(s) to plot" list box. Opens a popup window in which 
%                a list of the selected component(s) is presented. For changes to take place,
%                press the popup window 'OK' button, else press the 'Cancel' button. 
%                Uses the command line function cls_moveoutlier().
%  "Reassign selected component(s)"- [button] Moves selected component(s) from one cluster 
%                to another. The components that will reassign are the ones selected in the
%                "Select component(s) to plot" list box. Opens a popup window in which 
%                a list of possible clusters to which to move the selected component(s) is 
%                presented. For changes to take place, press the popup window 'OK' button, 
%                else press the 'Cancel' button. Uses the command line function cls_movecomp().
%  "Save STUDY set to disk" - [check box] Saves the STUDY set structure modified according 
%                to specified user edits to the disk. If no file name is entered will
%                overwrite the current STUDY set file. 
%
%  See also  pop_preclust(), pop_clust().         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN/INC/UCSD, October 11, 2004

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

function STUDY = pop_clustedit(varargin)
icadefs;
if ~isstr(varargin{1})
    if nargin < 2
        error('pop_clustedit(): You must provide ALLEEG and STUDY structures');
    end
    STUDY  = varargin{1};
    ALLEEG = varargin{2};
    clus_comps = 0; % the number of clustered components
    if nargin > 2 % load specific clusters
        cls = varargin{3}; %cluster numbers
        N = length(cls); %number of clusters
        % Check clusters are either from the same level (same parents) or are
        % all leaf clusters.
        % Check all have the same parent
        sameparent = 1;
        for k = 1: N
            % Assess the number of clustered components
            if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8))  & (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
                clus_comps = clus_comps + length(STUDY.cluster(cls(k)).comps);
            end
            if k == 1
                parent = STUDY.cluster(cls(k)).parent;
            else
                if isempty(parent) % if the first cluster was the parent cluster
                    parent = STUDY.cluster(cls(k)).parent;
                end
                % For any other case verify that all clusters have the same parents
                if ~(sum(strcmp(STUDY.cluster(cls(k)).parent, parent)) == length(parent)) % different parent
                    if ~strcmp(STUDY.cluster(cls(k)).parent,'manual') & ~strcmp(parent, 'manual') % if nither is an empty cluster (which was created manually)
                        sameparent = 0; % then the clusters have different parents
                    end
                end
            end
        end
        % If not same parent check if all leaf clusters 
        if ~sameparent
            for k = 1: N %check if all leaves
                 if ~isempty(STUDY.cluster(cls(k)).child) 
                     error('pop_clustedit(): All clusters must be from the same level \n         (i.e., have the same parents or not be child clusters)');
                 end
            end
        end
    else   % load leaf clusters
        sameparent = 1;
        cls = [];
        for k = 2:length(STUDY.cluster)
            if isempty(STUDY.cluster(k).child) 
                if isempty(cls)
                    parent = STUDY.cluster(k).parent;
                elseif ~isempty(STUDY.cluster(k).parent)  | ~isempty(parent) % if not both empty                          
                    % Check if all parents are the same
                    if ~(sum(strcmp(STUDY.cluster(k).parent, parent)) == length(parent)) % different parent
                        if ~strcmp(STUDY.cluster(k).parent,'manual') & ~strcmp(parent, 'manual') 
                            sameparent = 0;
                        end
                    end
                end
                cls = [ cls k];
                if ~strncmpi('Notclust',STUDY.cluster(k).name,8)
                     clus_comps = clus_comps + length(STUDY.cluster(k).comps);
                 end
            end
        end
        N = length(cls); %number of clusters
    end
    
    all_comps = length(STUDY.cluster(1).comps);
    % Hold a list of cluster names and the corresponding number of
    % components in the format: Cluster_name (Y ICs).
    num_cls = 0;
    for k = 1:N
        show_options{k+1} = [STUDY.cluster(cls(k)).name ' (' num2str(length(STUDY.cluster(cls(k)).comps))  ' ICs)'];
        if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8)) & (~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8))  & ...
                (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
            num_cls = num_cls + 1;
        end
    end
	show_options{1} = ['All ' num2str(num_cls) ' cluster centroids'];
    
    show_clust = [ 'pop_clustedit(''showclust'',gcf);'];
	plot_clus_maps = [ 'pop_clustedit(''plotclustmap'',gcf); ']; 
    plot_comp_maps = [ 'pop_clustedit(''plotcompmap'',gcf); ']; 
    plot_clus_ersps = ['pop_clustedit(''plotclustersp'',gcf); '];
    plot_comp_ersps = ['pop_clustedit(''plotcompersp'',gcf); '];
    plot_clus_itcs = ['pop_clustedit(''plotclustitc'',gcf); '];
    plot_comp_itcs = ['pop_clustedit(''plotcompitc'',gcf); '];
    plot_clus_spectra = ['pop_clustedit(''plotclustspec'',gcf); '];
    plot_comp_spectra = ['pop_clustedit(''plotcompspec'',gcf); '];
    plot_clus_erp = ['pop_clustedit(''plotclusterp'',gcf); '];
    plot_comp_erp = ['pop_clustedit(''plotcomperp'',gcf); '];
    plot_clus_dip = ['pop_clustedit(''plotclustdip'',gcf); '];
    plot_comp_dip = ['pop_clustedit(''plotcompdip'',gcf); '];
    plot_clus_sum = ['pop_clustedit(''plotclustsum'',gcf); '];
    plot_comp_sum = ['pop_clustedit(''plotcompsum'',gcf); '];
    rename_clust = ['pop_clustedit(''renameclust'',gcf);']; 
    move_comp = ['pop_clustedit(''movecomp'',gcf);'];
    move_outlier = ['pop_clustedit(''moveoutlier'',gcf);'];
    create_clus = ['pop_clustedit(''createclust'',gcf);'];
    reject_outliers = ['pop_clustedit(''rejectoutliers'',gcf);'];
    merge_clusters = ['pop_clustedit(''mergeclusters'',gcf);'];
    saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
    browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                  'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    
    fig_arg{1}{1} = ALLEEG;
    fig_arg{1}{2} = STUDY;
    fig_arg{1}{3} = cls;
    fig_arg{2} = N;
    
    cluster_on = '';
    if sameparent 
        k = 1;
        lm = 1;
        while (k ~= length(cls)) & (lm == 1) %  make sure to take clustered measures from clusters created by the clustering algorithm
            if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8))  & (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
                lm = k;
            end
            k = k +1;
        end
        for k = 1: length(STUDY.cluster(cls(lm)).preclust.preclustparams)
           cls_measure = STUDY.cluster(cls(lm)).preclust.preclustparams{k}{1};
           switch cls_measure
               case 'spec'
                   cls_measure = 'spectrum';
               case 'erp'
                   cls_measure = 'ERP';
               case 'dipoles'
                   cls_measure = 'dipole';
               case 'ersp'
                   cls_measure = 'ERSP';
               case 'itc'
                   cls_measure = 'ITC';
               case 'scalp'
                   cls_measure = 'map';
           end
           if k == 1
                cluster_on = strcat(cluster_on, cls_measure);
            elseif k ~= length(STUDY.cluster(cls(lm)).preclust.preclustparams)
                cluster_on = strcat(cluster_on, ',  ' ,cls_measure);
            else
                cluster_on = strcat(cluster_on, '& ' ,cls_measure);
           end
        end
    end
    
   [out_param userdat] = inputgui( { [4] [1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1] [0.2 1 1.5 0.3] [1] }, ...
	{ {'style' 'text' 'string' ['Study ''' STUDY.name ''': ' num2str(clus_comps) ' of ' num2str(all_comps) ' components ' ...
            fastif(sameparent, ['clustered on ' cluster_on ], 'from leaf clusters') ] ...
            'FontWeight' 'Bold' 'FontSize' 12 'HorizontalAlignment' 'center'} {} ...
    {'style' 'text' 'string' 'Select cluster to plot' 'FontWeight' 'Bold' } {} {'style' 'text' 'string' 'Select component(s) to plot' 'FontWeight' 'Bold'} ...
    {'style' 'listbox' 'string' show_options 'value' 1 'tag' 'clus_list' 'Callback' show_clust} {}  {'style' 'listbox' 'string' '' 'tag' 'clust_comp' 'max' 3 'min' 1} ... 
	{'style' 'pushbutton' 'string' 'Plot cluster properties' 'Callback' plot_clus_sum} {} ... 
    {'style' 'pushbutton' 'string' 'Plot component properties' 'Callback' plot_comp_sum 'enable' 'off'} ...
    {'style' 'pushbutton' 'string' 'Plot scalp maps' 'Callback' plot_clus_maps} {} {'style' 'pushbutton' 'string' 'Plot scalp map(s)' 'Callback' plot_comp_maps}...
	{'style' 'pushbutton' 'string' 'Plot ERSPs' 'Callback' plot_clus_ersps} {}  {'style' 'pushbutton' 'string' 'Plot ERSP(s)' 'Callback' plot_comp_ersps}...
	{'style' 'pushbutton' 'string' 'Plot ITCs' 'Callback' plot_clus_itcs} {}  {'style' 'pushbutton' 'string' 'Plot ITC(s)' 'Callback' plot_comp_itcs}...
    {'style' 'pushbutton' 'string' 'Plot dipoles' 'Callback' plot_clus_dip} {}  {'style' 'pushbutton' 'string' 'Plot dipole(s)' 'Callback' plot_comp_dip}...
	{'style' 'pushbutton' 'string' 'Plot spectra' 'Callback' plot_clus_spectra} {} {'style' 'pushbutton' 'string' 'Plot spectra' 'Callback' plot_comp_spectra} ...
    {'style' 'pushbutton' 'string' 'Plot ERPs' 'Callback' plot_clus_erp} {} {'style' 'pushbutton' 'string' 'Plot ERP(s)' 'Callback' plot_comp_erp} ...
    {} {'style' 'pushbutton' 'string' 'Create new cluster' 'Callback'  create_clus} {} ...
    {'style' 'pushbutton' 'string' 'Remove selected outlier component(s)' 'Callback' move_outlier} ...
    {'style' 'pushbutton' 'string' 'Rename selected cluster' 'Callback' rename_clust } {} ...
    {'style' 'pushbutton' 'string' 'Reassign selected component(s)' 'Callback' move_comp} ...
    {'style' 'pushbutton' 'string' 'Reject outlier components' 'Callback' reject_outliers} {} {} ...
    {'style' 'pushbutton' 'string' 'Merge clusters' 'Callback' merge_clusters} {} {} {}...
    {'style' 'checkbox'  'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} {'style' 'text' 'string' 'Save STUDY set to disk'} ...
    {'style' 'edit' 'string' '' 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
    {'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} {} },...
	'pophelp(''pop_clustoutput'')', 'View and edit current component clusters -- pop_clustedit()' , fig_arg, 'normal', ...
    [ 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] );
	
   if ~isempty(userdat)
       ALLEEG = userdat{1}{1};
       STUDY = userdat{1}{2};
       % If save updated STUDY to disk
        if out_param{3}
            if ~isempty(out_param{4})
                [filepath filename ext] = fileparts(out_param{4});
                a = sprintf('%s%s%s%s%s%s', 'STUDY = pop_savestudy(STUDY,' , '''filename'', ''', [filename ext], ''', ''filepath'', ''', filepath, ''');' );
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);              
                STUDY = pop_savestudy(STUDY, 'filename', [filename ext], 'filepath', filepath);
            else
                if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
                    a = sprintf('%s%s%s%s%s%s', 'STUDY = pop_savestudy(STUDY,' , '''filename'', ''', STUDY.filename, ''', ''filepath'', ''', STUDY.filepath, ''');' );
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);              
                    STUDY = pop_savestudy(STUDY, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
                else
                    a = sprintf('%s%s%s%s%s%s', 'STUDY = pop_savestudy(STUDY);' );
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);              
                    STUDY = pop_savestudy(STUDY);
                end
           end
       end

   end
 
else
    hdl = varargin{2};  %figure handle
    userdat = get(varargin{2}, 'userdat');    
    ALLEEG = userdat{1}{1};
    STUDY = userdat{1}{2};
    cls = userdat{1}{3};
    
    switch  varargin{1}
        
        case {'plotclustmap', 'plotclustersp','plotclustitc','plotclustspec', 'plotclusterp'}
            plotting_option = varargin{1};
            clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
            if (clus ~= 1 ) % specific cluster option
                if ~isempty(STUDY.cluster(cls(clus-1)).comps)
                    eval(['STUDY = cls_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ',''mode'',''comps'');'  ]);
                     % update Study history
                    a = ['STUDY = cls_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ',''mode'',''comps'');'  ];
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                end
            else % all clusters
                % All clusters does not include 'Notclust' 'ParentCluster' and 'Outliers' clusters. 
                tmpcls = [];
                for k = 1:length(cls) 
                    if ~strncmpi(STUDY.cluster(cls(k)).name,'Notclust',8) & ~strncmpi(STUDY.cluster(cls(k)).name,'Outliers',8) & ...
                            (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) & ~isempty(STUDY.cluster(cls(k)).comps)
                        tmpcls = [ tmpcls cls(k)];
                    end
                end
                eval(['STUDY = cls_' plotting_option '(STUDY,ALLEEG,''clusters'',[' num2str(tmpcls) '], ''mode'',''centroid'');'  ]);
                % update Study history
                a = ['STUDY = cls_' plotting_option '(STUDY,ALLEEG,''clusters'',['  num2str(tmpcls) '],''mode'',''centroid'');'  ];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            end
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); 
        case   'plotclustdip'  
            clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
            if (clus ~= 1 ) % specific cluster option
                if ~isempty(STUDY.cluster(cls(clus-1)).comps)
                    [STUDY] = cls_plotclustdip(STUDY, ALLEEG, 'clusters', cls(clus-1), 'mode', 'apart');
                    % update Study history
                    a = ['STUDY = cls_plotclustdip(STUDY, ALLEEG, ''clusters'','  num2str(cls(clus-1)) ',''mode'',''apart'' );'  ];
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                end
            else % all clusters
                % All clusters does not include 'Notclust' and 'Outliers' clusters. 
                tmpcls = [];
                for k = 1:length(cls) 
                    if ~strncmpi(STUDY.cluster(cls(k)).name,'Notclust',8) & ~strncmpi(STUDY.cluster(cls(k)).name,'Outliers',8) & ...
                            (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))  & ~isempty(STUDY.cluster(cls(k)).comps)                           
                        tmpcls = [ tmpcls cls(k)];
                    end
                end
                [STUDY] = cls_plotclustdip(STUDY, ALLEEG, 'clusters', tmpcls, 'mode', 'joined');
                % update Study history
                a = ['STUDY = cls_plotclustdip(STUDY, ALLEEG, ''clusters'',[' num2str(tmpcls)  '],''mode'',''joined'' );'  ];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            end
           userdat{1}{2} = STUDY;
           set(hdl, 'userdat',userdat);    
        case {'plotcompmap', 'plotcompersp','plotcompitc','plotcompspec', 'plotcomperp','plotcompdip'}
            plotting_option = varargin{1};
            clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
            comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
            if (clus ~= 1 ) %specific cluster
                if comp_ind(1) ~= 1  % check that not all comps in cluster are requested
                    eval(['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(clus-1)) ',[' num2str(comp_ind-1) '] );'  ]);
                    % update Study history
                    a = ['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(clus-1)) ',[' num2str(comp_ind-1) '] );'  ];
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                 else
                    eval(['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(clus-1)) ' );'  ]);
                    % update Study history
                    a = ['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(clus-1)) ' );'  ];
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                    if length(comp_ind) > 1 % plot specific components too
                        eval(['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(clus-1)) ',[' num2str(comp_ind(2:end)-1) '] );'  ]);
                        % update Study history
                        a = ['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(clus-1)) ',[' num2str(comp_ind(2:end)-1) '] );' ];
                        STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                    end
                end
            else % all clusters - plot average scalp map
               comp_list = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String');
               comp_name = comp_list(comp_ind);
               for ci = 1:length(comp_name)
                   num_comps = 0;
                   tmp = strfind(comp_name{ci},'''');
                   clust_name = comp_name{ci}(tmp(1)+1:tmp(end)-1);
                   for k = 1:length(cls)
                       if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) & ...
                            (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) 
                           if strcmpi(STUDY.cluster(cls(k)).name, clust_name)
                               cind = comp_ind(ci) - num_comps; % component index in the cluster
                               eval(['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(k)) ',[' num2str(cind) '] );'  ]);
                               % update Study history
                               a = ['STUDY = cls_' plotting_option '(STUDY,ALLEEG,'  num2str(cls(k)) ',[' num2str(cind) '] );' ];
                               STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                               break;
                           else
                               num_comps = num_comps + length(STUDY.cluster(cls(k)).comps);
                           end
                       end                       
                   end
               end
            end
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); 
            
        case 'showclust'
            cind = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
            N = userdat{2};
            count = 1;
            if cind ~= 1 %specific cluster
                len = length(STUDY.cluster(cls(cind-1)).comps);
                compid = cell(len+1,1);
                compid{1} = 'All components';
                % Convert from components numbering to the indexing form 'setXcomY'
                for l = 1:len % go over the components of the cluster
                    subject = STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls(cind-1)).sets(1,l))).subject;
                    compid{l+1} = [  subject ' IC' num2str(STUDY.cluster(cls(cind-1)).comps(1,l)) ];
                end
            else % All clusters accept 'Notclust' and 'Outliers'
                count = 1;
                for k = 1: length(cls)
                    if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) & ...
                            (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) 
                        for l = 1: length(STUDY.cluster(cls(k)).comps)
                            subject = STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cls(k)).sets(1,l))).subject;
                            compid{count} = [ 'Cluster ''' STUDY.cluster(cls(k)).name ''' comp. ' num2str(l) ' (' subject  ' IC' num2str(STUDY.cluster(cls(k)).comps(l)) ')'];
                            count = count +1;
                        end
                    end
                end
            end
           set(findobj('parent', hdl, 'tag', 'clust_comp'), 'String', compid, 'value', 1);
           
        case 'plotclustsum'
            clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
            if clus ~= 1 % specific cluster option
                [STUDY] = cls_plotclust(STUDY, ALLEEG, cls(clus-1));
                % update Study history
                a = ['STUDY = cls_plotclust(STUDY, ALLEEG, '  num2str(cls(clus-1)) ' );'  ];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            else % all clusters
                % All clusters does not include 'Notclust' and 'Outliers' clusters. 
                tmpcls = [];
                for k = 1:length(cls) 
                    if ~strncmpi(STUDY.cluster(cls(k)).name,'Notclust',8) & ~strncmpi(STUDY.cluster(cls(k)).name,'Outliers',8) & ...
                            (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) 
                        tmpcls = [tmpcls cls(k)];
                    end
                end
                [STUDY] = cls_plotclust(STUDY, ALLEEG, tmpcls);
                % update Study history
                a = ['STUDY = cls_plotclust(STUDY, ALLEEG, ['  num2str(tmpcls) '] );'  ];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            end
           userdat{1}{2} = STUDY;
           set(hdl, 'userdat',userdat);    
           
        case 'plotcompsum'
            comps_to_disp = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String');
            clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
            comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
            for ci = 1 : length(comp_ind)
            end
        case 'renameclust'
            clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
            clus_num = get(findobj('parent', hdl, 'tag', 'clus_list'), 'Value') -1;
            if clus_num == 0  % 'all clusters' option 
                return;
            end
            % Don't rename 'Notclust' and 'Outliers'  clusters.
            if strncmpi('Notclust',STUDY.cluster(cls(clus_num)).name,8) | strncmpi('Outliers',STUDY.cluster(cls(clus_num)).name,8) | ...
                    strncmpi('ParentCluster',STUDY.cluster(cls(clus_num)).name,13)
                warndlg2('The ParentCluster, Outliers, and Notclust clusters cannot be renamed');
                return;
			end
            old_name = STUDY.cluster(cls(clus_num)).name;
            rename_param  = inputgui( { [1] [1] [1]}, ...
                { {'style' 'text' 'string' ['Rename ' old_name] 'FontWeight' 'Bold'} {'style' 'edit' 'string' '' 'tag' 'clus_rename' } {} }, ...
            [], 'Rename cluster - from pop_clustedit()' );
            if ~isempty(rename_param) %if not canceled
                new_name = rename_param{1};
                STUDY = cls_renameclust(STUDY, ALLEEG, cls(clus_num), new_name);
                new_name = STUDY.cluster(cls(clus_num)).name;
                % update Study history
                a = ['STUDY = cls_renameclust(STUDY, ALLEEG, ' num2str(cls(clus_num)) ', ' new_name ');'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                clus_name_list{clus_num+1} = [new_name ' (' num2str(length(STUDY.cluster(cls(clus_num)).comps))  ' ICs)'];
                set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                set(findobj('parent', hdl, 'tag', 'clus_rename'), 'String', '');
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update STUDY
            end
            
        case 'movecomp'
            old_clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value') -1;
            comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
            if old_clus == 0 % 'all clusters' option 
                return;
            end
            % Don't reassign components of 'Notclust' or the 'ParentCluster'.
            if strncmpi('Notclust',STUDY.cluster(cls(old_clus)).name,8) | strncmpi('ParentCluster',STUDY.cluster(cls(old_clus)).name,13)  
                warndlg2('Cannot reassign components of ''ParentCluster'' or of ''Notclust''.');
                return;
			end
            old_name = STUDY.cluster(cls(old_clus)).name;
            ncomp = length(comp_ind); % number of selected components
            optionalcls =[];
            for k = 1:length(cls)
                if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8)) & (~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8)) & ...
                        (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))  & (k~= old_clus)
                    optionalcls = [optionalcls cls(k)];
                end
            end                    
            reassign_param  = inputgui( { [1] [1] [1]}, ...
                { {'style' 'text' 'string' ['Reassign ' fastif(ncomp >1, [num2str(length(comp_ind)) ' currently selected components '], 'currently selected component') ...
                            ' from ' old_name ' to the cluster selected below'] 'FontWeight' 'Bold'} ...
                  {'style' 'listbox' 'string' {STUDY.cluster(optionalcls).name} 'tag' 'new_clus'} {} }, ...
                  [], 'Reassign cluster - from pop_clustedit()' ,[] , 'normal', [1 3 1] );
            if ~isempty(reassign_param) %if not canceled
                new_clus = reassign_param{1};
                comp_to_disp = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String');      
                if strcmp(comp_to_disp{comp_ind(1)},'All components')
                    warndlg2('Cannot move all the components of the cluster - abort move components', 'Aborting move components');
                    return;
                end
                STUDY = cls_movecomp(STUDY, ALLEEG,  cls(old_clus), optionalcls(new_clus), comp_ind - 1);                
                % update Study history
                a = ['STUDY = cls_movecomp(STUDY, ALLEEG, ' num2str(cls(old_clus)) ', ' num2str(optionalcls(new_clus)) ', [' num2str(comp_ind - 1) ']);'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                newind = find(cls == optionalcls(new_clus));
                clus_name_list{newind+1} = [STUDY.cluster(optionalcls(new_clus)).name ' (' num2str(length(STUDY.cluster(optionalcls(new_clus)).comps))  ' ICs)'];
                clus_name_list{old_clus+1} = [STUDY.cluster(cls(old_clus)).name ' (' num2str(length(STUDY.cluster(cls(old_clus)).comps))  ' ICs)'];
                set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 
                pop_clustedit('showclust',hdl);
            end          
            
        case 'moveoutlier'
            old_clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value') -1;
            comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
            if ~isempty(find(comp_ind ==1))
                warndlg2('Cannot remove all the cluster components');
                return;
            end
            if old_clus == 0 % 'all clusters' option 
                return;
            end
            if strncmpi('Notclust',STUDY.cluster(cls(old_clus)).name,8) | strncmpi('ParentCluster',STUDY.cluster(cls(old_clus)).name,13)    % There are no outliers to 'Notclust'
                warndlg2('Cannot reassign components of ''Notclust'' or ''ParentCluster''.');
                return;
            end
            comp_list = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String'); 
            ncomp = length(comp_ind);
            old_name = STUDY.cluster(cls(old_clus)).name;            
            if strncmpi('Outliers',STUDY.cluster(cls(old_clus)).name,8)  % There are no outliers of 'Outliers'
                warndlg2('Cannot use ''Outliers'' clusters for this option.');
                return;
			end
            reassign_param  = inputgui( { [1] [1] [1]}, ...
                { {'style' 'text' 'string' ['Remove ' fastif(ncomp >1, [num2str(length(comp_ind)) ' currently selected components below '], 'currently selected component below ') ...
                            'from ' old_name ' to its outlier cluster?'] 'FontWeight' 'Bold'} ...
                  {'style' 'listbox' 'string' {comp_list{comp_ind}} 'tag' 'new_clus'} {} }, ...
                  [], 'Remove outliers - from pop_clustedit()' ,[] , 'normal', [1 3 1] );
            if ~isempty(reassign_param) %if not canceled
                STUDY = cls_moveoutlier(STUDY, ALLEEG,  cls(old_clus), comp_ind - 1);
                clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                outlier_clust = cls_findoutlierclust(STUDY,cls(old_clus)); %find the outlier cluster for this cluster
                oind = find(cls == outlier_clust); % the outlier clust index (if already exist) in the cluster list GUI
                if ~isempty(oind) % the outlier clust is already presented in the cluster list GUI
                    clus_name_list{oind+1} = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                elseif outlier_clust == length(STUDY.cluster) % update the list with the Outlier cluster (if didn't exist before)
                    clus_name_list{end+1} = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                    userdat{2} = userdat{2} + 1; % update N, number of clusters in edit window 
                    cls(end +1) = length(STUDY.cluster); % update the GUI clusters list with the outlier cluster
                    userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                end
                clus_name_list{old_clus+1} = [STUDY.cluster(cls(old_clus)).name ' (' num2str(length(STUDY.cluster(cls(old_clus)).comps))  ' ICs)'];
                set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                % update Study history
                a = ['STUDY = cls_moveoutlier(STUDY, ALLEEG, ' num2str(cls(old_clus)) ', [' num2str(comp_ind - 1) ']);'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 
                pop_clustedit('showclust',hdl);    
            end
            
        case 'rejectoutliers'
            clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'Value') -1;
            if clus
                cls_name = STUDY.cluster(cls(clus)).name;
                % Cannot reject outliers from 'Notclust', 'ParentCluster' and 'Outlier' clusters
                if strncmpi('Notclust',cls_name,8) | strncmpi('ParentCluster', cls_name,13) | ...
                        strncmpi('Outliers',cls_name,8)
                    warndlg2('Cannot reject outliers of ''Notclust'' or ''Outliers'' or ''ParentCluster'' clusters.');
                    return;
			    end
                clusters = cls(clus);
            else
                cls_name = 'All clusters';
                clusters = [];
                for k = 1:length(cls)
                     if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) & ...
                             ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)  
                        clusters = [ clusters cls(k)];
                    end
                end
            end
            reject_param  = inputgui( { [1] [2 1 3] [1]}, ...
                { {'style' 'text' 'string' ['Reject ' cls_name  ' outliers ' ] 'FontWeight' 'Bold'} ...
                   {'style' 'text' 'string' 'Move outlier components that are more than'} {'style' 'edit' 'string' '3' 'tag' 'outliers_std' } ...
                   {'style' 'text' 'string' [ 'standard deviations devs from the ' cls_name  ' centroid to an outlier cluster.']} {} }, ...
                  [], 'Reject outliers - from pop_clustedit()' );
            if ~isempty(reject_param) %if not canceled
                ostd = reject_param{1}; % the requested outlier std
                [STUDY] = cls_rejectoutliers(STUDY, ALLEEG, clusters, str2num(ostd));  
                % update Study history
                a = ['STUDY = cls_rejectoutliers(STUDY, ALLEEG, [ ' num2str(clusters) ' ], ' ostd ');'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                for k = 1:length(clusters)
                    outlier_clust = cls_findoutlierclust(STUDY,clusters(k)); %find the outlier cluster for this cluster
                    oind = find(cls == outlier_clust); % the outlier clust index (if already exist) in the cluster list GUI
                    if ~isempty(oind) % the outlier clust is already presented in the cluster list GUI
                        clus_name_list{oind+1} = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                    else % update the list with the outlier cluster 
                        clus_name_list{end+1} = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                        userdat{2} = userdat{2} + 1; % update N, number of clusters in edit window 
                        cls(end +1) = outlier_clust; % update the GUI clusters list with the outlier cluster
                        userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                    end
                    clsind = find(cls == clusters(k));
                    clus_name_list{clsind+1} = [STUDY.cluster(clusters(k)).name ' (' num2str(length(STUDY.cluster(clusters(k)).comps))  ' ICs)'];
                    set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                end
                % If outlier cluster doesn't exist in the GUI window add it 
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 
                pop_clustedit('showclust',hdl);
            end
            
        case 'createclust'
            create_param  = inputgui( { [1] [1 1] [1]}, ...
                { {'style' 'text' 'string' 'Create new empty cluster' 'FontWeight' 'Bold'} ...
                   {'style' 'text' 'string' 'Enter cluster name:'} {'style' 'edit' 'string' '' } {} }, ...
                  [], 'Create new empty cluster - from pop_clustedit()' );
            if ~isempty(create_param) %if not canceled
                clus_name = create_param{1}; % the name of the new cluster
                [STUDY] = cls_createclust(STUDY, ALLEEG, clus_name); 
                % Update cluster list
                clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                clus_name_list{end+1} = [STUDY.cluster(end).name ' (0 ICs)']; %update the cluster list with the new cluster
                % update the first option on the GUI list : 'All 10 cluster centroids'
                % with the new number of cluster centroids
                ti = strfind(clus_name_list{1},'cluster'); %get the number of clusters centroid 
                cent = num2str(str2num(clus_name_list{1}(5:ti-2))+1); % new number of centroids
                clus_name_list{1} = [clus_name_list{1}(1:4) cent clus_name_list{1}(ti-1:end)]; %update list
                set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                % update Study history
                if isempty(clus_name)
                    a = ['STUDY = cls_createclust(STUDY, ALLEEG);'];
                else
                    a = ['STUDY = cls_createclust(STUDY, ALLEEG, ' clus_name ');'];
                end                
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                userdat{1}{2} = STUDY;
                userdat{2} = userdat{2} + 1; % update N, the number of cluster options in edit window 
                cls(end +1) = length(STUDY.cluster); % update the GUI clusters list with the new cluster
                userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                set(hdl, 'userdat',userdat); %update STUDY, cls and N
            end
            
        case 'mergeclusters'
            clus_names = get(findobj('parent', hdl, 'tag', 'clus_list'), 'string') ;
            optionalcls =[];
            for k = 2:length(clus_names)
                if (~strncmpi('Notclust',clus_names{k},8)) & (~strncmpi('Outliers',clus_names{k},8)) & ...
                        (~strncmpi('ParentCluster',clus_names{k},13))
                    optionalcls = [optionalcls k];
                end
            end           
            reassign_param  = inputgui( { [1] [1] [1] [2 1] [1]}, ...
                { {'style' 'text' 'string' 'Select clusters to Merge' 'FontWeight' 'Bold'} ...
                  {'style' 'listbox' 'string' clus_names(optionalcls) 'tag' 'new_clus' 'max' 3 'min' 1} {} ...
                  {'style' 'text' 'string' 'Optional, enter a name for the merged cluster:' 'FontWeight' 'Bold'} ...
                  {'style' 'edit' 'string' ''} {} }, ...
                  [], 'Merge clusters - from pop_clustedit()' ,[] , 'normal', [1 3 1 1 1] );
              if ~isempty(reassign_param)
                  cls_mrg = cls(optionalcls(reassign_param{1})-1);
                  name = reassign_param{2};
                  allleaves = 1;
                  N = userdat{2};
                  for k = 1: N %check if all leaves
                      if ~isempty(STUDY.cluster(cls(k)).child) 
                          allleaves = 0;
                      end
                  end                     
                  [STUDY] = cls_mergeclust(STUDY, ALLEEG, cls_mrg, name); 
                  % 
                  % update Study history
                  % 
                  if isempty(name)
                      a = ['STUDY = cls_mergeclust(STUDY, ALLEEG, [' num2str(cls_mrg) ']);'];
                  else
                      a = ['STUDY = cls_mergeclust(STUDY, ALLEEG, [' num2str(cls_mrg) '], ' name ');'];
                  end                  
                  STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                  userdat{1}{2} = STUDY;
                  %
                  % Replace the merged clusters with the one new merged cluster 
                  % in the GUI if all clusters are leaves
                  %
                  if allleaves                      
                    %
                    % Update cluster list
                    %
                    clus_names{end+1} = [STUDY.cluster(end).name ' (' num2str(length(STUDY.cluster(end).comps))  ' ICs)']; 
                    %
                    % update the cluster list with the new cluster
                    %
                    clus_names([optionalcls(reassign_param{1})]) = [];
                    cls = setdiff(cls, cls_mrg); % remove from the GUI clusters list the merged clusters
                    cls(end+1) = length(STUDY.cluster); % update the GUI clusters list with the new cluster
                    N  = length(cls);
                    %
                    % update the first option on the GUI list : 'All 10 cluster centroids'
                    % with the new number of cluster centroids
                    %
                    ti = strfind(clus_names{1},'cluster'); %get the number of clusters centroid 
                    cent = num2str(str2num(clus_names{1}(5:ti-2))+1- length(cls_mrg)); % new number of centroids
                    clus_names{1} = [clus_names{1}(1:4) cent clus_names{1}(ti-1:end)]; %update list
                    set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_names);
                    %
                    % update Study history
                    %
                    userdat{2} = N; % update N, the number of cluster options in edit window 
                    userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                  end
                  set(hdl, 'userdat',userdat); %update information (STUDY)     
                  pop_clustedit('showclust',hdl);
              end
    end
end
