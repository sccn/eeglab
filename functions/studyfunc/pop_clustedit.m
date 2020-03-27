% pop_clustedit() - graphic user interface (GUI)-based function with editing and plotting 
%                options for visualizing and manipulating an input STUDY structure. 
%                Only component measures (e.g., dipole locations, scalp maps, spectra, 
%                ERPs, ERSPs, ITCs) that have been computed and saved in the study EEG 
%                datasets can be visualized. These can be computed during pre-clustering 
%                using the GUI-based function pop_preclust() or the equivalent command
%                line functions std_preclust(). To use dipole locations for clustering, 
%                they must first be stored in the EEG dataset structures using dipfit(). 
%                Supported cluster editing functions include new cluster creation, cluster 
%                merging, outlier rejection, and cluster renaming. Components can also be 
%                moved from one cluster to another or to the outlier cluster. 
% Usage:    
%                >> STUDY = pop_clustedit(STUDY, ALLEEG, clusters, addui, addgeom);   
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
%   addui      - [struct] additional uicontrols entries for the graphic
%                interface. Must contains the fiels "uilist", "geometry".
%
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user edits,
%                if any. Plotted cluster measure means (maps, ERSPs, etc.) are added to 
%                the STUDY structure after they are first plotted to allow quick replotting.  
%
% Graphic interface buttons:
%  "Select cluster to plot" - [list box] Displays available clusters to plot (format is
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
%                Uses the command line function std_propplot().
%  "Plot scalp maps"  - [button] Displays the scalp maps of cluster components.
%                If applied to a cluster, scalp maps of the cluster components
%                are plotted along with the cluster mean scalp map in one figure. 
%                If "All # cluster centroids" option is selected, all cluster scalp map
%                means are plotted in the same figure. If applied to components, displays
%                the scalp maps of the specified cluster components in separate figures.
%                Uses the command line functions std_topoplot().
%  "Plot ERSPs" - [button] Displays the cluster component ERSPs. 
%                If applied to a cluster, component ERSPs are plotted in one figure  
%                (per condition) with the cluster mean ERSP. If "All # cluster centroids" 
%                option is selected, plots all average ERSPs of the clusters in one figure 
%                per condition. If applied to components, display the ERSP images of specified 
%                cluster components in separate figures, using one figure for all conditions.
%                Uses the command line functions std_erspplot().
%  "Plot ITCs" - [button] Same as  "Plot ERSPs" but with ITC.
%                Uses the command line functions std_itcplot().
%  "Plot dipoles" - [button] Displays the dipoles of the cluster components.
%                If applied to a cluster, plots the cluster component dipoles (in blue) 
%                plus the average cluster dipole (in red). If "All # cluster centroids" option 
%                is selected, all cluster plots are displayed in one figure each cluster in 
%                a separate subplot. If applied to components, displays the ERSP images of the
%                specified cluster. For specific components displays components dipole (in blue) 
%                plus the average cluster dipole (in Red) in separate figures. 
%                Uses the command line functions std_dipplot().
%  "Plot spectra" - [button] Displays the cluster component spectra.   
%                If applied to a cluster, displays component spectra plus the average cluster 
%                spectrum in bold. For a specific cluster, displays the cluster component 
%                spectra plus the average cluster spectrum (in bold) in one figure per condition.
%                If the "All # cluster centroids" option is selected, displays the average 
%                spectrum of all clusters in the same figure, with spectrum for different 
%                conditions (if any) plotted in different colors.  
%                If applied to components, displays the spectrum of specified cluster 
%                components in separate figures using one figure for all conditions.  
%                Uses the command line functions std_specplot().
%  "Plot ERPs" - [button] Same as "Plot spectra" but for ERPs.
%                Uses the command line functions std_erpplot().
%  "Plot ERPimage" - [button] Same as "Plot ERP" but for ERPimave.
%                Uses the command line functions std_erpimplot().
%  "Create new cluster" - [button] Creates a new empty cluster.
%                Opens a popup window in which a name for the new cluster can be entered.
%                If no name is given the default name is 'Cls #', where '#' is the next
%                available cluster number. For changes to take place, press the popup 
%                window 'OK' button, else press the 'Cancel' button. After the empty 
%                cluster is created, components can be moved into it using, 
%                'Reassign selected component(s)' (see below). Uses the command line 
%                function std_createclust().
%  "Rename selected cluster" - [button] Renames a cluster using the selected (mnemonic) name. 
%                Opens a popup window in which a new name for the selected cluster can be 
%                entered. For changes to take place, press the popup window 'OK' button, 
%                else press the 'Cancel' button. Uses the command line function std_renameclust().
%  "Reject outlier components" - [button] rejects outlier components to an outlier cluster.
%                Opens a popup window to specify the outlier threshold. Move outlier 
%                components that are more than x standard deviations devs from the 
%                cluster centroid to an outlier cluster. For changes to take place, 
%                press the popup window 'OK' button, else press the 'Cancel' button. 
%                Uses the command line function std_rejectoutliers().
%  "Merge clusters" - [button] Merges several clusters into one cluster.
%                Opens a popup window in which the clusters to merge may be specified 
%                An optional name can be given to the merged cluster. If no name is given, 
%                the default name is 'Cls #', where '#' is the next available cluster number.   
%                For changes to take place, press the popup window 'OK' button, else press
%                the 'Cancel' button. Uses the command line function std_mergeclust().
%  "Remove selected outlier component(s)" - [button] Moves selected component(s) to the 
%                outlier cluster. The components that will be moved are the ones selected 
%                in the "Select component(s) to plot" list box. Opens a popup window in which 
%                a list of the selected component(s) is presented. For changes to take place,
%                press the popup window 'OK' button, else press the 'Cancel' button. 
%                Uses the command line function std_moveoutlier().
%  "Reassign selected component(s)" - [button] Moves selected component(s) from one cluster 
%                to another. The components that will reassign are the ones selected in the
%                "Select component(s) to plot" list box. Opens a popup window in which 
%                a list of possible clusters to which to move the selected component(s) is 
%                presented. For changes to take place, press the popup window 'OK' button, 
%                else press the 'Cancel' button. Uses the command line function std_movecomp().
%  "Save STUDY set to disk" - [check box] Saves the STUDY set structure modified according 
%                to specified user edits to the disk. If no file name is entered will
%                overwrite the current STUDY set file. 
%
% See also:  pop_preclust(), pop_clust().         
%
% Authors: Arnaud Delorme, Hilit Serby, Scott Makeig, SCCN/INC/UCSD, October 11, 2004

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

function [STUDY, com] = pop_clustedit(varargin)
icadefs;

if nargin < 2
    help pop_clustedit;
    return;
end

if ~ischar(varargin{1})
    STUDY  = varargin{1};
    STUDY.etc.erpparams.topotime    = NaN; % [] for channels and NaN for components
    STUDY.etc.specparams.topofreq   = NaN; % NaN -> GUI disabled
    STUDY.etc.erspparams.topotime   = NaN;
    STUDY.etc.erspparams.topofreq   = NaN;
    STUDY.etc.erpimparams.topotime  = NaN;
    STUDY.etc.erpimparams.topotrial = NaN;
    STUDY.tmphist = '';
    ALLEEG = varargin{2};
    clus_comps = 0; % the number of clustered components
    if nargin > 2 && ~isempty(varargin{3}) % load specific clusters
        cls = varargin{3}; %cluster numbers
        N = length(cls); %number of clusters
        
        % Check clusters are either from the same level (same parents) or are
        % all leaf clusters.
        % Check all input clusters have the same parent
        
        sameparent = 1;
        for k = 1: N
            % Assess the number of clustered components
            if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8))  && (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
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
                    if ~strcmp(STUDY.cluster(cls(k)).parent,'manual') && ~strcmp(parent, 'manual') 
                        % if nither is an empty cluster (which was created manually)
                        sameparent = 0; % then the clusters have different parents
                    end
                end
            end
        end
        
        % If not same parent check if all leaf clusters 
        % ---------------------------------------------
        if ~sameparent
            for k = 1: N %check if all leaves
                 if ~isempty(STUDY.cluster(cls(k)).child) 
                     error([ 'pop_clustedit(): All clusters must be from the same level \n' ...
                             '         (i.e., have the same parents or not be child clusters)' ]);
                 end
            end
        end

        % ploting text etc ...
        % --------------------
        num_cls = 0;
        for k = 1:N
            show_options{k+1} = [STUDY.cluster(cls(k)).name ' (' num2str(length(STUDY.cluster(cls(k)).comps))  ' ICs)'];
            if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8)) && (~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8))  && ...
                    (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
                num_cls = num_cls + 1;
            end
        end
        show_options{1} = ['All ' num2str(num_cls) ' cluster centroids'];
        
    else   % load leaf clusters
        
        sameparent = 1;
        cls = [];
        for k = 2:length(STUDY.cluster)
            if isempty(STUDY.cluster(k).child) 
                if isempty(cls)
                    parent = STUDY.cluster(k).parent;
                elseif ~isempty(STUDY.cluster(k).parent)  || ~isempty(parent) % if not both empty                          
                    % Check if all parents are the same
                    if ~(sum(strcmp(STUDY.cluster(k).parent, parent)) == length(parent)) % different parent
                        if ~strcmp(STUDY.cluster(k).parent,'manual') && ~strcmp(parent, 'manual') 
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
        
    
        % Plot clusters hierarchically
        % ----------------------------
        num_cls = 0;
        cls = 1:length(STUDY.cluster);
        N = length(cls); %number of clusters
        
        show_options{1} = [STUDY.cluster(1).name ' (' num2str(length(STUDY.cluster(1).comps))  ' ICs)'];
        cls(1) = 1;
        count = 2;
        for index1 = 1:length(STUDY.cluster(1).child)
            
            indclust1 = strmatch( STUDY.cluster(1).child(index1), { STUDY.cluster.name }, 'exact');
            show_options{count} = ['   ' STUDY.cluster(indclust1).name ' (' num2str(length(STUDY.cluster(indclust1).comps))  ' ICs)'];
            cls(count) = indclust1;
            count = count+1;
            
            for index2 = 1:length( STUDY.cluster(indclust1).child )
                indclust2 = strmatch( STUDY.cluster(indclust1).child(index2), { STUDY.cluster.name }, 'exact');
                show_options{count} = ['      ' STUDY.cluster(indclust2).name ' (' num2str(length(STUDY.cluster(indclust2).comps))  ' ICs)'];
                cls(count) = indclust2;
                count = count+1;
                    
                for index3 = 1:length( STUDY.cluster(indclust2).child )
                    indclust3 = strmatch( STUDY.cluster(indclust2).child(index3), { STUDY.cluster.name }, 'exact');
                    show_options{count} = ['         ' STUDY.cluster(indclust3).name ' (' num2str(length(STUDY.cluster(indclust3).comps))  ' ICs)'];
                    cls(count) = indclust3;
                    count = count+1;
                end
            end
        end
        show_options = { ['All cluster centroids'] show_options{:} }; 
    end
    all_comps = length(STUDY.cluster(1).comps);
    
    show_clust      = [ 'pop_clustedit(''showclust'',gcf);'];
    show_comps      = [ 'pop_clustedit(''showcomplist'',gcf);'];
	plot_clus_maps  = [ 'pop_clustedit(''topoplot'',gcf); ']; 
    plot_comp_maps  = [ 'pop_clustedit(''plotcomptopo'',gcf); ']; 
    plot_clus_ersps = ['pop_clustedit(''erspplot'',gcf); '];
    plot_comp_ersps = ['pop_clustedit(''plotcompersp'',gcf); '];
    plot_clus_itcs  = ['pop_clustedit(''itcplot'',gcf); '];
    plot_comp_itcs  = ['pop_clustedit(''plotcompitc'',gcf); '];
    plot_clus_erpim = ['pop_clustedit(''erpimageplot'',gcf); '];
    plot_comp_erpim = ['pop_clustedit(''plotcomperpimage'',gcf); '];
    plot_clus_spectra = ['pop_clustedit(''specplot'',gcf); '];
    plot_comp_spectra = ['pop_clustedit(''plotcompspec'',gcf); '];
    plot_clus_erp = ['pop_clustedit(''erpplot'',gcf); '];
    plot_comp_erp = ['pop_clustedit(''plotcomperp'',gcf); '];
    plot_clus_dip = ['pop_clustedit(''dipplot'',gcf); '];
    plot_comp_dip = ['pop_clustedit(''plotcompdip'',gcf); '];
    plot_clus_sum = ['pop_clustedit(''plotsum'',gcf); '];
    plot_comp_sum = ['pop_clustedit(''plotcompsum'',gcf); '];
    rename_clust  = ['pop_clustedit(''renameclust'',gcf);']; 
    move_comp     = ['pop_clustedit(''movecomp'',gcf);'];
    move_outlier  = ['pop_clustedit(''moveoutlier'',gcf);'];
    create_clus   = ['pop_clustedit(''createclust'',gcf);'];
    reject_outliers = ['pop_clustedit(''rejectoutliers'',gcf);'];
    merge_clusters = ['pop_clustedit(''mergeclusters'',gcf);'];
    dip_opt        = ['pop_clustedit(''dip_opt'',gcf);'];
    erp_opt        = ['pop_clustedit(''erp_opt'',gcf);'];
    spec_opt       = ['pop_clustedit(''spec_opt'',gcf);'];
    ersp_opt       = ['pop_clustedit(''ersp_opt'',gcf);'];
    erpim_opt      = ['pop_clustedit(''erpim_opt'',gcf);'];
    stat_opt       = ['pop_clustedit(''stat_opt'',gcf);'];
    saveSTUDY      = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
    browsesave     = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    
    % Create default ERSP / ITC time/freq. paramters 
    % ----------------------------------------------
    if isempty(ALLEEG)
        error('STUDY contains no datasets');
    end

    % enable buttons
    % --------------
    filename = fullfile(STUDY.datasetinfo(1).filepath, STUDY.datasetinfo(1).subject);
    if ~isempty(dir([filename '*.icaspec'])),   spec_enable = 'on'; else spec_enable  = 'off'; end
    if ~isempty(dir([filename '*.icaerp'] )) ,   erp_enable = 'on'; else erp_enable   = 'off'; end
    if ~isempty(dir([filename '*.icaerpim'] )), erpim_enable = 'on'; else erpim_enable = 'off'; end
    if ~isempty(dir([filename '*.icatimef'])) ,   ersp_enable = 'on'; else ersp_enable  = 'off'; end
    filename = fullfile( ALLEEG(1).filepath, ALLEEG(1).filename(1:end-4));
    if ~isempty(dir([filename '*.icatopo'])), scalp_enable = 'on'; else scalp_enable = 'off'; end
    
    if isfield(ALLEEG(1).dipfit, 'model'), dip_enable   = 'on'; else dip_enable   = 'off'; end
    
    % userdata below
    % --------------
    fig_arg{1}{1} = ALLEEG;
    fig_arg{1}{2} = STUDY;
    fig_arg{1}{3} = cls;
    fig_arg{2}    = N;
        
    str_name       = sprintf('STUDY ''%s'' - ''%s'' component clusters', STUDY.name, STUDY.design(STUDY.currentdesign).name);
    if length(str_name) > 80, str_name = [ str_name(1:80) '...''' ]; end
    if length(cls) > 1, vallist = 1; else vallist = 2; end
    geomline = [1 0.35 1];
    geometry = { [0.8 3] [1] geomline geomline geomline geomline geomline geomline geomline geomline ...
                 geomline [1] geomline geomline }; %geomline };
    geomvert = [ 1 .5 1 3 1 1 1 1 1 1 1 1 1 1 ]; %1];
    uilist   = { ...
        {'style' 'text'       'string' 'Select design:' 'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} ...
        {'style' 'popupmenu'  'string' { STUDY.design.name } 'FontWeight' 'Bold' 'tag' 'design' 'value' STUDY.currentdesign  } ...
        { } ...
        {'style' 'text'       'string' 'Select cluster to plot' 'FontWeight' 'Bold' } {} ...
        {'style' 'text'       'string' 'Select component to plot        ' 'FontWeight' 'Bold'} ...
        {'style' 'listbox'    'string' show_options 'value' vallist 'tag' 'clus_list' 'Callback' show_clust } ...
        {'style' 'pushbutton' 'enable'   'on'       'string' 'STATS' 'Callback' stat_opt }  ...
        {'style' 'listbox'    'string' '' 'tag' 'clust_comp' 'max' 2 'min' 1 'callback'    show_comps } ... 
        {'style' 'pushbutton' 'enable' scalp_enable 'string' 'Plot scalp maps' 'Callback' plot_clus_maps} {} ...
        {'style' 'pushbutton' 'enable' scalp_enable 'string' 'Plot scalp map(s)' 'Callback' plot_comp_maps}...
        {'style' 'pushbutton' 'enable'   dip_enable 'string' 'Plot dipoles' 'Callback' plot_clus_dip}  ...
        {'style' 'pushbutton' 'enable'   dip_enable 'string' 'Params' 'Callback' dip_opt }  ...
        {'style' 'pushbutton' 'enable'   dip_enable 'string' 'Plot dipole(s)' 'Callback' plot_comp_dip}...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Plot ERPs' 'Callback' plot_clus_erp} ...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Params' 'Callback' erp_opt }  ...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Plot ERP(s)' 'Callback' plot_comp_erp} ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Plot spectra' 'Callback' plot_clus_spectra} ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Params' 'Callback' spec_opt } ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Plot spectra' 'Callback' plot_comp_spectra} ...
        {'style' 'pushbutton' 'enable' erpim_enable 'string' 'Plot ERPimage' 'Callback' plot_clus_erpim} ...
        {'style' 'pushbutton' 'enable' erpim_enable 'string' 'Params' 'Callback' erpim_opt } ...
        {'style' 'pushbutton' 'enable' erpim_enable 'string' 'Plot ERPimage(s)' 'Callback' plot_comp_erpim} ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ERSPs' 'Callback' plot_clus_ersps} ...
        {'vertexpand' 2.1 'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Params' 'Callback' ersp_opt } ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ERSP(s)' 'Callback' plot_comp_ersps} ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ITCs' 'Callback' plot_clus_itcs} { } ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ITC(s)' 'Callback' plot_comp_itcs} ...
        ... % {} {}... %{'style' 'pushbutton' 'string' 'Plot cluster properties' 'Callback' plot_clus_sum 'enable' 'off'} {} ... 
        ... % {'style' 'pushbutton' 'string' 'Plot component properties' 'Callback' plot_comp_sum 'enable' 'on'} ... % nima, was off
        {} ...
        {'style' 'pushbutton' 'string' 'Create new cluster' 'Callback'  create_clus} {} ...
        {'style' 'pushbutton' 'string' 'Reassign selected component(s)' 'Callback' move_comp} ...
        {'style' 'pushbutton' 'string' 'Rename selected cluster' 'Callback' rename_clust } {} ...
        {'style' 'pushbutton' 'string' 'Remove selected outlier comps.' 'Callback' move_outlier} };
%        {'style' 'pushbutton' 'string' 'Merge clusters' 'Callback' merge_clusters } {} ...
%        {'style' 'pushbutton' 'string' 'Auto-reject outlier components' 'Callback' reject_outliers } };
   
   % additional UI given on the command line
   % ---------------------------------------
   if nargin > 3
       addui = varargin{4};
       if ~isfield(addui, 'uilist')
           error('Additional GUI definition (argument 4) requires the field "uilist"');
       end
       if ~isfield(addui, 'geometry')
           addui.geometry = mat2cell(ones(1,length(addui.uilist)));
       end
       uilist = { uilist{:}, addui.uilist{:} };
       geometry = { geometry{:} addui.geometry{:} };
       geomvert = [ geomvert ones(1,length(addui.geometry)) ];
   end
   
   [out_param, userdat] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''pop_clustoutput'')', ...
                                   'title', 'View and edit current component clusters -- pop_clustedit()' , 'userdata', fig_arg, ...
                                   'geomvert', geomvert, 'eval', show_clust );
	
   if ~isempty(userdat)
       ALLEEG = userdat{1}{1};
       STUDY = userdat{1}{2};
   end
   
   % history
   % -------
   com = STUDY.tmphist;
   STUDY = rmfield(STUDY, 'tmphist');
 
else
    hdl = varargin{2};  %figure handle
    userdat = get(varargin{2}, 'userdat');    
    ALLEEG  = userdat{1}{1};
    STUDY   = userdat{1}{2};
    cls     = userdat{1}{3};
    design  = get(findobj('parent', hdl, 'tag', 'design')      , 'value');
	if ~std_checkdesign(STUDY, design)
        return;
    end
    
    clus     = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
    comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
    if clus == 1 && length(cls) == 1
        warndlg2('No cluster', 'No cluster');
        return;
    end
    
    try
        switch  varargin{1}

            case {'plotcomptopo', 'plotcompersp','plotcompitc','plotcompspec', 'plotcomperp', 'plotcompdip', 'plotcomperpimage'}
                plotting_option = varargin{1};
                plotting_option = [ plotting_option(9:end) 'plot' ];
                if (clus ~= 1 ) %specific cluster
                    if comp_ind(1) ~= 1  % check that not all comps in cluster are requested
                        subject = STUDY.datasetinfo( STUDY.cluster(cls(clus-1)).sets(1,comp_ind-1)).subject;
                        a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ', ''comps'', ' num2str(comp_ind-1) ', ''design'', ' int2str(design) ' );' ];
                        eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                     else
                        a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ', ''design'', ' int2str(design) ', ''plotsubjects'', ''on'' );' ];
                        eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                    end
                else
                   comp_list = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String');
                   comp_name = comp_list(comp_ind);
                   for ci = 1:length(comp_name)
                       num_comps = 0;
                       tmp = strfind(comp_name{ci},'''');
                       clust_name = comp_name{ci}(tmp(1)+1:tmp(end)-1);
                       for k = 1:length(cls)
                           if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) && ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) && ...
                                (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) 
                               if strcmpi(STUDY.cluster(cls(k)).name, clust_name)
                                   cind = comp_ind(ci) - num_comps; % component index in the cluster
                                   subject = STUDY.datasetinfo( STUDY.cluster(cls(k)).sets(1,cind)).subject;
                                   a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(k)) ', ''design'', ' int2str(design) ', ''comps'',' num2str(cind) ' );' ];
                                   eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
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

            case {'topoplot', 'erspplot', 'itcplot', 'specplot', 'erpplot', 'dipplot', 'erpimageplot' }
                plotting_option = varargin{1};
                plotting_option = [ plotting_option(1:end-4) 'plot' ];
                if (clus ~= 1 ) % specific cluster option
                    if ~isempty(STUDY.cluster(cls(clus-1)).comps)
                        a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ', ''design'', ' int2str(design) ');' ];
                        eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                    end
                else % all clusters
                    % All clusters does not include 'Notclust' 'ParentCluster' and 'Outliers' clusters. 
                    tmpcls = [];
                    for k = 1:length(cls) 
                        if ~strncmpi(STUDY.cluster(cls(k)).name,'Notclust',8) && ~strncmpi(STUDY.cluster(cls(k)).name,'Outliers',8) && ...
                                (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) && ~isempty(STUDY.cluster(cls(k)).comps)
                            tmpcls = [ tmpcls cls(k)];
                        end
                    end
                    a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'',['  num2str(tmpcls) '], ''design'', ' int2str(design) ');' ];
                    %if strcmpi(plotting_option, 'dipplot'), a = [a(1:end-2) ',''mode'', ''together'');' ]; end
                    eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 

            case 'dip_opt' % save the list of selected chaners
                [STUDY com] = pop_dipparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'erp_opt' % save the list of selected chaners
                [STUDY com] = pop_erpparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'stat_opt' % save the list of selected chaners
                [STUDY com] = pop_statparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'spec_opt' % save the list of selected channels
                [STUDY com] = pop_specparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'erpim_opt' % save the list of selected channels
                [STUDY com] = pop_erpimparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'ersp_opt' % save the list of selected channels
                [STUDY com] = pop_erspparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

           case 'showcomplist' % save the list of selected clusters
                clust = get(findobj('parent', hdl, 'tag', 'clus_list') , 'value');
                comp  = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'value');
                N = userdat{2};
                count = 1;
                if clust ~= 1 %specific cluster
                    STUDY.cluster(cls(clust-1)).selected = comp;
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

           case 'showclust'
                cind = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
                N = userdat{2};
                count = 1;
                selected = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'value');
                if cind ~= 1 %specific cluster
                    len = length(STUDY.cluster(cls(cind-1)).comps);
                    compid = cell(len+1,1);
                    compid{1} = 'All components';
                    % Convert from components numbering to the indexing form 'setXcomY'
                    for l = 1:len % go over the components of the cluster 
                        if ~isnan(STUDY.cluster(cls(cind-1)).sets(1,l))
                            subject = STUDY.datasetinfo(STUDY.cluster(cls(cind-1)).sets(1,l)).subject;
                            compid{l+1} = [  subject ' IC' num2str(STUDY.cluster(cls(cind-1)).comps(1,l)) ];
                        end
                    end
                    if isfield(STUDY.cluster, 'selected')
                        if ~isempty(STUDY.cluster(cls(cind-1)).selected)
                            selected = min(STUDY.cluster(cls(cind-1)).selected, 1+length(STUDY.cluster(cls(cind-1)).comps(1,:)));
                            STUDY.cluster(cls(cind-1)).selected = selected;
                        end
                    end

                else % All clusters accept 'Notclust' and 'Outliers'
                    count = 1;
                    for k = 1: length(cls)
                        if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) && ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) && ...
                                (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) 
                            for l = 1: length(STUDY.cluster(cls(k)).comps)
                                if ~isnan(STUDY.cluster(cls(k)).sets(1,l))
                                    subject = STUDY.datasetinfo(STUDY.cluster(cls(k)).sets(1,l)).subject;   % This line chokes on NaNs. TF 2007.05.31
                                    compid{count} = [ '''' STUDY.cluster(cls(k)).name ''' comp. ' ...
                                        num2str(l) ' (' subject  ' IC' num2str(STUDY.cluster(cls(k)).comps(l)) ')'];
                                    count = count +1;
                                end
                            end
                        end
                    end
                end
                if selected > length(compid), selected = 1; end
                set(findobj('parent', hdl, 'tag', 'clust_comp'), 'value', selected, 'String', compid);

            case 'plotsum'
                if clus ~= 1 % specific cluster option
                    [STUDY] = std_propplot(STUDY, ALLEEG, 'cluster', cls(clus-1));
                    % update Study history
                    a = ['STUDY = std_propplot(STUDY, ALLEEG, ''cluster'', '  num2str(cls(clus-1)) ' );'  ];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                else % all clusters
                    % All clusters does not include 'Notclust' and 'Outliers' clusters. 
                    tmpcls = [];
                    for k = 1:length(cls) 
                        if ~strncmpi(STUDY.cluster(cls(k)).name,'Notclust',8) && ~strncmpi(STUDY.cluster(cls(k)).name,'Outliers',8) && ...
                                (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) 
                            tmpcls = [tmpcls cls(k)];
                        end
                    end
                    [STUDY] = std_propplot(STUDY, ALLEEG, 'cluster', tmpcls);
                    % update Study history
                    a = ['STUDY = std_propplot(STUDY, ALLEEG, ''cluster'', ['  num2str(tmpcls) '] );'  ];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                end
               userdat{1}{2} = STUDY;
               set(hdl, 'userdat',userdat);    

            case 'plotcompsum'
                for ci = 1 : length(comp_ind)
                    % place holder for component properties % nima
                end

            case 'renameclust'
                STUDY.saved = 'no';
                clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                clus_num = get(findobj('parent', hdl, 'tag', 'clus_list'), 'Value') -1;
                if clus_num == 0  % 'all clusters' option 
                    return;
                end
                % Don't rename 'Notclust' and 'Outliers'  clusters.
                if strncmpi('Notclust',STUDY.cluster(cls(clus_num)).name,8) || strncmpi('Outliers',STUDY.cluster(cls(clus_num)).name,8) || ...
                        strncmpi('ParentCluster',STUDY.cluster(cls(clus_num)).name,13)
                    warndlg2('The ParentCluster, Outliers, and Notclust clusters cannot be renamed');
                    return;
                end
                old_name = STUDY.cluster(cls(clus_num)).name;
                rename_param  = inputgui( { [1] [1] [1]}, ...
                    { {'style' 'text' 'string' ['Rename ' old_name] 'FontWeight' 'Bold'} {'style' 'edit' 'string' '' 'tag' 'clus_rename' } {} }, ...
                '', 'Rename cluster - from pop_clustedit()' );
                if ~isempty(rename_param) %if not canceled
                    new_name = rename_param{1};
                    STUDY = std_renameclust(STUDY, ALLEEG, cls(clus_num), new_name);
                    % update Study history
                    a = ['STUDY = std_renameclust(STUDY, ALLEEG, ' num2str(cls(clus_num)) ', ' STUDY.cluster(cls(clus_num)).name  ');'];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                    % Renaming cluster in list
                    new_name = [ STUDY.cluster(cls(clus_num)).name ' (' num2str(length(STUDY.cluster(cls(clus_num)).comps))  ' ICs)'];
                    clus_name_list{clus_num+1} = renameclust( clus_name_list{clus_num+1}, new_name);
                    % Renaming Outlier cluster if exist
                    outlier_clust = std_findoutlierclust(STUDY,cls(clus_num));
                    if outlier_clust ~= 0
                        new_outliername = [ STUDY.cluster(cls(outlier_clust)).name ' (' num2str(length(STUDY.cluster(cls(outlier_clust)).comps))  ' ICs)'];
                        clus_name_list{outlier_clust+1} = renameclust( clus_name_list{outlier_clust+1}, new_outliername);
                    end
                    set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                    set(findobj('parent', hdl, 'tag', 'clus_rename'), 'String', '');
                    userdat{1}{2} = STUDY;
                    set(hdl, 'userdat',userdat); %update STUDY
                end

            case 'movecomp'
                STUDY.saved = 'no';
                old_clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value') -1;
                comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
                if old_clus == 0 % 'all clusters' option 
                    return;
                end
                % Don't reassign components of 'Notclust' or the 'ParentCluster'.
                if strncmpi('ParentCluster',STUDY.cluster(cls(old_clus)).name,13)  
                    warndlg2('Cannot reassign components of ''ParentCluster''.');
                    return;
                end
                old_name = STUDY.cluster(cls(old_clus)).name;
                ncomp = length(comp_ind); % number of selected components
                optionalcls =[];
                for k = 1:length(cls)
                    if (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))  && (k~= old_clus)
                        optionalcls = [optionalcls cls(k)];
                    end
                end                    
                reassign_param  = inputgui( { [1] [1] [1]}, ...
                    { {'style' 'text' 'string' strvcat(['Reassign ' fastif(ncomp >1, [num2str(length(comp_ind)) ' currently selected components'], ...
                                                                  'currently selected component') ], ...
                                [' from ' old_name ' to the cluster selected below']) 'FontWeight' 'Bold'} ...
                      {'style' 'listbox' 'string' {STUDY.cluster(optionalcls).name} 'tag' 'new_clus'} {} }, ...
                      '', 'Reassign cluster - from pop_clustedit()' ,[] , 'normal', [2 3 1] );
                if ~isempty(reassign_param) %if not canceled
                    new_clus = reassign_param{1};
                    comp_to_disp = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String');      
                    if strcmp(comp_to_disp{comp_ind(1)},'All components')
                        warndlg2('Cannot move all the components of the cluster - abort move components', 'Aborting move components');
                        return;
                    end
                    STUDY = std_movecomp(STUDY, ALLEEG,  cls(old_clus), optionalcls(new_clus), comp_ind - 1);                
                    % update Study history
                    a = ['STUDY = std_movecomp(STUDY, ALLEEG, ' num2str(cls(old_clus)) ', ' num2str(optionalcls(new_clus)) ', [' num2str(comp_ind - 1) ']);'];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
                    newind = find(cls == optionalcls(new_clus));

                    % update GUI
                    % ----------
                    clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');     
                    newname = [STUDY.cluster(optionalcls(new_clus)).name ' (' num2str(length(STUDY.cluster(optionalcls(new_clus)).comps))  ' ICs)'];
                    clus_name_list{newind+1} = renameclust( clus_name_list{newind+1}, newname);
                    newname = [STUDY.cluster(cls(old_clus)).name ' (' num2str(length(STUDY.cluster(cls(old_clus)).comps))  ' ICs)'];
                    clus_name_list{old_clus+1} = renameclust( clus_name_list{old_clus+1}, newname);
                    set( findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                    userdat{1}{2} = STUDY;
                    set(hdl, 'userdat',userdat); 
                    pop_clustedit('showclust',hdl);
                end          

            case 'moveoutlier'
                STUDY.saved = 'no';
                old_clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value') -1;
                comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value'); 
                if ~isempty(find(comp_ind ==1))
                    warndlg2('Cannot remove all the cluster components');
                    return;
                end
                if old_clus == 0 % 'all clusters' option 
                    return;
                end
                if strncmpi('Notclust',STUDY.cluster(cls(old_clus)).name,8) || strncmpi('ParentCluster',STUDY.cluster(cls(old_clus)).name,13)    % There are no outliers to 'Notclust'
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
                      '', 'Remove outliers - from pop_clustedit()' ,[] , 'normal', [1 3 1] );
                if ~isempty(reassign_param) %if not canceled
                    STUDY = std_moveoutlier(STUDY, ALLEEG,  cls(old_clus), comp_ind - 1);
                    clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                    outlier_clust = std_findoutlierclust(STUDY,cls(old_clus)); %find the outlier cluster for this cluster
                    oind = find(cls == outlier_clust); % the outlier clust index (if already exist) in the cluster list GUI
                    if ~isempty(oind) % the outlier clust is already presented in the cluster list GUI
                        newname = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                        clus_name_list{oind+1} = renameclust( clus_name_list{oind+1}, newname);
                    elseif outlier_clust == length(STUDY.cluster) % update the list with the Outlier cluster (if didn't exist before)
                        clus_name_list{end+1} = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                        userdat{2} = userdat{2} + 1; % update N, number of clusters in edit window 
                        cls(end +1) = length(STUDY.cluster); % update the GUI clusters list with the outlier cluster
                        userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                    end
                    newname = [STUDY.cluster(cls(old_clus)).name ' (' num2str(length(STUDY.cluster(cls(old_clus)).comps))  ' ICs)'];
                    clus_name_list{old_clus+1} = renameclust(clus_name_list{old_clus+1}, newname);
                    set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                    % update Study history
                    a = ['STUDY = std_moveoutlier(STUDY, ALLEEG, ' num2str(cls(old_clus)) ', [' num2str(comp_ind - 1) ']);'];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
                    userdat{1}{2} = STUDY;
                    set(hdl, 'userdat',userdat); 
                    pop_clustedit('showclust',hdl);    
                end

            case 'rejectoutliers'
                STUDY.saved = 'no';
                clus = get(findobj('parent', hdl, 'tag', 'clus_list'), 'Value') -1;
                if clus
                    std_name = STUDY.cluster(cls(clus)).name;
                    % Cannot reject outliers from 'Notclust', 'ParentCluster' and 'Outlier' clusters
                    if strncmpi('Notclust',std_name,8) || strncmpi('ParentCluster', std_name,13) || ...
                            strncmpi('Outliers',std_name,8)
                        warndlg2('Cannot reject outliers of ''Notclust'' or ''Outliers'' or ''ParentCluster'' clusters.');
                        return;
                    end
                    clusters = cls(clus);
                else
                    std_name = 'All clusters';
                    clusters = [];
                    for k = 1:length(cls)
                         if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) && ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) && ...
                                 ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)  
                            clusters = [ clusters cls(k)];
                        end
                    end
                end
                reject_param  = inputgui( { [1] [1] [4 1 2] [1]}, ...
                    { {'style' 'text' 'string' ['Reject "' std_name  '" outliers ' ] 'FontWeight' 'Bold'} {} ...
                       {'style' 'text' 'string' 'Move outlier components that are more than'} {'style' 'edit' 'string' '3' 'tag' 'outliers_std' } ...
                      {'style' 'text' 'string' 'standard deviations' } ...
                      {'style' 'text' 'string' [ 'from the "' std_name  '" centroid to an outlier cluster.']} }, ...
                      '', 'Reject outliers - from pop_clustedit()' );
                if ~isempty(reject_param) %if not canceled
                    ostd = reject_param{1}; % the requested outlier std
                    [STUDY] = std_rejectoutliers(STUDY, ALLEEG, clusters, str2num(ostd));  
                    % update Study history
                    a = ['STUDY = std_rejectoutliers(STUDY, ALLEEG, [ ' num2str(clusters) ' ], ' ostd ');'];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
                    clus_name_list = get(findobj('parent', hdl, 'tag', 'clus_list'), 'String');
                    for k = 1:length(clusters)
                        outlier_clust = std_findoutlierclust(STUDY,clusters(k)); %find the outlier cluster for this cluster
                        oind = find(cls == outlier_clust); % the outlier clust index (if already exist) in the cluster list GUI
                        if ~isempty(oind) % the outlier clust is already presented in the cluster list GUI
                            newname = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                            clus_name_list{oind+1} = renameclust( clus_name_list{oind+1}, newname);
                        else % update the list with the outlier cluster 
                            clus_name_list{end+1} = [STUDY.cluster(outlier_clust).name ' (' num2str(length(STUDY.cluster(outlier_clust).comps))  ' ICs)'];
                            userdat{2} = userdat{2} + 1; % update N, number of clusters in edit window 
                            cls(end +1) = outlier_clust; % update the GUI clusters list with the outlier cluster
                            userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                        end
                        clsind = find(cls == clusters(k));
                        newname = [STUDY.cluster(clusters(k)).name ' (' num2str(length(STUDY.cluster(clusters(k)).comps))  ' ICs)'];
                        clus_name_list{clsind+1} = renameclust( clus_name_list{clsind+1}, newname);
                        set(findobj('parent', hdl, 'tag', 'clus_list'), 'String', clus_name_list);
                    end
                    % If outlier cluster doesn't exist in the GUI window add it 
                    userdat{1}{2} = STUDY;
                    set(hdl, 'userdat',userdat); 
                    pop_clustedit('showclust',hdl);
                end

            case 'createclust'
                STUDY.saved = 'no';
                create_param  = inputgui( { [1] [1 1] [1]}, ...
                    { {'style' 'text' 'string' 'Create new empty cluster' 'FontWeight' 'Bold'} ...
                       {'style' 'text' 'string' 'Enter cluster name:'} {'style' 'edit' 'string' '' } {} }, ...
                      '', 'Create new empty cluster - from pop_clustedit()' );
                if ~isempty(create_param) %if not canceled
                    clus_name = create_param{1}; % the name of the new cluster
                    [STUDY] = std_createclust(STUDY, ALLEEG, 'name', clus_name); 
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
                        a = ['STUDY = std_createclust(STUDY, ALLEEG);'];
                    else
                        a = ['STUDY = std_createclust(STUDY, ALLEEG, ''name'', ' clus_name ');'];
                    end                
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                    userdat{1}{2} = STUDY;
                    userdat{2} = userdat{2} + 1; % update N, the number of cluster options in edit window 
                    cls(end +1) = length(STUDY.cluster); % update the GUI clusters list with the new cluster
                    userdat{1}{3} = cls;  % update cls, the cluster indices in edit window
                    set(hdl, 'userdat',userdat); %update STUDY, cls and N
                end

            case 'mergeclusters'
                STUDY.saved = 'no';
                clus_names = get(findobj('parent', hdl, 'tag', 'clus_list'), 'string') ;
                optionalcls =[];
                for k = 2:length(clus_names)
                    if (~strncmpi('Notclust',clus_names{k},8)) && (~strncmpi('Outliers',clus_names{k},8)) && ...
                            (~strncmpi('ParentCluster',clus_names{k},13))
                        optionalcls = [optionalcls k];
                    end
                end           
                reassign_param  = inputgui( { [1] [1] [1] [2 1] [1]}, ...
                    { {'style' 'text' 'string' 'Select clusters to Merge' 'FontWeight' 'Bold'} ...
                      {'style' 'listbox' 'string' clus_names(optionalcls) 'tag' 'new_clus' 'max' 3 'min' 1} {} ...
                      {'style' 'text' 'string' 'Optional, enter a name for the merged cluster:' 'FontWeight' 'Bold'} ...
                      {'style' 'edit' 'string' ''} {} }, ...
                      '', 'Merge clusters - from pop_clustedit()' ,[] , 'normal', [1 3 1 1 1] );
                  if ~isempty(reassign_param)
                      std_mrg = cls(optionalcls(reassign_param{1})-1);
                      name = reassign_param{2};
                      allleaves = 1;
                      N = userdat{2};
                      for k = 1: N %check if all leaves
                          if ~isempty(STUDY.cluster(cls(k)).child) 
                              allleaves = 0;
                          end
                      end                     
                      [STUDY] = std_mergeclust(STUDY, ALLEEG, std_mrg, name); 
                      % 
                      % update Study history
                      % 
                      if isempty(name)
                          a = ['STUDY = std_mergeclust(STUDY, ALLEEG, [' num2str(std_mrg) ']);'];
                      else
                          a = ['STUDY = std_mergeclust(STUDY, ALLEEG, [' num2str(std_mrg) '], ' name ');'];
                      end                  
                      STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
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
                        cls = setdiff_bc(cls, std_mrg); % remove from the GUI clusters list the merged clusters
                        cls(end+1) = length(STUDY.cluster); % update the GUI clusters list with the new cluster
                        N  = length(cls);
                        %
                        % update the first option on the GUI list : 'All 10 cluster centroids'
                        % with the new number of cluster centroids
                        %
                        ti = strfind(clus_names{1},'cluster'); %get the number of clusters centroid 
                        cent = num2str(str2num(clus_names{1}(5:ti-2))+1- length(std_mrg)); % new number of centroids
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
    catch
        eeglab_error;
    end       
end

function newname = renameclust(oldname, newname);
    
    tmpname = deblank(oldname(end:-1:1));
    strpos  = strfind(oldname, tmpname(end:-1:1));
    
    newname = [ oldname(1:strpos-1) newname ];

    
