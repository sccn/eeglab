% pop_preclust() - prepare STUDY components' location and activity measures for later clustering.
%                  Collect information in an interactive pop-up query window. To pre-cluster
%                  from the commandline, use eeg_preclust(). After data entry into the pop window,
%                  Selected measures (one or more from options: ERP, dipole locations, spectra,
%                  scalp maps, ERSP, and ITC) are computed for each dataset in the STUDY 
%                  set, unless they already present. After all requested measures are computed 
%                  and saved in the STUDY datasets, a PCA  matrix (by runica() with 'pca' option) 
%                  is constructed (this is the feature reduction step). This matrix will be used 
%                  as input to the clustering  algorithm in pop_clust(). pop_preclust() allows 
%                  selection of a subset of components to use in the clustering. This subset 
%                  may be a user-specified component subset, components with dipole model residual 
%                  variance lower than a defined threshold (see dipfit()), or components from 
%                  an already existing cluster (for hierarchical clustering). The EEG datasets
%                  in the ALLEEG structure are updated, and updated EEG sets are saved to disk.
%                  Calls eeg_preclust().
% Usage:    
%                >> [STUDY, ALLEEG] = pop_preclust(STUDY, ALLEEG); % pop up interactive window
%                >> [STUDY, ALLEEG] = pop_preclust(STUDY, ALLEEG, clustind); % sub-cluster 
%
% Graphic interface:
% Inputs:
%   STUDY        - an EEGLAB STUDY set (containing loaded EEG structures)
%   ALLEEG       - ALLEEG data structure, can also be an EEG dataset structure.
%   clustind     - a (single) cluster index for sub-clustering (hierarchical clustering) --
%                  for example to cluster a mu component cluster into left mu and right mu 
%                  sub-clusters. Should be empty for top level (whole STUDY) clustering 
%                  {default: []}
% Outputs:
%   STUDY        - the input STUDY set with added pre-clustering data, for use by pop_clust() 
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures modified by adding preprocessing 
%                  data (pointers to float files that hold ERSP, spectrum, etc. information).
%
% Authors: Arnaud Delorme, Hilit Serby & Scott Makeig, SCCN, INC, UCSD, May 13, 2004-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, May 13,2004, hilit@sccn.ucsd.edu
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

% $pop_preclust.m,v$

function [STUDY, ALLEEG, com] = pop_preclust(varargin)

com = '';
if ~isstr(varargin{1}) %intial settings
    if length(varargin) < 2
        error('pop_preclust(): needs both ALLEEG and STUDY structures');
    end
    STUDY  = varargin{1};
    ALLEEG = varargin{2};
    if length(varargin) >= 3
        if length(varargin{3}) > 1
            error('pop_preclust(): To cluster components from several clusters, merge them first!');
        end
        cluster_ind = varargin{3};
    else
        cluster_ind = [];
    end
    
    % Create default ERSP / ITC time/freq. paramters 
    ersp.f = '3 50';
    ersp.c = '3 0.5';
    ersp.p = '4';
    ersp.a = '0.01';
    seti = STUDY.datasetinfo(1).index; %first dataset in ALLEEG that is part of STUDY
    [time_range, winsize] = compute_ERSP_times(str2num(ersp.c),  ALLEEG(seti).srate, ...
        [ALLEEG(seti).xmin ALLEEG(seti).xmax]*1000, str2num(ersp.f(1)),str2num(ersp.p)); 
    ersp.t = [num2str(round(time_range(1))) ' ' num2str(round(time_range(2))) ];
    erspparams_str = ['                                                             ''frange'', [' ...
            ersp.f '], ''cycles'', [' ersp.c '], ''alpha'', ' ersp.a ', ''padratio'', ' ersp.p ', ''tlimits'', [' ersp.t ']'];

    
    scalp_options = {'Use channel values' 'Use Laplacian values' 'Use Gradient values'} ;
    
    % cluster text
    % ------------
    % load leaf clusters
    cls = [];
    for k = 1:length(STUDY.cluster)
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
        end
    end
    N = length(cls); %number of clusters
    show_options = {};
    num_cls = 0;
    for k = 2:N
        show_options{k-1} = [STUDY.cluster(cls(k)).name ' (' num2str(length(STUDY.cluster(cls(k)).comps))  ' ICs)'];
        if (~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8)) & (~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8))  & ...
                (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
            num_cls = num_cls + 1;
        end
    end

    % callbacks
    % ---------
    show_clust      = [ 'pop_preclust(''showclust'',gcf);'];
    show_comps      = [ 'pop_preclust(''showcomplist'',gcf);'];
    help_spectopo =  ['pophelp(''spectopo'')'];         
	set_spectra = ['pop_preclust(''setspec'',gcf);']; 
    set_erp = ['pop_preclust(''seterp'',gcf);']; 
    set_scalp = ['pop_preclust(''setscalp'',gcf);']; 
    set_dipole = ['pop_preclust(''setdipole'',gcf);'];
    set_ersp = ['pop_preclust(''setersp'',gcf);']; 
    set_itc = ['pop_preclust(''setitc'',gcf);']; 
    help_clusteron = ['pophelp(''cls_helpselecton'');']; 
    help_ersp = ['pophelp(''pop_timef'')'];
    preclust_PCA = ['pop_preclust(''preclustOK'',gcf);'];           
    ersp_params = ['pop_preclust(''erspparams'',gcf);']; 
    ersp_edit =  ['pop_preclust(''erspedit'',gcf);']; 
    saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
    browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_preclust()''); ' ... 
                  'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    str_name   = ['Precompute measures to cluster components from study ''' STUDY.name '''' ];
    str_time   = [ num2str(round(ALLEEG(seti).xmin*1000)) '  ' num2str(round(ALLEEG(seti).xmax*1000)) ];
    
    gui_spec = { ...
    {'style' 'text'       'string' str_name 'FontSize' 12 'FontWeight' 'Bold' 'horizontalalignment' 'left'} ...
	{'style' 'text'       'string' 'Select cluster leaves from which to refine clustering (remove cluster info (edit study) to restart from scratch)' } {} ...
    {'style' 'listbox'    'string' show_options 'value' 1 'tag' 'clus_list' 'Callback' show_clust 'max' 2 } {}  {} ...
    {'style' 'text'       'string' 'Pre-compute or load           Dim.       Norm.  Rel. weight' 'FontWeight' 'Bold'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'spectra_on' 'value' 0 'Callback' set_spectra 'userdata' '1'}  ...
	{'style' 'text'       'string' 'spectra' 'FontWeight' 'Bold' 'horizontalalignment' 'center' } ...
	{'style' 'edit'       'string' '10' 'tag' 'spectra_PCA' 'enable' 'off' 'userdata' 'specP'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'spectra_norm' 'value' 1 'enable' 'off' 'userdata' 'specP' } ...
    {'style' 'edit'       'string' '1' 'tag' 'spectra_weight' 'enable' 'off' 'userdata' 'specP'} ...
    {'style' 'text'       'string' 'Frequency range [Hz]' 'tag' 'spectra_freq_txt' 'userdata' 'spec' 'enable' 'off' } ...
	{'style' 'edit'       'string' '3  30'  'tag' 'spectra_freq_edit' 'userdata' 'spec' 'enable' 'off' } ...
    {'style' 'checkbox'   'string' '' 'tag' 'erp_on' 'value' 0 'Callback' set_erp 'userdata' '1'}  ...
	{'style' 'text'       'string' 'ERPs' 'FontWeight' 'Bold' 'horizontalalignment' 'center' } ...
    {'style' 'edit'       'string' '10' 'tag' 'erp_PCA' 'enable' 'off' 'userdata' 'erpP'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'erp_norm' 'value' 1 'enable' 'off' 'userdata' 'erpP' } ...
    {'style' 'edit'       'string' '1' 'tag' 'erp_weight' 'enable' 'off' 'userdata' 'erpP'} ...
    {'style' 'text'       'string' 'Latency window [ms]' 'tag' 'erp_time_txt' 'userdata' 'erp' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_time 'tag' 'erp_time_edit' 'userdata' 'erp' 'enable' 'off' } ...
	{'style' 'checkbox'   'string' '' 'tag' 'dipole_on' 'value' 0 'Callback' set_dipole 'userdata' '1'} ...
	{'style' 'text'       'string' 'dipoles' 'FontWeight' 'Bold' 'HorizontalAlignment' 'center' } ...
	{'style' 'text'       'string' '3' 'enable' 'off' 'userdata' 'dipoleP' } ...
	{'style' 'checkbox'   'string' '' 'tag' 'locations_norm' 'value' 1 'enable' 'off' 'userdata' 'dipoleP'}  ...
	{'style' 'edit'       'string' '10' 'tag' 'locations_weight' 'enable' 'off' 'userdata' 'dipoleP'} {} {} ...
    {'style' 'checkbox'   'string' '' 'tag' 'scalp_on' 'value' 0 'Callback' set_scalp 'userdata' '1'} ...
	{'style' 'text'       'string' 'scalp maps' 'FontWeight' 'Bold' 'HorizontalAlignment' 'center' } ...
	{'style' 'edit'       'string' '10' 'tag' 'scalp_PCA' 'enable' 'off' 'userdata' 'scalpP'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'scalp_norm' 'value' 1 'enable' 'off' 'userdata' 'scalpP'}   ...
    {'style' 'edit'       'string' '1' 'tag' 'scalp_weight' 'enable' 'off' 'userdata' 'scalpP'} ...
    {'style' 'popupmenu'  'string' scalp_options 'value' 1 'tag' 'scalp_choice' 'enable' 'off' 'userdata' 'scalp' } ...
    {'style' 'checkbox'   'string' 'Absolute values' 'value' 1 'tag'  'scalp_absolute' 'enable' 'off' 'userdata' 'scalp' } ...
    {'style' 'checkbox'   'string' '' 'tag' 'ersp_on' 'value' 0 'Callback' set_ersp 'userdata' '1'} ...
	{'style' 'text'       'string' 'ERSPs' 'FontWeight' 'Bold' 'horizontalalignment' 'center' } ...
    {'style' 'edit'       'string' '10' 'tag' 'ersp_PCA' 'enable' 'off' 'userdata' 'erspP'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'ersp_norm' 'value' 1 'enable' 'off' 'userdata' 'erspP'}   ...
	{'style' 'edit'       'string' '1' 'tag' 'ersp_weight' 'enable' 'off' 'userdata' 'erspP'} ...
    {'style' 'pushbutton' 'string' 'Time/freq. parameters' 'tag' 'ersp_push' 'value' 1 'enable' 'off' 'userdata' 'ersp' 'Callback' ersp_params} ...
    {'style' 'edit'       'string' erspparams_str 'tag' 'ersp_params' 'enable' 'off' 'userdata' 'ersp' 'Callback' ersp_edit}...
    {'style' 'checkbox'   'string' '' 'tag' 'itc_on' 'value' 0 'Callback' set_itc 'userdata' '1'} ...
	{'style' 'text'       'string' 'ITCs'  'FontWeight' 'Bold' 'horizontalalignment' 'center' } ...
    {'style' 'edit'       'string' '10' 'tag' 'itc_PCA' 'enable' 'off' 'userdata' 'itcP'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'itc_norm' 'value' 1 'enable' 'off' 'userdata' 'itcP'}   ...
	{'style' 'edit'       'string' '1' 'tag' 'itc_weight' 'enable' 'off' 'userdata' 'itcP'} ...
    {'style' 'pushbutton' 'string' 'Time/freq. parameters' 'tag' 'itc_push' 'value' 1 'enable' 'off' 'userdata' 'itc' 'Callback' ersp_params} ...
    {'style' 'edit'       'string' erspparams_str 'tag' 'itc_params' 'enable' 'off' 'userdata' 'itc' 'Callback' ersp_edit} {} ...
    {'style' 'checkbox'   'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} ...
	{'style' 'text'       'string' 'Save STUDY set to disk' 'FontWeight' 'Bold'  } ...
    {'style' 'edit'       'string' fullfile(STUDY.filepath, STUDY.filename) 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
	{'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} };
  
	%{'style' 'checkbox'   'string' '' 'tag' 'preclust_PCA'  'Callback' preclust_PCA 'value' 0} ...
	%{'style' 'text'       'string' 'Do not prepare data matrix for clustering at this time.' 'FontWeight' 'Bold'  } {} ...

    fig_arg{1} = { ALLEEG STUDY cls };
    fig_arg{2} = ersp;
    geomline = [0.45 2 1 0.45 1 3 3 ];
    geometry = { [1] [1] [1 1 1] [1] ...
                 [3] geomline geomline geomline geomline geomline geomline [1] [0.32 2 5 0.8] };
    geomvert = [ 1 1 4 1 1 1 1 1 1 1 1 1 1 1 1 1 ];

    if isempty(show_options)
        gui_spec(2:6) = { {} ...
            { 'style' 'text'      'string' [ 'Among the pre-selected components (Edit study),' ...
                            'remove those which dipole res. var, exceed' ] 'tag' 'dipole_select_on' }  ...
            {'style' 'edit'       'string' '0.15' 'horizontalalignment' 'center' 'tag' 'dipole_rv'} ...
            {'style' 'text'       'string'  '(empty=all)'} {} };
        geometry{3} = [2.5 0.25 0.4];
        geomvert(3) = 1;
    end;
    
	[preclust_param, userdat2, strhalt, os] = inputgui( 'geometry', geometry, 'uilist', gui_spec, 'geomvert', geomvert, ...
                                                      'helpcom', ' pophelp(''eeg_preclust'')', ...
                                                      'title', 'Select and compute component measures for later clustering -- pop_preclust()', ...
                                                      'userdata', fig_arg);	
	if isempty(preclust_param), return; end;
    
    options = { STUDY, ALLEEG };
    
    % precluster on what?
    % -------------------
    if ~isempty(show_options)
        options{3} = os.clus_list+2; % hierarchical clustering
        options{4} = [];
    else
        options{3} = [];
        options{4} = [];
        if ~isempty(str2num(os.dipole_rv)) %dipole information is used for component selection
            options{end+1} = { 'dipselect' 'rv' os.dipole_rv };
        end
    end;

    %if ~(os.preclust_PCA) %create PCA data for clustering
    %preclust_command = '[STUDY ALLEEG] = eeg_createdata(STUDY, ALLEEG, ';
    %end
    
    % Spectrum option is on
    % --------------------
    if os.spectra_on== 1 
        options{end+1} = {  'spec' 'npca' str2num(os.spectra_PCA) 'norm' os.spectra_norm ...
                            'weight' str2num(os.spectra_weight)  'freqrange' os.spectra_freq_edit };
    end
    
    % ERP option is on
    % ----------------
    if os.erp_on == 1 
        options{end} = { 'erp' 'npca' str2num(os.erp_PCA) 'norm' os.erp_norm ...
                         'weight' str2num(os.erp_weight) 'timewindow' eval( [ '[' os.erp_time_edit ']' ]) };
    end
    
    % Scalp maps option is on
    % ----------------------
    if os.scalp_on == 1 
        if os.scalp_absolute %absolute maps
            abso = 1;
        else abso = 0;
        end
        if (os.scalp_choice == 2)  %Laplacian scalp maps
            options{end+1} = { 'scalpLaplac' 'npca' str2num(os.scalp_PCA) 'norm' os.scalp_norm ...
                               'weight' str2num(os.scalp_weight) 'abso' abso };
        elseif (os.scalp_choice == 3)  %Gradient scalp maps
            options{end+1} = { 'scalpGrad' 'npca' str2num(os.scalp_PCA) 'norm' os.scalp_norm, ...
                               'weight' str2num(os.scalp_weight) 'abso' abso };
        elseif (os.scalp_choice == 1) %scalp map case
            options{end+1} = { 'scalp' 'npca' str2num(os.scalp_PCA) 'norm' os.scalp_norm, ...
                               'weight' str2num(os.scalp_weight) 'abso' abso };
        end
    end
    
    % Dipole option is on
    % -------------------
    if os.dipole_on == 1 
        options{end+1} = { 'dipoles' 'norm' os.locations_norm 'weight' str2num(os.locations_weight) };
    end
    
    % ERSP option is on
    % -----------------
    if os.ersp_on  == 1 
        ersp = userdat2{2};
        options{end+1} = { 'ersp' 'npca' str2num(os.ersp_PCA) 'freqrange' num2str(ersp.f) ...
                           'cycles' num2str(ersp.c) 'alpha' num2str(ersp.a) 'padratio' num2str(ersp.p) ...
                           'timewindow' num2str(ersp.t) 'norm' os.ersp_norm ...
                           'weight' str2num(os.ersp_weight) };
    end
    
    % ITC option is on 
    % ----------------
    if os.itc_on  == 1 
        ersp = userdat2{2};
        options{end+1} = { 'itc' 'npca' str2num(os.ersp_PCA) 'freqrange' num2str(ersp.f) ...
                           'cycles' num2str(ersp.c) 'alpha' num2str(ersp.a) 'padratio' num2str(ersp.p) ...
                           'timewindow' num2str(ersp.t) 'norm' os.ersp_norm ...
                           'weight' str2num(os.ersp_weight) };
    end       
    
    % evaluate command
    % ----------------
    [STUDY ALLEEG] = eeg_preclust(STUDY, ALLEEG, options{:});
    com = sprintf('%s\n[STUDY ALLEEG] = eeg_preclust(STUDY, ALLEEG, %s);',  STUDY.history, vararg2str(options));
    
    % save updated STUDY to the disk
    % ------------------------------
    if os.saveSTUDY == 1 
        if ~isempty(os.studyfile)
            [filepath filename ext] = fileparts(os.studyfile);
            STUDY.filename = [ filename ext ];
            STUDY.filepath = filepath;
        end;
        STUDY = pop_savestudy(STUDY, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
        com = sprintf('%s\nSTUDY = pop_savestudy(STUDY, %s);',  com, ...
                      vararg2str( { 'filename', STUDY.filename, 'filepath', STUDY.filepath }));
    end
else
    hdl = varargin{2}; %figure handle
    userdat = get(varargin{2}, 'userdat');    
    ALLEEG  = userdat{1}{1};
    STUDY   = userdat{1}{2};
    cls     = userdat{1}{3};
    N       = length(cls);

    switch  varargin{1}
               
        case 'setspec'
            set_spec =  get(findobj('parent', hdl, 'tag', 'spectra_on'), 'value'); 
            set(findobj('parent', hdl, 'userdata', 'spec'), 'enable', fastif(set_spec,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
                set(findobj('parent', hdl, 'userdata', 'specP'), 'enable', 'off');
            else
                set(findobj('parent', hdl, 'userdata', 'specP'), 'enable', fastif(set_spec,'on','off'));
            end
        case 'seterp'
            set_erp =  get(findobj('parent', hdl, 'tag', 'erp_on'), 'value'); 
            set(findobj('parent', hdl, 'userdata', 'erp'), 'enable', fastif(set_erp,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
                set(findobj('parent', hdl, 'userdata', 'erpP'), 'enable', 'off');
            else
                set(findobj('parent', hdl, 'userdata', 'erpP'), 'enable', fastif(set_erp,'on','off'));
            end
        case 'setscalp'
            set_scalp =  get(findobj('parent', hdl, 'tag', 'scalp_on'), 'value'); 
            set(findobj('parent', hdl, 'userdata', 'scalp'), 'enable', fastif(set_scalp,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
                set(findobj('parent', hdl, 'userdata', 'scalpP'), 'enable', 'off');
            else
                set(findobj('parent', hdl, 'userdata', 'scalpP'), 'enable', fastif(set_scalp,'on','off'));
            end
        case 'setdipole'
            set_dipole =  get(findobj('parent', hdl, 'tag', 'dipole_on'), 'value'); 
            set(findobj('parent', hdl, 'userdata', 'dipole'), 'enable', fastif(set_dipole,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
               set(findobj('parent', hdl, 'userdata', 'dipoleP'), 'enable','off');
           else
               set(findobj('parent', hdl, 'userdata', 'dipoleP'), 'enable', fastif(set_dipole,'on','off'));
           end
        case 'setersp'
            set_ersp =  get(findobj('parent', hdl, 'tag', 'ersp_on'), 'value'); 
            set(findobj('parent', hdl,'userdata', 'ersp'), 'enable', fastif(set_ersp,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
                set(findobj('parent', hdl,'userdata', 'erspP'), 'enable', 'off');
            else
                set(findobj('parent', hdl,'userdata', 'erspP'), 'enable', fastif(set_ersp,'on','off'));
            end
            set_itc =  get(findobj('parent', hdl, 'tag', 'itc_on'), 'value'); 
            set(findobj('parent', hdl,'tag', 'ersp_push'), 'enable', fastif(set_itc,'off','on'));
            set(findobj('parent', hdl,'tag', 'ersp_params'), 'enable', fastif(set_itc,'off','on'));
             if  (set_itc & (~set_ersp) )
                set(findobj('parent', hdl,'tag', 'itc_push'), 'enable', 'on');
                set(findobj('parent', hdl,'tag', 'itc_params'), 'enable', 'on');
            end
       case 'setitc'
            set_itc =  get(findobj('parent', hdl, 'tag', 'itc_on'), 'value'); 
            set(findobj('parent', hdl,'userdata', 'itc'), 'enable', fastif(set_itc,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
                set(findobj('parent', hdl,'userdata', 'itcP'), 'enable','off');
            else
                set(findobj('parent', hdl,'userdata', 'itcP'), 'enable', fastif(set_itc,'on','off'));
            end
            set_ersp = get(findobj('parent', hdl, 'tag', 'ersp_on'), 'value'); 
            set(findobj('parent', hdl,'tag', 'itc_push'), 'enable', fastif(set_ersp,'off','on'));
            set(findobj('parent', hdl,'tag', 'itc_params'), 'enable', fastif(set_ersp,'off','on'));
            if  (set_ersp & (~set_itc) )
                set(findobj('parent', hdl,'tag', 'ersp_push'), 'enable', 'on');
                set(findobj('parent', hdl,'tag', 'ersp_params'), 'enable', 'on');
            end
        case 'erspedit'
            ersp_params = get(findobj('parent', hdl, 'tag', 'ersp_params'), 'string');
            if ~isempty( ersp_params)
                ersp_str = [ 'ersp_p = struct( ' ersp_params ');' ];
                eval(ersp_str);
                seti = STUDY.datasetinfo(1).index; %first dataset in ALLEEG that is part of STUDY
                [time_range, winsize] = compute_ERSP_times(ersp_p.cycles,  ALLEEG(seti).srate,...
                    [ALLEEG(seti).xmin ALLEEG(seti).xmax]*1000, ersp_p.frange(1),ersp_p.padratio); 
                if isfield(ersp_p, 'tlimits')
                    if ersp_p.tlimits(1) < time_range(1)
                        ersp_p.tlimits(1) = time_range(1);
                    end
                    if ersp_p.tlimits(2) > time_range(2)
                        ersp_p.tlimits(2) = time_range(2);
                    end
                else
                    ersp_p.tlimits = time_range;
                end
                ersp.c = ersp_p.cycles;
                ersp.a = ersp_p.alpha;
                ersp.p = ersp_p.padratio;
                ersp.t = ersp_p.tlimits;
                ersp.f = ersp_p.frange;
                userdat{2} = ersp;
                set(findobj('parent', hdl, 'tag', 'ersp_params'), 'string', ...
                    ['                                                             ''frange'', [' num2str([ersp.f]) '], ''cycles'', [' ...
                        num2str([ersp.c]) '], ''alpha'', ' num2str(ersp.a) ', ''padratio'', ' num2str(ersp.p) ', ''tlimits'', [' num2str([ersp.t]) ']']);
                set(findobj('parent', hdl, 'tag', 'itc_params'), 'string', ...
                    ['                                                             ''frange'', [' num2str([ersp.f]) '], ''cycles'', [' ...
                    num2str([ersp.c]) '], ''alpha'', ' num2str(ersp.a) ', ''padratio'', ' num2str(ersp.p) ', ''tlimits'', [' num2str([ersp.t]) ']']);                       
                set(hdl, 'userdat',userdat); 
            end
        case 'erspparams'
            ersp = userdat{2};
            ERSP_timewindow =  ['get_ersptime(ALLEEG, STUDY, gcf);']; 
            [ersp_paramsout, erspuserdat, strhalt, erspstruct] = inputgui( { [1] [3 1] [3 1] [3 1] [3 1] [3 1] [1]}, ...
                    { {'style' 'text' 'string' 'ERSP and ITC time/freq. parameters' 'FontWeight' 'Bold'} ...
                    {'style' 'text' 'string' 'Frequency range [Hz]' 'tag' 'ersp_freq' } ...
                    {'style' 'edit' 'string' ersp.f 'tag' 'ersp_f' 'Callback' ERSP_timewindow } ...
                    {'style' 'text' 'string' 'Wavelet cycles (see help timef())' 'tag' 'ersp_cycle' } ...
                    {'style' 'edit' 'string' ersp.c 'tag' 'ersp_c' 'Callback' ERSP_timewindow} ...    
                    {'style' 'text' 'string' 'Significance level (< 0.1)' 'tag' 'ersp_alpha' } ...
                    {'style' 'edit' 'string' ersp.a 'tag' 'ersp_a'} ...
                    {'style' 'text' 'string' 'timef() padratio' 'tag' 'ersp_pad' } ...
                    {'style' 'edit' 'string' ersp.p 'tag' 'ersp_p' 'Callback'  ERSP_timewindow} ...
					{'style' 'text' 'string' 'Desired time window within the indicated latency range [ms]' 'tag' 'ersp_trtxt' } ...
					{'style' 'edit' 'string' ersp.t 'tag' 'ersp_timewindow' 'Callback'  ERSP_timewindow} {} }, ...
 	                'pophelp(''pop_timef'')', 'Select clustering ERSP and ITC time/freq. parameters -- pop_preclust()');    
            if ~isempty(ersp_paramsout)
                ersp.f = erspstruct(1).ersp_f;
                ersp.c = erspstruct(1).ersp_c;
                ersp.p = erspstruct(1).ersp_p;
                ersp.a = erspstruct(1).ersp_a;
                ersp.t = erspstruct(1).ersp_timewindow;
                userdat{2} = ersp;
                set(findobj('parent', hdl, 'tag', 'ersp_params'), 'string', ...
                    ['                                                             ''frange'', [' ersp.f '], ''cycles'', [' ...
                     ersp.c '], ''alpha'', ' ersp.a ', ''padratio'', ' ersp.p ', '' tlimits'', [' ersp.t ']']);
                set(findobj('parent', hdl, 'tag', 'itc_params'), 'string', ...
                    ['                                                             ''frange'', [' ersp.f '], ''cycles'', [' ...
                     ersp.c '], ''alpha'', ' ersp.a ', ''padratio'', ' ersp.p ', '' tlimits'', [' ersp.t ']']);
                set(hdl, 'userdat',userdat); 
            end
       case 'preclustOK'
           set_PCA =  get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value'); 
           set_ersp =  get(findobj('parent', hdl, 'tag', 'ersp_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'erspP'), 'enable', fastif(~set_PCA & set_ersp,'on','off'));
           set_itc =  get(findobj('parent', hdl, 'tag', 'itc_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'itcP'), 'enable', fastif(~set_PCA & set_itc,'on','off'));
           set_erp =  get(findobj('parent', hdl, 'tag', 'erp_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'erpP'), 'enable', fastif(~set_PCA & set_erp,'on','off'));
           set_spec =  get(findobj('parent', hdl, 'tag', 'spectra_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'specP'), 'enable', fastif(~set_PCA & set_spec,'on','off'));
           set_scalp =  get(findobj('parent', hdl, 'tag', 'scalp_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'scalpP'), 'enable', fastif(~set_PCA & set_scalp,'on','off'));
           set_dipole =  get(findobj('parent', hdl, 'tag', 'dipole_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'dipoleP'), 'enable', fastif(~set_PCA & set_dipole,'on','off'));
           set(findobj('parent', hdl,'tag', 'chosen_component'), 'enable', fastif(~set_PCA,'on','off'));
           set(findobj('parent', hdl,'tag', 'dipole_rv'), 'enable', fastif(~set_PCA,'on','off'));
           set(findobj('parent', hdl,'tag', 'compcls_str'), 'enable', fastif(~set_PCA,'on','off'));
    end
end
STUDY.saved = 'no';
