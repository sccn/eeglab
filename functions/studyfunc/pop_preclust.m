% pop_preclust() - prepare STUDY components' location and activity measures for later clustering.
%                  Collect information in an interactive pop-up query window. To pre-cluster
%                  from the commandline, use std_preclust(). After data entry into the pop window,
%                  selected measures (one or more from options: ERP, dipole locations, spectra,
%                  scalp maps, ERSP, and ITC) are computed for each dataset in the STUDY 
%                  set, unless they already present. After all requested measures are computed 
%                  and saved in the STUDY datasets, a PCA  matrix (by runica() with 'pca' option) 
%                  is constructed (this is the feature reduction step). This matrix will be used 
%                  as input to the clustering  algorithm in pop_clust(). pop_preclust() allows 
%                  selection of a subset of components to cluster. This subset may either be 
%                  user-specified, all components with dipole model residual variance lower than 
%                  a defined threshold (see dipfit()), or components from an already existing cluster 
%                  (for hierarchical clustering). The EEG datasets in the ALLEEG structure are 
%                  updated; then the updated EEG sets are saved to disk.  Calls std_preclust().
% Usage:    
%                >> [STUDY, ALLEEG] = pop_preclust(STUDY, ALLEEG); % pop up interactive window
%                >> [STUDY, ALLEEG] = pop_preclust(STUDY, ALLEEG, clustind); % sub-cluster 
%
% Inputs:
%   STUDY        - STUDY set structure containing (loaded) EEG dataset structures
%   ALLEEG       - ALLEEG vector of EEG structures, else a single EEG dataset.
%   clustind     - (single) cluster index to sub-cluster, Hhierarchical clustering may be
%                  useful, for example, to separate a bilteral cluster into left and right 
%                  hemisphere sub-clusters. Should be empty for whole STUDY (top level) clustering 
%                  {default: []}
% Outputs:
%   STUDY        - the input STUDY set with added pre-clustering data for use by pop_clust() 
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures modified by adding 
%                  pre-clustering data (pointers to .mat files that hold cluster measure information).
%
% Authors: Arnaud Delorme, Hilit Serby & Scott Makeig, SCCN, INC, UCSD, May 13, 2004-
%
% See also: std_preclust()

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, May 13,2004, hilit@sccn.ucsd.edu
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

function [STUDY, ALLEEG, com] = pop_preclust(varargin)

com = '';

if ~ischar(varargin{1}) %initial settings
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
    
    scalp_options = {'Use channel values' 'Use Laplacian values' 'Use Gradient values'} ;
    
    if isempty(ALLEEG)
        error('STUDY contains no datasets');
    end
         
    if any(isnan(STUDY.cluster(1).sets(:)))
        warndlg2( [ 'FOR CLUSTERING, YOU MAY ONLY USE DIPOLE OR SCALP MAP CLUSTERING.' 10 ...
                    'This is because some datasets do not have ICA pairs. Look for NaN values in ' 10 ...
                    'STUDY.cluster(1).sets which indicate missing datasets. Each column in this ' 10 ...
                    'array indicate datasets with common ICA decompositions' ]);
    end
    if length(STUDY.design(STUDY.currentdesign).cases.value) ~= length(STUDY.subject)
        warndlg2( [ 'GO BACK TO THE DESIGN INTERFACE AND SELECT A DESIGN THAT ' 10 ...
                    'INCLUDES ALL DATASETS. Some subjects or datasets have been excluded' 10 ...
                    'in the current design. ICA clusters are common to all designs, so all' 10 ...
                    'all datasets must be included for clustering. After clustering, you' 10 ...
                    'will then be able to select a different design (and keep the clustered' 10 ...
                    'components) if you wish to exclude a subject or group of dataset.' ]);
        return;
    end
    
    % cluster text
    % ------------
    % load leaf clusters
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

    % callbacks
    % ---------
    erspparams_str = [ '''cycles'', [3 0.5], ''padratio'', 1' ];
    specparams_str = '';
    show_clust      = [ 'pop_preclust(''showclust'',gcf);'];
    show_comps      = [ 'pop_preclust(''showcomplist'',gcf);'];
    help_spectopo =  ['pophelp(''spectopo'')'];         
	set_spectra  = ['pop_preclust(''setspec'',gcf);']; 
    set_erp      = ['pop_preclust(''seterp'',gcf);']; 
    set_scalp    = ['pop_preclust(''setscalp'',gcf);']; 
    set_dipole   = ['pop_preclust(''setdipole'',gcf);'];
    set_moment   = ['pop_preclust(''setmoment'',gcf);'];
    set_ersp     = ['pop_preclust(''setersp'',gcf);']; 
    set_itc      = ['pop_preclust(''setitc'',gcf);']; 
    set_secpca   = ['pop_preclust(''setsec'',gcf);']; 
    
    set_mpcluster   = ['tmp_preclust(''mpcluster'',gcf);']; % nima
    
    help_clusteron = ['pophelp(''std_helpselecton'');']; 
    help_ersp    = ['pophelp(''pop_timef'')'];
    preclust_PCA = ['pop_preclust(''preclustOK'',gcf);'];           
    ersp_params  = ['pop_preclust(''erspparams'',gcf);']; 
    ersp_edit    = ['pop_preclust(''erspedit'',gcf);']; 
    test_ersp    = ['pop_precomp(''testersp'',gcf);']; 
    itc_edit     = 'set(findobj(gcbf, ''tag'', ''ersp_params''), ''string'', get(gcbo, ''string''));';
    ersp_edit    = 'set(findobj(gcbf, ''tag'', ''itc_params'' ), ''string'', get(gcbo, ''string''));';
    
    saveSTUDY  = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
    browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_preclust()''); ' ... 
                  'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    str_name   = ['Build pre-clustering matrix for STUDY set:  ' STUDY.name '' ];
    str_time   = '';
    str_lpfilt = '20';
    help_secpca = [ 'warndlg2(strvcat(''This is the final number of dimensions (otherwise use the sum'',' ...
                    '''of dimensions for all the selected options). See tutorial for more info''), ''Final number of dimensions'');' ];
    cb_help = [ 'warndlg2( strvcat(''When clustering on both time- and location-based measures, take care'',' ... 
                                   '''when performing statistical tests on cluster measures to avoid so-called'',' ...
                                   '''''''double dipping'''' problems in interpreting the results: Adequate nonparametric'',' ...
                                   '''statistical testing should be performed to test their robustness.''));' ];
        
    gui_spec = { ...
    {'style' 'text'       'string' str_name 'FontWeight' 'Bold' 'horizontalalignment' 'left'} ...
    {'style' 'text'       'string' 'Only measures that have been precomputed may be used for clustering'} ...
    {'style' 'text'       'string' 'Mixing time-based and location-based measures might result in statistical double-dipping'} ...
    {'style' 'pushbutton' 'string' 'Help' 'callback' cb_help } ...
    {'style' 'text'       'string' 'Time-based info             PCA              Weight' 'FontWeight' 'Bold'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'spectra_on' 'value' 0 'Callback' set_spectra 'userdata' '1'}  ...
	{'style' 'text'       'string' 'spectra' 'horizontalalignment' 'center' } ...
	{'style' 'edit'       'string' '3' 'tag' 'spectra_PCA' 'enable' 'off' 'userdata' 'specP'} { } ...
    {'style' 'edit'       'string' '1' 'tag' 'spectra_weight' 'enable' 'off' 'userdata' 'specP'} ...
    {'style' 'text'       'string' 'Freq.range [Hz]' 'tag' 'spectra_freq_txt' 'userdata' 'spec' 'enable' 'off' } ...
	{'style' 'edit'       'string' '3  25'  'tag' 'spectra_freq_edit' 'userdata' 'spec' 'enable' 'off' } { } { } ...
    ...
    {'style' 'checkbox'   'string' '' 'tag' 'erp_on' 'value' 0 'Callback' set_erp 'userdata' '1'}  ...
	{'style' 'text'       'string' 'ERPs' 'horizontalalignment' 'center' } ...
    {'style' 'edit'       'string' '3' 'tag' 'erp_PCA' 'enable' 'off' 'userdata' 'erpP'} { } ...
    {'style' 'edit'       'string' '1' 'tag' 'erp_weight' 'enable' 'off' 'userdata' 'erpP'} ...
    {'style' 'text'       'string' 'Time range [ms]' 'tag' 'erp_time_txt' 'userdata' 'erp' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_time 'tag' 'erp_time_edit' 'userdata' 'erp' 'enable' 'off' }...
    {'style' 'text'       'string' 'Lowpass [Hz]' 'tag' 'erp_filter_txt' 'userdata' 'erp' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_lpfilt 'tag' 'erp_filter_edit' 'userdata' 'erp' 'enable' 'off' }...
    ...
    {'style' 'checkbox'   'string' '' 'tag' 'ersp_on' 'value' 0 'Callback' set_ersp 'userdata' '1'} ...
	{'style' 'text'       'string' 'ERSPs' 'horizontalalignment' 'center' } ...
    {'style' 'edit'       'string' '3' 'tag' 'ersp_PCA' 'enable' 'off' 'userdata' 'erspP'} { }   ...
	{'style' 'edit'       'string' '1' 'tag' 'ersp_weight' 'enable' 'off' 'userdata' 'erspP'} ...
    {'style' 'text'       'string' 'Time range [ms]' 'tag' 'ersp_time_txt' 'userdata' 'ersp' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_time 'tag' 'ersp_time_edit' 'userdata' 'ersp' 'enable' 'off' } ...
    {'style' 'text'       'string' 'Freq. range [Hz]' 'tag' 'ersp_time_txt' 'userdata' 'ersp' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_time 'tag' 'ersp_freq_edit' 'userdata' 'ersp' 'enable' 'off' } ...
    ...
    {'style' 'checkbox'   'string' '' 'tag' 'itc_on' 'value' 0 'Callback' set_itc 'userdata' '1'} ...
	{'style' 'text'       'string' 'ITCs' 'horizontalalignment' 'center' } ...
    {'style' 'edit'       'string' '3' 'tag' 'itc_PCA' 'enable' 'off' 'userdata' 'itcP'} { }   ...
	{'style' 'edit'       'string' '1' 'tag' 'itc_weight' 'enable' 'off' 'userdata' 'itcP'} ...
    {'style' 'text'       'string' 'Time range [ms]' 'tag' 'itc_time_txt' 'userdata' 'itcP' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_time 'tag' 'itc_time_edit' 'userdata' 'itcP' 'enable' 'off' } ...
    {'style' 'text'       'string' 'Freq. range [Hz]' 'tag' 'itc_time_txt' 'userdata' 'itcP' 'enable' 'off' } ...
	{'style' 'edit'       'string' str_time 'tag' 'itc_freq_edit' 'userdata' 'itcP' 'enable' 'off' } ...
    ...
    {'style' 'text'       'string' 'Location-based info      PCA              Weight' 'FontWeight' 'Bold'} ...
	{'style' 'checkbox'   'string' '' 'tag' 'dipole_on' 'value' 0 'Callback' set_dipole 'userdata' '1'} ...
	{'style' 'text'       'string' 'dipole locations' 'HorizontalAlignment' 'center' } ...
	{'style' 'text'       'string' '3' 'enable' 'off' 'userdata' 'dipoleP' } { }  ...
	{'style' 'edit'       'string' '1' 'tag' 'locations_weight' 'enable' 'off' 'userdata' 'dipoleP'} {} {} {} {} ...
    ...
	{'style' 'checkbox'   'string' '' 'tag' 'moment_on' 'value' 0 'Callback' set_moment 'userdata' '1'} ...
	{'style' 'text'       'string' 'dipole orient.' 'HorizontalAlignment' 'center' } ...
	{'style' 'text'       'string' '3' 'enable' 'off' 'userdata' 'dipoleM' } { }  ...
	{'style' 'edit'       'string' '1' 'tag' 'moment_weight' 'enable' 'off' 'userdata' 'dipoleM'} ...
    {'style' 'text'       'string' 'Amplitude & polarity is ignored'} {} {} {} ...
    ...
    {'style' 'checkbox'   'string' '' 'tag' 'scalp_on' 'value' 0 'Callback' set_scalp 'userdata' '1'} ...
	{'style' 'text'       'string' 'scalp maps' 'HorizontalAlignment' 'center' } ...
	{'style' 'edit'       'string' '3' 'tag' 'scalp_PCA' 'enable' 'off' 'userdata' 'scalpP'} { }   ...
    {'style' 'edit'       'string' '1' 'tag' 'scalp_weight' 'enable' 'off' 'userdata' 'scalpP'} ...
    {'style' 'popupmenu'  'string' scalp_options 'value' 1 'tag' 'scalp_choice' 'enable' 'off' 'userdata' 'scalp' } {} ...
    {'style' 'checkbox'   'string' 'Absolute values' 'value' 1 'tag'  'scalp_absolute' 'enable' 'off' 'userdata' 'scalp' } {} ...
    ...
    {} };
  

%    {'link2lines' 'style'  'text'   'string' '' } {} {} {} ...
%    {'style' 'text'       'string' 'Time/freq. parameters' 'tag' 'ersp_push' 'value' 1 'enable' 'off' 'userdata' 'ersp' 'Callback' ersp_params} ...
%    {'style' 'edit'       'string' erspparams_str 'tag' 'ersp_params' 'enable' 'off' 'userdata' 'ersp' 'Callback' ersp_edit}...
%    {'style' 'text'       'string' 'Time/freq. parameters' 'tag' 'itc_push' 'value' 1 'enable' 'off' 'userdata' 'itc' 'Callback' ersp_params} ...
%    {'style' 'edit'       'string' erspparams_str 'tag' 'itc_params' 'enable' 'off' 'userdata' 'itc' 'Callback' itc_edit}%    {'style' 'checkbox'   'string' '' 'tag' 'preclust_PCA'  'Callback' preclust_PCA 'value' 0} ...
%    {'style' 'text'       'string' 'Do not prepare dataset for clustering at this time.' 'FontWeight' 'Bold'  } {} ...

    fig_arg{1} = { ALLEEG STUDY cls };
    geomline = [0.5 2 1 0.5 1 2 1 2 1 ];
    geometry = { [1] [1] [1 0.2] ...
                 [3] geomline geomline  geomline geomline [1] geomline [0.5 2 1 0.5 1 5.4 0.2 0.2 0.2 ] [0.5 2 1 0.5 1 2.9 .1 2.9 .1 ] [1] };
    geomvert = [ 1 1 1 1 1 1 1 1 1 1 1 1 0.5 ];

    %if length(show_options) < 3
    %    gui_spec(2:6) = { {} ...
    %        { 'style' 'text'      'string' [ 'Among the pre-selected components (Edit study),' ...
    %                        'remove those which dipole res. var, exceed' ] 'tag' 'dipole_select_on' }  ...
    %        {'style' 'edit'       'string' '0.15' 'horizontalalignment' 'center' 'tag' 'dipole_rv'} ...
    %        {'style' 'text'       'string'  '(empty=all)'} {} };
    %    geometry{3} = [2.5 0.25 0.4];
    %    geomvert(3) = 1;
    %end
    disp('Even though the graphics have changed, this function remains 100% backward compatible.');
	[preclust_param, userdat2, strhalt, os] = inputgui( 'geometry', geometry, 'uilist', gui_spec, 'geomvert', geomvert, ...
                                                      'helpcom', ' pophelp(''std_preclust'')', ...
                                                      'title', 'Select and compute component measures for later clustering -- pop_preclust()', ...
                                                      'userdata', fig_arg);	
	if isempty(preclust_param), return; end
    
    options = { STUDY, ALLEEG, 1};
    
    % Spectrum option is on
    % --------------------
    if os.spectra_on== 1 
        options{end+1} = {  'spec' 'npca' str2num(os.spectra_PCA) ...
                            'weight' str2num(os.spectra_weight)  'freqrange' str2num(os.spectra_freq_edit) };
    end
    
    % ERP option is on
    % ----------------
    if os.erp_on == 1 
        options{end+1} = { 'erp' 'npca' str2num(os.erp_PCA) ...
                         'weight' str2num(os.erp_weight) 'timewindow' str2num(os.erp_time_edit), 'erpfilter', os.erp_filter_edit };
    end
    
    % Scalp maps option is on
    % ----------------------
    if os.scalp_on == 1 
        if os.scalp_absolute %absolute maps
            abso = 1;
        else abso = 0;
        end
        if (os.scalp_choice == 2)  %Laplacian scalp maps
            options{end+1} = { 'scalpLaplac' 'npca' str2num(os.scalp_PCA) ...
                               'weight' str2num(os.scalp_weight) 'abso' abso };
        elseif (os.scalp_choice == 3)  %Gradient scalp maps
            options{end+1} = { 'scalpGrad' 'npca' str2num(os.scalp_PCA) ...
                               'weight' str2num(os.scalp_weight) 'abso' abso };
        elseif (os.scalp_choice == 1) %scalp map case
            options{end+1} = { 'scalp' 'npca' str2num(os.scalp_PCA) ...
                               'weight' str2num(os.scalp_weight) 'abso' abso };
        end
    end
    
    % Dipole option is on
    % -------------------
    if os.dipole_on == 1 
        options{end+1} = { 'dipoles' 'weight' str2num(os.locations_weight) };
    end
    if os.moment_on == 1 
        options{end+1} = { 'moments' 'weight' str2num(os.moment_weight) };
    end
    
    % ERSP option is on
    % -----------------
    if os.ersp_on  == 1 
        options{end+1} = { 'ersp' 'npca' str2num(os.ersp_PCA) 'freqrange' str2num(os.ersp_freq_edit) ...
                          'timewindow' str2num(os.ersp_time_edit) 'weight' str2num(os.ersp_weight) };
    end
    
    % ITC option is on 
    % ----------------
    if os.itc_on  == 1 
        options{end+1} = { 'itc' 'npca' str2num(os.itc_PCA) 'freqrange' str2num(os.itc_freq_edit) 'timewindow' ...
                           str2num(os.itc_time_edit) 'weight' str2num(os.itc_weight) };
    end       
    
    % evaluate command
    % ----------------
    if length(options) == 3
        warndlg2('No measure selected: aborting.'); 
        return; 
    end
    
    [STUDY ALLEEG] = std_preclust(options{:});
    com = sprintf('[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, %s);', vararg2str(options(3:end)));
    
    % save updated STUDY to the disk
    % ------------------------------
%     if os.saveSTUDY == 1 
%         if ~isempty(os.studyfile)
%             [filepath filename ext] = fileparts(os.studyfile);
%             STUDY.filename = [ filename ext ];
%             STUDY.filepath = filepath;
%         end
%         STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
%         com = sprintf('%s\nSTUDY = pop_savestudy(STUDY, ALLEEG, %s);',  com, ...
%                       vararg2str( { 'filename', STUDY.filename, 'filepath', STUDY.filepath }));
%     end
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
            
        case 'mpcluster' % nima
            mpclust =  get(findobj('parent', hdl, 'tag', 'mpclust'), 'value');
            if mpclust
                set(findobj('parent', hdl, 'tag', 'spectra_PCA'), 'visible','off');
                set(findobj('parent', hdl, 'tag', 'spectra_weight'), 'visible','off');
                set(findobj('parent', hdl, 'tag',  'erp_PCA' ), 'visible','off');
                set(findobj('parent', hdl, 'tag','erp_weight' ), 'visible','off');
                set(findobj('parent', hdl, 'tag','locations_weight'), 'visible','off');
                set(findobj('parent', hdl, 'tag', 'scalp_PCA'), 'visible','off');
                set(findobj('parent', hdl, 'tag','scalp_weight'), 'visible','off');
                set(findobj('parent', hdl, 'tag','ersp_PCA'), 'visible','off');
                set(findobj('parent', hdl, 'tag','ersp_weight' ), 'visible','off');
                set(findobj('parent', hdl, 'tag', 'itc_PCA'), 'visible','off');
                set(findobj('parent', hdl, 'tag','itc_weight'), 'visible','off');

                set(findobj('parent', hdl, 'tag','sec_PCA'), 'visible','off');
                set(findobj('parent', hdl, 'tag','sec_on'), 'visible','off');
                set(findobj('parent', hdl, 'userdata' ,'dipoleP'), 'visible','off');
                set(findobj('parent', hdl, 'userdata' ,'dipoleM'), 'visible','off');
                set(findobj('parent', hdl, 'string','Final dimensions'), 'visible','off');
                set(findobj('parent', hdl, 'tag','finalDimHelp' ), 'visible','off');
                set(findobj('parent', hdl, 'tag','spectra_freq_txt'), 'visible','off');
                set(findobj('parent', hdl, 'tag','spectra_freq_edit'), 'visible','off');

                %% these are made invisible for now,  but in future we might use them in the new method
                set(findobj('parent', hdl, 'tag','erp_time_txt'), 'visible','off');
                set(findobj('parent', hdl, 'tag','erp_time_edit'), 'visible','off');
                set(findobj('parent', hdl, 'tag','scalp_choice'), 'visible','off');
                set(findobj('parent', hdl, 'tag','scalp_absolute'), 'visible','off');
                set(findobj('parent', hdl, 'tag','ersp_time_txt'), 'visible','off');
                set(findobj('parent', hdl, 'tag','ersp_time_edit'), 'visible','off');
                set(findobj('parent', hdl, 'tag','ersp_freq_edit'), 'visible','off');
                set(findobj('parent', hdl, 'tag','itc_time_txt'), 'visible','off');
                set(findobj('parent', hdl, 'tag','itc_time_edit'), 'visible','off');
                set(findobj('parent', hdl, 'tag','itc_freq_edit'), 'visible','off');

                set(findobj('parent', hdl, 'string','Measures                         Dims.   Rel. Wt.'), 'string','Measures');
            else
                set(findobj('parent', hdl, 'tag', 'spectra_PCA'), 'visible','on');
                set(findobj('parent', hdl, 'tag', 'spectra_weight'), 'visible','on');
                set(findobj('parent', hdl, 'tag',  'erp_PCA' ), 'visible','on');
                set(findobj('parent', hdl, 'tag','erp_weight' ), 'visible','on');
                set(findobj('parent', hdl, 'tag','locations_weight'), 'visible','on');
                set(findobj('parent', hdl, 'tag', 'scalp_PCA'), 'visible','on');
                set(findobj('parent', hdl, 'tag','scalp_weight'), 'visible','on');
                set(findobj('parent', hdl, 'tag','ersp_PCA'), 'visible','on');
                set(findobj('parent', hdl, 'tag','ersp_weight' ), 'visible','on');
                set(findobj('parent', hdl, 'tag', 'itc_PCA'), 'visible','on');
                set(findobj('parent', hdl, 'tag','itc_weight'), 'visible','on');

                set(findobj('parent', hdl, 'tag','sec_PCA'), 'visible','on');
                set(findobj('parent', hdl, 'tag','sec_on'), 'visible','on');
                set(findobj('parent', hdl, 'userdata' ,'dipoleP'), 'visible','on');
                set(findobj('parent', hdl, 'userdata' ,'dipoleM'), 'visible','on');
                set(findobj('parent', hdl, 'string','Final dimensions'), 'visible','on');
                set(findobj('parent', hdl, 'tag','finalDimHelp' ), 'visible','on');
                set(findobj('parent', hdl, 'tag','spectra_freq_txt'), 'visible','on');
                set(findobj('parent', hdl, 'tag','spectra_freq_edit'), 'visible','on');

                %% these are made invisible for now,  but in future we might use them in the new method
                set(findobj('parent', hdl, 'tag','erp_time_txt'), 'visible','on');
                set(findobj('parent', hdl, 'tag','erp_time_edit'), 'visible','on');
                set(findobj('parent', hdl, 'tag','scalp_choice'), 'visible','on');
                set(findobj('parent', hdl, 'tag','scalp_absolute'), 'visible','on');
                set(findobj('parent', hdl, 'tag','ersp_time_txt'), 'visible','on');
                set(findobj('parent', hdl, 'tag','ersp_time_edit'), 'visible','on');
                set(findobj('parent', hdl, 'tag','ersp_freq_edit'), 'visible','on');
                set(findobj('parent', hdl, 'tag','itc_time_txt'), 'visible','on');
                set(findobj('parent', hdl, 'tag','itc_time_edit'), 'visible','on');
                set(findobj('parent', hdl, 'tag','itc_freq_edit'), 'visible','on');

                set(findobj('parent', hdl, 'string','Measures to Cluster on:'), 'string','Load                                  Dims.   Rel. Wt.');
                set(findobj('parent', hdl, 'string','Measures'), 'string', 'Measures                         Dims.   Rel. Wt.');
            end

                
%             set_mpcluster =  get(findobj('parent', hdl, 'tag', 'spectra_on'), 'value'); 
%             set(findobj('parent', hdl, 'userdata', 'spec'), 'enable', fastif(set_spec,'on','off'));
%             PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
%             if PCA_on
%                 set(findobj('parent', hdl, 'userdata', 'specP'), 'enable', 'off');
%             else
%                 set(findobj('parent', hdl, 'userdata', 'specP'), 'enable', fastif(set_spec,'on','off'));
%             end         
            
            
            
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
        case 'setmoment'
            moment_on =  get(findobj('parent', hdl, 'tag', 'moment_on'), 'value'); 
            set(findobj('parent', hdl, 'userdata', 'moment'), 'enable', fastif(moment_on,'on','off'));
            PCA_on = get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value');
            if PCA_on
               set(findobj('parent', hdl, 'userdata', 'dipoleM'), 'enable','off');
           else
               set(findobj('parent', hdl, 'userdata', 'dipoleM'), 'enable', fastif(moment_on,'on','off'));
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
             if  (set_itc && (~set_ersp) )
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
            if  (set_ersp && (~set_itc) )
                set(findobj('parent', hdl,'tag', 'ersp_push'), 'enable', 'on');
                set(findobj('parent', hdl,'tag', 'ersp_params'), 'enable', 'on');
            end
        case 'setsec'
            set_sec =  get(findobj('parent', hdl, 'tag', 'sec_on'), 'value'); 
            set(findobj('parent', hdl,'userdata', 'sec'), 'enable', fastif(set_sec,'on','off'));
        case 'erspparams'
            ersp = userdat{2};
            [ersp_paramsout, erspuserdat, strhalt, erspstruct] = inputgui( { [1] [3 1] [3 1] [3 1] [3 1] [3 1] [1]}, ...
                    { {'style' 'text' 'string' 'ERSP and ITC time/freq. parameters' 'FontWeight' 'Bold'} ...
                    {'style' 'text' 'string' 'Frequency range [Hz]' 'tag' 'ersp_freq' } ...
                    {'style' 'edit' 'string' ersp.f 'tag' 'ersp_f' 'Callback' ERSP_timewindow } ...
                    {'style' 'text' 'string' 'Wavelet cycles (see >> help timef())' 'tag' 'ersp_cycle' } ...
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
                     ersp.c '], ''alpha'', ' ersp.a ', ''padratio'', ' ersp.p ', ''tlimits'', [' ersp.t ']']);
                set(findobj('parent', hdl, 'tag', 'itc_params'), 'string', ...
                    ['                                                             ''frange'', [' ersp.f '], ''cycles'', [' ...
                     ersp.c '], ''alpha'', ' ersp.a ', ''padratio'', ' ersp.p ', ''tlimits'', [' ersp.t ']']);
                set(hdl, 'userdat',userdat); 
            end
       case 'preclustOK'
           set_PCA =  get(findobj('parent', hdl, 'tag', 'preclust_PCA'), 'value'); 
           set_ersp =  get(findobj('parent', hdl, 'tag', 'ersp_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'erspP'), 'enable', fastif(~set_PCA && set_ersp,'on','off'));
           set_itc =  get(findobj('parent', hdl, 'tag', 'itc_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'itcP'), 'enable', fastif(~set_PCA && set_itc,'on','off'));
           set_erp =  get(findobj('parent', hdl, 'tag', 'erp_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'erpP'), 'enable', fastif(~set_PCA && set_erp,'on','off'));
           set_spec =  get(findobj('parent', hdl, 'tag', 'spectra_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'specP'), 'enable', fastif(~set_PCA && set_spec,'on','off'));
           set_scalp =  get(findobj('parent', hdl, 'tag', 'scalp_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'scalpP'), 'enable', fastif(~set_PCA && set_scalp,'on','off'));
           set_dipole =  get(findobj('parent', hdl, 'tag', 'dipole_on'), 'value'); 
           set(findobj('parent', hdl,'userdata', 'dipoleP'), 'enable', fastif(~set_PCA && set_dipole,'on','off'));
           set(findobj('parent', hdl,'tag', 'chosen_component'), 'enable', fastif(~set_PCA,'on','off'));
           set(findobj('parent', hdl,'tag', 'dipole_rv'), 'enable', fastif(~set_PCA,'on','off'));
           set(findobj('parent', hdl,'tag', 'compstd_str'), 'enable', fastif(~set_PCA,'on','off'));
    end
end
STUDY.saved = 'no';
