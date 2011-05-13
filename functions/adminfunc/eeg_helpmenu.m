% eeg_helpmenu() - Call the help file for EEGLAB menus       
%       
% Usage:
%   To call a menu of eeglab functions
%       >> eeg_helpmenu;
%   To open a help window with the specific called file.
%   Equivalent to calling >> pophelp( 'filename' );
%        >> eeg_helpmenu( 'filename' );
%       
% Author: Arnaud Delorme CNL / Salk Institute 2001     
%       
% See also: eeglab()       
       
% Copyright (C) 2001 Arnaud Delorme Salk Institute arno@salk.edu     
%       
% This program is free software; you can redistribute it and/or modify       
% it under the terms of the GNU General Public License as published by       
% the Free Software Foundation; either version 2 of the License or      
% (at your option) any later version.       
%       
% This program is distributed in the hope that it will be useful       
% but WITHOUT ANY WARRANTY; without even the implied warranty of       
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
% GNU General Public License for more details.       
%       
% You should have received a copy of the GNU General Public License       
% along with this program; if not write to the Free Software      
% Foundation Inc. 59 Temple Place Suite 330 Boston MA  02111-1307  USA  
       
function eeg_helpmenu( filename );       
       
allmenus = { ...       
'File' ''      ...
'   Import data' ''  ...    
'      From ASCII/float file or Matlab array' 'pop_importdata'     ... %       Import Matlab data array      
'      From continuous or seg. EGI .RAW file' 'pop_readegi'        ... %       Read continuous or seg. .EGI or .RAW datafile      
'      From Multiple seg. EGI .RAW files'     'pop_readsegegi'     ... %       Read multiple seg. .EGI or .RAW datafile
'      From BCI2000 ASCII file'               'pop_loadbci'        ... %       Read .BCI datafile (Snapmaster)      
'      From Snapmaster .SMA file'             'pop_snapread'       ... %       Read .SMA datafile (Snapmaster)      
'      From Neuroscan .CNT file'              'pop_loadcnt'        ... %       Read .CNT datafile (Neuroscan continuous)      
'      From Neuroscan .EEG file'              'pop_loadeeg'        ... %       Read .EEG datafile (Neuroscan epochs)      
'      From Biosemi .BDF file'                'pop_biosig'         ... %       Read .BDF datafile
'      From .EDF file'                        'pop_biosig'         ... %       Read .EDF datafile
'      From ANT EEProbe .CNT file'            'pop_loadeep'        ... %       Read .CNT datafile (ANT EEProbe)
'      From ANT EEProbe .AVR file'            'pop_loadeep_avg'    ... %       Read .AVR datafile (ANT EEProbe)
'      From .BDF file (backup function)'      'pop_readbdf'        ... %       Read .BDF datafile
'      From Brain Vis. Rec. .VHDR file'       'pop_loadbv'         ... %       Read .VHDR datafile (Brain Vision Rec.)
'      From Brain Vis. Anal. Matlab file'     'pop_loadbva'        ... %       Read Matlab datafile (Brain Vision Analysis)
'      From CTF Folder (MEG)'                 'ctf_folder'         ... %       Import MEG datafiles from chosen CTF folder.
'      From ERPSS .RAW or .RDF file'          'pop_read_erpss'     ...
'      From INStep .ASC file'                 'pop_loadascinstep'  ... %       Read .ASC datafile (INStep)
'      From 4D .M4D pdf file'                 'pop_read4d'         ... %       Read .M4D datafile (4D pdf)
'      From other formats using FILEIO'       'pop_fileio'         ... %       Import data using FILEIO
'      From other formats using BIOSIG'       'pop_biosig'         ...
'   Import epoch info' ''      ...
'      From Matlab array or ASCII file'       'pop_importepoch'    ... %       Import Matlab array or ASCII file 
'      From Neuroscan .DAT file'              'pop_loaddat'        ... %       Import .DAT info file (Neuroscan epochs)      
'   Import event info' ''      ...
'      From  Matlab array or ASCII file'      'pop_importevent'    ... %       Import Matlab array or ASCII file      
'      From data channel'                     'pop_chanevent'      ... %       Import events from channels      
'      From Presentation .LOG file'           'pop_importpres'     ... %       Import .LOG event file (Presentation)  
'      From Neuroscan .EV2 file'              'pop_importev2'      ... %       Import .EV2 even file (Neuroscan)
'   Export ' '' ...
'      Data and ICA to text file'             'pop_export'         ... %       Export Data and ICA to text file
'      Weight matrix to text file'            'pop_expica'         ... %       Export weight matrix to text file
'      Inverse weight matrix to text file'    'pop_expica'         ... %       Export inverse weight matrix to text file
'      Data to EDF/BDF/GDF file'              'pop_writeeeg'       ... %       Export Data to EDF/BDF/GDF data formats
'      Write Brain Vis. exchange format file' '            '       ... 
'   Load existing dataset'                    'pop_loadset'        ... %       Load dataset      
'   Save current dataset'                     'pop_saveset'        ... %       Save dataset      
'   Save datasets'                            'pop_saveset'        ... %       Save dataset      
'   Clear dataset(s)'                         'pop_delset'         ... %       Clear dataset(s)      
'   Create Study ' '' ...
'      Using all loaded datasets'             'pop_study'          ... %       Create study using currently loaded datasets
'      Browse for datasets'                   'pop_study'          ... %       Find datasets to use in Study
'   Load existing study'                      'pop_loadstudy'      ... %       Load study      
'   Save current study'                       'pop_savestudy'      ... %       Save study      
'   Save current study as'                    'pop_savestudy'      ... %       Save current study as      
'   Clear study'                              'pop_delset'         ... %       Clear study      
'   Maximize memory'                          'pop_editoptions'    ... %       Maximize memory     
'   History Scripts' '' ...
'      Save dataset history script'           'pop_saveh'          ... %       Save current dataset history
'      Save session history script'           'pop_saveh'          ... %       Save current session history
'      Run Script'                            'pop_runscript'      ... %       Run History script     
'   Quit' ''      ...
'Edit' ''      ...
'   Dataset info'                             'pop_editset'        ... %       Edit dataset info      
'   Event fields'                             'pop_editeventfield' ... %       Edit event fields      
'   Event values'                             'pop_editeventvals'  ... %       Edit event values      
'   About this dataset'                       'pop_comments'       ... %       About this dataset      
'   Channel locations'                        'pop_chanedit'       ... %       Edit channels      
'   Select data'                              'pop_select'         ... %       Select data      
'   Select data using events'                 'pop_rmdat'          ... %       Select data using events
'   Select epochs or events'                  'pop_selectevent'    ... %       Select epochs or events      
'   Copy current dataset'                     'pop_copyset'        ... %       Copy current dataset      
'   Append datasets'                          'pop_mergeset'       ... %       Append another dataset      
'   Delete dataset(s)'                        'pop_delset'         ... %       Delete dataset(s)      
'Tools' ''      ...
'   Change sampling rate'                     'pop_resample'       ... %       Change sampling rate      
'   Filter the data' '' ...
'      Basic FIR filter'                      'pop_eegfilt'        ... %       Filter the data      
'      Short IIR filter'                      'pop_iirfilt'        ... %       Short IIR filter
'      ERPLAB Butterworth Digital Filter'     'pop_ERPLAB_butter1' ... %       ERPLAB Butterworth filter
'      ERPLAB Polynomial Detrending'          'pop_ERPLAB_polydetrend' ... % Polynomial detrending
'   Re-referencing'                           'pop_reref'          ... %       Average reference      
'   Interpolate electrodes'                   'pop_interp'         ...
'   Reject continuous data by eye'            'pop_eegplot'        ... %       Reject continuous data      
'   Extract epochs'                           'pop_epoch'          ... %       Extract epochs      
'   Remove baseline'                          'pop_rmbase'         ... %       Remove baseline
'   Run ICA'                                  'pop_runica'         ... %       Run ICA      
'   Remove components'                        'pop_subcomp'        ... %       Remove components 
'   Automatic channel rejection'              'pop_rejchan'        ... %       Reject channels
'   Automatic epoch rejection'                'pop_autorej'        ... %       Reject epochs
'   Reject data epochs' ''  ...    
'      Reject data (all methods)'             'pop_rejmenu'        ... %       Reject data (all methods)      
'      Reject by inspection'                  'pop_eegplot'        ... %       Reject by inspection      
'      Reject extreme values'                 'pop_eegthresh'      ... %       Reject extreme values      
'      Reject by linear trend/variance'       'pop_rejtrend'       ... %       Reject by linear trend/variance     
'      Reject by probability'                 'pop_jointprob'      ... %       Reject by probability      
'      Reject by kurtosis'                    'pop_rejkurt'        ... %       Reject by kurtosis      
'      Reject by spectra'                     'pop_rejspec'        ... %       Reject by spectra      
'      Reject marked epochs'                  'pop_rejepoch'       ... %       Reject labeled epochs      
'   Reject using ICA' '' ...     
'      Reject components by map'              'pop_selectcomps'    ... %       Reject components by map      
'      Reject data (all methods)'             'pop_rejmenu'        ... %       Reject data (all methods)      
'      Reject by inspection'                  'pop_eegplot'        ... %       Reject by inspection      
'      Reject extreme values'                 'pop_eegthresh'      ... %       Reject extreme values      
'      Reject by linear trend/variance'       'pop_rejtrend'       ... %       Reject by linear trend/variance      
'      Reject by probability'                 'pop_jointprob'      ... %       Reject by probability      
'      Reject by kurtosis'                    'pop_rejkurt'        ... %       Reject by kurtosis      
'      Reject by spectra'                     'pop_rejspec'        ... %       Reject by spectra      
'      Reject marked epochs'                  'pop_rejepoch'       ... %       Reject labeled epochs      
'   Run AMICA'  ''  ...
'      Run AMICA in parallel'                 'pop_runica'         ... %       Run AMICA in parallel (plugin)
'      Run AMICA locally'                     'pop_runica'         ... %       Run AMICA locally (plugin)
'      Load AMICA component'                  'pop_loadmodout'     ... %       Load AMICA component (plugin)
'   Run CICA'                                 ''                   ...
'   Locate dipoles using DIPFIT 2.x' '' ...
'      Head model and settings'               'pop_dipfit_settings'... %       DIPFIT Plugin settings
'      Coarse fit (grid scan)'                ''                   ...
'      Fine fit (iterative)'                  'pop_dipfit_nonlinear'... %      Fine fit (iterative)
'      Autofit (coarse fit, fine fit & plot)' ''                   ... 
'      Plot component dipoles'                ''                   ...
'   Peak detection using EEG toolbox'         ''                   ...
'   FMRIB Tools' '' ...
'      FASTR: Remove FMRI gradient artifacts' 'pop_fmrib_fastr'    ... %       Remove FMRI gradient artifacts
'      Detect QRS events'                     'pop_fmrib_qrsdetect'... %       Detect QRS events
'      Remove pulse artifacts'                'pop_fmrib_pas'      ... %       Remove pulse artifacts
'   BEM4 SCCN beta'  '' ...
'      Forward Problem Solution'              'pop_forward problem_solution' ... % Not currently working
'   Locate dipoles using LORETA' '' ...
'      Export components to LORETA'           ''                   ...
'      Import dipoles from LORETA'            ''                   ...
'      Plot dipoles on LORETA head'           ''                   ...
'Plot' ''      ...
'   Channel locations' '' ...     
'      By name'                               'topoplot'           ... %       Electrodes      
'      By number'                             'topoplot'           ... %       Numbers      
'   Channel data (scroll)'                    'pop_eegplot'        ... %       Scroll EEG data      
'   Channel spectra and maps'                 'pop_spectopo'       ... %       Channel spectra and maps      
'   Channel properties'                       'pop_prop'           ... %       Component properties      
'   Channel ERP image'                        'pop_erpimage'       ... %       Channel ERP image      
'   Channel ERPs' ''      ...
'      With scalp maps'                       'pop_timtopo'        ... %       ERP and scalp maps      
'      In scalp array'                        'pop_plottopo'       ... %       ERP in scalp array      
'      In rect. array'                        'pop_plotdata'       ... %       ERP in rect. array      
'   ERP maps' ''      ...
'      As 2-D scalp maps'                     'pop_topoplot'       ... %       ERP scalp maps      
'      As 3-D head plots'                     'pop_headplot'       ... %       ERP head plots      
'   Sum/Compare ERPs'                         'pop_comperp'        ... %       Compare ERPs      
'   Component activations (scroll)'           'pop_eegplot'        ... %       Scroll component activations      
'   Component spectra and maps'               'pop_spectopo'       ... %       Channel spectra and maps with components      
'   Component maps' ''      ...
'      As 2-D scalp maps'                     'pop_topoplot'       ... %       Component scalp maps      
'      As 3-D head plots'                     'pop_headplot'       ... %       Component head plots      
'   Component properties'                     'pop_prop'           ... %       Component properties      
'   Component ERP image'                      'pop_erpimage'       ... %       Component ERP image      
'   Component ERPs' ''      ...
'      With component maps'                   'pop_envtopo'        ... %       Largest ERP components      
'      With comp. maps (compare)'             'pop_envtopo'        ... %       Largest ERP components
'      In rectangular array'                  'pop_plotdata'       ... %       Comp. ERP time courses  
'   Sum/Compare comp. ERPs'                   'pop_comperp'        ... %       Compare ERPs 
'   Data statistics' ''      ...
'      Channel statistics'                    'pop_signalstat'     ...
'      Component statistics'                  'pop_signalstat'     ...
'      Event statistics'                      'pop_eventstat'      ...
'   Time-frequency transforms' ''  ...
'      Channel time-frequency'                'pop_timef'          ... %       Channel time-frequency      
'      Channel cross-coherence'               'pop_crossf'         ... %       Channel cross-coherence      
'      Component time-frequency'              'pop_timef'          ... %       Component time-frequency EEGLAB'));     
'      Component cross-coherence'             'pop_crossf'         ... %       Component cross-coherence  
'   Cluster dataset ICs'                      'pop_miclust'        ... %       Cluster ICs
'Study' ''     ...
'   Edit study info'                          'pop_study'          ...  
'   Precompute channel measures'              'pop_precomp'        ...
'   Plot channel measures'                    'pop_chanplot'       ... 
'   Precompute component measures'            'pop_precomp'        ... 
'   Build preclustering array'                'pop_preclust'       ... 
'   Cluster components'                       'pop_clust'          ... 
'   Edit/plot clusters'                       'pop_clustedit'      ... 
'Datasets' ''   ...
'   Current/Active Datasets (listed as selectable items)' ''       ... 
'   Select multiple datasets'                  ''                  ...
'   Select the study set'                      ''                  ... %      Only selectable if study is loaded
'Help' ''    ...
'   About EEGLAB'                              'eeglab'            ...
'   About EEGLAB Help'                         'eeg_helphelp'      ... 
'   EEGLAB menus'                              'eeg_helpmenu'      ... 
'   EEGLAB functions'                          ''                  ...
'       Toolbox functions'                     ''                  ...
'       Signal processing functions'           'eeg_helpsigproc'   ... 
'       Interactive pop_functions'             'eeg_helppop'       ...
'   EEGLAB advanced'                           ''                  ...
'       Dataset structure'                     'eeg_checkset'      ...
'       Admin functions'                       'eeg_helpadmin'     ... 
'   Web tutorial'                              ''                  ...
'   Email EEGLAB'                              ''                  ...
                           
};       

text  = allmenus( 1:2:length(allmenus));
command = allmenus( 2:2:length(allmenus));

for index = 1:length(command)
	if ~isempty(command{index})
		command{index} = [ 'pophelp(''' command{index} ''');' ];
	end;
end;

textgui( text, command,'fontsize', 15, 'fontname', 'times', 'linesperpage', 18, 'title',strvcat( 'Interactive (pop_) functions', '(Click on blue text for help)'));
try 
	icadefs; set(gcf, 'COLOR', BACKCOLOR);
	h = findobj('parent', gcf, 'style', 'slider');
	set(h, 'backgroundcolor', GUIBACKCOLOR);
catch, end;

if nargin > 0
	close(gcf);
    pophelp(filename);
	%command = allmenus( 2:2:length(allmenus));
	%for index = 1:length(text)
	%	if isempty(allmenus{index*2})
	%		tmpparam{index} = {  '' allmenus{(index-1)*2+1} };
	%	else
	%		tmpparam{index} = {  [ allmenus{index*2} '.m' ] allmenus{(index-1)*2+1} };
	%	end;
	%end;
 	%makehtml( tmpparam, '/home/www/eeglab/allfunctions', 'outputfile', 'tempmenu.html', 'mainonly', 'on');
    %}
end;

return;


% $$$ [textmenu nblines l] = getallmenus(findobj('tag', 1) ~= 32 )
% $$$ fontsize{ index } = fontsize{ index }+1;
% $$$ fontweight{ index } ='bold';
% $$$ 
% $$$ fontsize(1:length( textmenu )) = { 15 }; 6) ~= 32 fontweight{ index } ='bold'; end;     
% $$$ fontweight(1:length( textmenu )) = { 'normal' }; :));      
% $$$ for index = 1:length( textmenu )       
% $$$ if textmenu(index       
% $$$ if textmenu(index       
% $$$ a = deblank(textmenu(index       
% $$$ if index > length( commands )       
% $$$ commands = { commands{:} [] }; 1:length(a)) = a;      
% $$$ end;       
% $$$ if ~isempty(commands{index})       
% $$$ a = [ a ' -- ' commands{index} '()'];       
% $$$ textmenu(index       
% $$$ commands{index} = [ ''pophelp(''' commands{index} ''');' ]; ...      
% $$$ end; click on blue text below)'  ' ...    
% $$$ end; :));      
% $$$        
% $$$ textmenu = strvcat('Functions called through the EEGLAB menu'       
% $$$ (For help message       
% $$$ textmenu(1:end-1 :) commands ...    
% $$$ fontsize   = { 16 15 15 fontsize{:} }; fontsize fontweight' fontweight fontname' times' linesperpage' 17 );
% $$$ fontweight = { 'bold' 'normal' 'normal' fontweight{:} };       
% $$$ commands   = { '' ' ' '' commands{:} };       
% $$$ textgui(textmenu(1:end-1       
% $$$ fontsize'       
% $$$ return;       
