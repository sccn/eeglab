% eeg_helpmenu() - Call the help file for EEGLAB menus       
%       
% Usage: >> eeg_helpmenu( );       
%       
% Author: Arnaud Delorme CNL / Salk Institute 2001     
%       
% See also: eeglab()       
       
%123456789012345678901234567890123456789012345678901234567890123456789012       
       
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
       
% $Log: not supported by cvs2svn $
% Revision 1.13  2002/11/12 01:33:18  arno
% add BDF read file command
%
% Revision 1.12  2002/09/09 02:00:02  arno
% menu has changed
%
% Revision 1.11  2002/09/05 01:34:13  arno
% menu editing
%
% Revision 1.10  2002/09/04 23:21:18  arno
% updating menu help
%
% Revision 1.9  2002/08/16 16:27:28  arno
% new version with HTML menu generation
%      
% Revision 1.8  2002/07/29 18:37:55  arno       
% updating the menus       
%       
% Revision 1.7  2002/07/29 17:06:21  arno       
% adding menus       
%       
% Revision 1.6  2002/04/11 03:36:59  arno       
% change load/save sets menus       
%       
% Revision 1.5  2002/04/06 01:57:38  arno       
% title editing       
%       
% Revision 1.4  2002/04/06 01:20:04  arno       
% more editing title mainly      
%       
% Revision 1.3  2002/04/06 01:09:19  arno       
% increasing font size       
%       
% Revision 1.2  2002/04/06 00:38:13  arno       
% details on fonts size and formating      
%       
% Revision 1.1  2002/04/05 17:46:04  jorn       
% Initial revision       
%       
% 01-25-02 reformated help & license -ad       
% 03-13-02 updated event function calls -ad       
% 04-01-02 complete remodelling -ad       
       
function eeg_helpmenu( filename );       
       
allmenus = { ...       
'File' ''      ...
'   Import data' ''  ...    
'      From ASCII/float file or Matlab array' 'pop_editset'        ... %       Import Matlab data array      
'      From EGI .RAW file' 'pop_readegi'        ... %       Read .BCI datafile (Snapmaster)      
'      From BCI2000 ASCII file' 'pop_loadbci'        ... %       Read .BCI datafile (Snapmaster)      
'      From Snapmaster .SMA file' 'pop_snapread'       ... %       Read .SMA datafile (Snapmaster)      
'      From .BDF and Biosemi .EDF file' 'pop_readedf'        ... %       Read .ELP datafile (Neuroscan continuous)      
'      From Neuroscan .CNT file' 'pop_loadcnt'        ... %       Read .CNT datafile (Neuroscan continuous)      
'      From Neuroscan .EEG file' 'pop_loadeeg'        ... %       Read .EEG datafile (Neuroscan epochs)      
'   Import epoch info' ''      ...
'      From Matlab array or ASCII file' 'pop_importepoch'    ... %       Import Matlab array or ASCII file      
'      From Neuroscan .DAT file' 'pop_loaddat'        ... %       Import .DAT info file (Neuroscan epochs)      
'   Import event info' ''      ...
'      From  Matlab array or ASCII file' 'pop_importevent'    ... %       Import Matlab array or ASCII file      
'      From data channel' 'pop_chanevent'      ... %       Import events from channels      
'      From Presentation .LOG file' 'pop_importpres'     ... %       Import .LOG event file (Presentation)      
'   Load existing dataset' 'pop_loadset'        ... %    Load dataset      
'   Save current dataset' 'pop_saveset'        ... %    Save dataset      
'   Save datasets' 'pop_saveset'        ... %    Save dataset      
'   Clear dataset(s)' 'pop_delset'         ... %    Clear dataset(s)      
'   Maximize memory' 'pop_editoptions'    ... %    Maximize memory      
'   Save history' 'pop_saveh'          ... %    Save history      
'   Quit' ''      ...
'Edit' ''      ...
'   Dataset info' 'pop_editset'        ... %    Edit dataset info      
'   Event fields' 'pop_editeventfield' ... %    Edit event fields      
'   Event values' 'pop_editeventvals'  ... %    Edit event values      
'   About this dataset' 'pop_comments'        ... %    About this dataset      
'   Channel locations' 'pop_chanedit'       ... %    Edit channels      
'   Select data' 'pop_select'         ... %    Select data      
'   Select events' 'pop_selectevent'    ... %    Select events      
'   Copy current dataset' 'pop_copyset'        ... %    Copy current dataset      
'   Append datasets' 'pop_mergeset'       ... %    Append another dataset      
'   Delete dataset(s)' 'pop_delset'         ... %    Delete dataset(s)      
'Tools' ''      ...
'   Change sampling rate' 'pop_resample'       ... %    Change sampling rate      
'   Filter the data' 'pop_eegfilt'        ... %    Filter the data      
'   Average reference' 'pop_averef'         ... %    Average reference      
'   Reject continuous data' 'pop_eegplot'        ... %    Reject continuous data      
'   Extract epochs' 'pop_epoch'          ... %    Extract epochs      
'   Remove baseline' 'pop_rmbase'         ... %    Remove baseline      
'   Reject data epochs' ''  ...    
'      Reject data (all methods)' 'pop_rejmenu'        ... %       Reject data (all methods)      
'      Reject by inspection' 'pop_eegplot'        ... %       Reject by inspection      
'      Reject extreme values' 'pop_eegthresh'      ... %       Reject extreme values      
'      Reject flat line data' 'pop_rejtrend'       ... %       Reject flat line data      
'      Reject by probability' 'pop_jointprob'      ... %       Reject by probability      
'      Reject by kurtosis' 'pop_rejkurt'        ... %        Reject by kurtosis      
'      Reject by spectra' 'pop_rejspec'        ... %       Reject by spectra      
'      Reject marked epochs' 'pop_rejepoch'      ... %       Reject labeled epochs      
'   Run ICA' 'pop_runica'         ... %    Run ICA      
'   Remove components' 'pop_subcomp'        ... %    Remove components      
'   Reject using ICA' '' ...     
'      Reject components by map' 'pop_selectcomps'    ... %       Reject components by map      
'      Reject data (all methods)' 'pop_rejmenu'        ... %       Reject data (all methods)      
'      Reject by inspection' 'pop_eegplot'        ... %       Reject by inspection      
'      Reject extreme values' 'pop_eegthresh'      ... %       Reject extreme values      
'      Reject flat line activity' 'pop_rejtrend'       ... %       Reject flat line activity      
'      Reject by probability' 'pop_jointprob'      ... %       Reject by probability      
'      Reject by kurtosis' 'pop_rejkurt'        ... %       Reject by kurtosis      
'      Reject by spectra' 'pop_rejspec'        ... %       Reject by spectra      
'      Reject marked epochs' 'pop_rejepoch'      ... %       Reject labeled epochs      
'Plot' ''      ...
'   Channel locations' '' ...     
'      By name' 'topoplot'           ... %       Electrodes      
'      By number' 'topoplot'           ... %       Numbers      
'   EEG data (scroll)' 'pop_eegplot'        ... %    Scroll EEG data      
'   Channel spectra and maps' 'pop_spectopo'       ... %    Channel spectra and maps      
'   Channel properties' 'pop_prop'       ... %    Component properties      
'   Channel ERP image' 'pop_erpimage'       ... %    Channel ERP image      
'   Channel ERPs' ''      ...
'      With scalp maps' 'pop_timtopo'        ... %       ERP and scalp maps      
'      In scalp array' 'pop_plottopo'       ... %       ERP in scalp array      
'      In rect. array' 'pop_plotdata'       ... %       ERP in rect. array      
'   ERP maps' ''      ...
'      As 2-D scalp maps' 'pop_topoplot'       ... %       ERP scalp maps      
'      As 3-D head plots' 'pop_headplot'       ... %       ERP head plots      
'   Compare ERPs' 'pop_compareerps'    ... %    Compare ERPs      
'   Component activations (scroll)' 'pop_eegplot'        ... %    Scroll component activations      
'   Component spectra and maps' 'pop_spectopo'       ... %    Channel spectra and maps with components      
'   Component maps' ''      ...
'      As 2-D scalp maps' 'pop_topoplot'       ... %       Component scalp maps      
'      As 3-D head plots' 'pop_headplot'       ... %       Component head plots      
'   Component properties' 'pop_prop'       ... %    Component properties      
'   Component ERP image' 'pop_erpimage'       ... %    Component ERP image      
'   Component ERPs' ''      ...
'      With component maps)' 'pop_envtopo'        ... %       Largest ERP components      
'      In rectangular array' 'pop_plotdata'       ... %       Comp. ERP time courses      
'   Data statistics' ''      ...
'      Channel statistics' 'pop_signalstat'      ...
'      Component statistics' 'pop_signalstat'      ...
'      Event statistics' 'pop_eventstat'      ...
'   Time-frequency transforms' ''  ...
'      Channel time-frequency' 'pop_timef'          ... %       Channel time-frequency      
'      Channel cross-coherence' 'pop_crossf'         ... %       Channel cross-coherence      
'      Component time-frequency' 'pop_timef'          ... %       Component time-frequency EEGLAB'));     
'      Component cross-coherence' 'pop_crossf'         ... %       Component cross-coherence      
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
	command = allmenus( 2:2:length(allmenus));
	for index = 1:length(text)
		if isempty(allmenus{index*2})
			tmpparam{index} = {  '' allmenus{(index-1)*2+1} };
		else
			tmpparam{index} = {  [ allmenus{index*2} '.m' ] allmenus{(index-1)*2+1} };
		end;
	end;
 	makehtml( tmpparam, '/home/www/eeglab/allfunctions', 'outputfile', 'tempmenu.html', 'mainonly', 'on');
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
