% eeg_helpmenu() - Call the help file for EEGLAB menus
%
% Usage: >> eeg_helpmenu( );
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.6  2002/04/11 03:36:59  arno
% change load/save sets menus
%
% Revision 1.5  2002/04/06 01:57:38  arno
% title editing
%
% Revision 1.4  2002/04/06 01:20:04  arno
% more editing, title mainly
%
% Revision 1.3  2002/04/06 01:09:19  arno
% increasing font size
%
% Revision 1.2  2002/04/06 00:38:13  arno
% details on fonts, size and formating
%
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 
% 03-13-02 updated event function calls -ad
% 04-01-02 complete remodelling -ad

function eeg_helpmenu();

commands = { ...
       ''                   ... % File                                           
       ''                   ... %    Import data                                 
       'pop_editset'        ... %       Import Matlab data array                 
       'pop_loadbci'        ... %       Read .BCI datafile (Snapmaster)          
       'pop_snapread'       ... %       Read .SMA datafile (Snapmaster)          
       'pop_loadelp'        ... %       Read .ELP datafile (Neuroscan continuous)
       'pop_loadcnt'        ... %       Read .CNT datafile (Neuroscan continuous)
       'pop_loadeeg'        ... %       Read .EEG datafile (Neuroscan epochs)    
       ''                   ... %    Import epoch info                           
       'pop_importepoch'    ... %       Import Matlab array or ASCII file        
       'pop_loaddat'        ... %       Import .DAT info file (Neuroscan epochs) 
       ''                   ... %    Import event info                           
       'pop_importevent'    ... %       Import Matlab array or ASCII file        
       'pop_importpres'     ... %       Import .LOG event file (Presentation)    
       'pop_loadset'        ... %    Load dataset                                
       'pop_saveset'        ... %    Save dataset                                
       'pop_saveset'        ... %    Save dataset                                
       'pop_delset'         ... %    Clear dataset(s)                            
       'pop_editoptions'    ... %    Maximize memory                             
       'pop_saveh'          ... %    Save history                                
       ''                   ... %    Quit                                        
       ''                   ... % Edit                                           
       'pop_editset'        ... %    Edit dataset info                           
       'pop_editeventfield' ... %    Edit event fields                           
       'pop_editeventvals'  ... %    Edit event values                           
       'pop_comment'        ... %    About this dataset                          
       'pop_select'         ... %    Select data                                 
       'pop_selectevent'    ... %    Select events                               
       'pop_copyset'        ... %    Copy current dataset                        
       'pop_mergeset'       ... %    Append another dataset                      
       'pop_delset'         ... %    Delete dataset(s)                           
       ''                   ... % Tools                                          
       'pop_resample'       ... %    Change sampling rate                        
       'pop_eegfilt'        ... %    Filter the data                             
       'pop_averef'         ... %    Average reference                           
       'pop_eegplot'        ... %    Reject continuous data                      
       'pop_epoch'          ... %    Extract epochs                              
       'pop_rmbase'         ... %    Remove baseline                             
       ''                   ... %    Reject data epochs                          
       'pop_rejmenu'        ... %       Reject data (all methods)                
       'pop_eegplot'        ... %       Reject by inspection                     
       'pop_eegthresh'      ... %       Reject extreme values                    
       'pop_rejtrend'       ... %       Reject flat line data                    
       'pop_jointprob'      ... %       Reject by probability                    
       'pop_rejkurt'        ... %        Reject by kurtosis                       
       'pop_rejspec'        ... %       Reject by spectra                        
       'pop_rejepochs'      ... %       Reject labeled epochs                    
       'pop_runica'         ... %    Run ICA                                     
       'pop_subcomp'        ... %    Remove components                           
       ''                   ... %    Reject using ICA                            
       'pop_selectcomps'    ... %       Reject components by map                 
       'pop_rejmenu'        ... %       Reject data (all methods)                
       'pop_eegplot'        ... %       Reject by inspection                     
       'pop_eegthresh'      ... %       Reject extreme values                    
       'pop_rejtrend'       ... %       Reject flat line activity                
       'pop_jointprob'      ... %       Reject by probability                    
       'pop_rejkurt'        ... %       Reject by kurtosis                       
       'pop_rejspec'        ... %       Reject by spectra                        
       'pop_rejepochs'      ... %       Reject labeled epochs                    
       ''                   ... % Plot                                           
       ''                   ... %    Plot channel locations                      
       'topoplot'           ... %       Electrodes                   
       'topoplot'           ... %       Numbers                      
       'pop_eegplot'        ... %    Scroll EEG data                             
       'pop_spectopo'       ... %    Channel spectra and maps                    
       'pop_erpimage'       ... %    Channel ERP image                           
       ''                   ... %    ERP plots                                   
       'pop_timtopo'        ... %       ERP and scalp maps                       
       'pop_plottopo'       ... %       ERP in scalp array                       
       'pop_plotdata'       ... %       ERP in rect. array                       
       ''                   ... %    Scalp plots                                 
       'pop_topoplot'       ... %       ERP scalp maps                           
       'pop_headplot'       ... %       ERP head plots                           
       'pop_compareerps'    ... %    Compare ERPs                                
       'pop_eegplot'        ... %    Scroll component activations                
       ''                   ... %    Component scalp plots                       
       'pop_topoplot'       ... %       Component scalp maps                     
       'pop_headplot'       ... %       Component head plots                     
       'pop_compprop'       ... %    Component properties                        
       'pop_erpimage'       ... %    Component ERP image                         
       ''                   ... %    Comp. ERP                      
       'pop_envtopo'        ... %       Largest ERP components                      
       'pop_plotdata'       ... %       Comp. ERP time courses                      
       ''                   ... %    Time-freq. plots                            
       'pop_timef'          ... %       Channel time-frequency                   
       'pop_crossf'         ... %       Channel cross-coherence                  
       'pop_timef'          ... %       Component time-frequency                 
       'pop_crossf'         ... %       Component cross-coherence                
 };

[textmenu nblines l] = getallmenus(findobj('tag', 'EEGLAB'));
fontsize(1:length( textmenu )) = { 15 };
fontweight(1:length( textmenu )) = { 'normal' };
for index = 1:length( textmenu )
    if textmenu(index,1) ~= 32, fontsize{ index } = fontsize{ index }+1; fontweight{ index } ='bold'; end;
    if textmenu(index,6) ~= 32, fontweight{ index } ='bold'; end;
    a = deblank(textmenu(index, :));
    if index > length( commands )
        commands = { commands{:} [] };
    end;
    if ~isempty(commands{index})
        a = [ a ' -- ' commands{index} '()'];
        textmenu(index, 1:length(a)) = a;
        commands{index} = [ 'pophelp(''' commands{index} ''');' ];
    end;
end;

textmenu = strvcat('Functions called through the EEGLAB menu', ...
		   '(For help message, click on blue text below)', ' ', ...
		   textmenu(1:end-1,:));
fontsize   = { 16 15 15 fontsize{:} };
fontweight = { 'bold' 'normal' 'normal' fontweight{:} };
commands   = { '' ' ' '' commands{:} };
textgui(textmenu(1:end-1,:), commands, ...
            'fontsize', fontsize, 'fontweight', fontweight, 'fontname', 'times', 'linesperpage', 17 );
return;
