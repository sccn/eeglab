% pop_eegplot() - Visually inspect EEG data in a scrolling display.
%                 Perform rejection or marking for rejection of manually 
%                 (and/or previously manually) selected data portions 
%                 (= stretches of continuous data or whole data epochs).
% Usage:
%   >> pop_eegplot( EEG ) % Scroll EEG channel data. Allow marking for rejection via
%                         % button 'Update Marks' but perform no actual data rejection.
%                         % Do not show or use marks from previous visual inspections
%                         % or from semi-auotmatic rejection.
%   >> pop_eegplot( EEG, typerej, superpose, reject );
%
% Graphic interface:
%   "Add to previously marked rejections" - [edit box] Can be either YES or NO. 
%                    Command line eqivalent: 'superpose'.
%   "Reject marked trials" - [edit box] Can be either YES or NO. Command line
%                    equivalent 'reject'.
% Inputs:
%   EEG        - input EEG dataset
%   typerej    - type of rejection 0 = EEG independent components; 
%                                  1 = EEG data channels. {default: 1, data channels}
%   superpose  - 0 = Show new marks only: Do not color the background of data portions 
%                    previously marked for rejection by visual inspection. Mark data 
%                    portions for rejection by first coloring them (by dragging the left 
%                    mouse button) and then pressing the 'Update Marks' or 'Reject' 
%                    button (see 'reject' below). Previous markings from visual inspection 
%                    will be lost.
%                1 = Show data portions previously marked by visual inspection plus 
%                    data portions selected in this window for rejection (by dragging 
%                    the left mouse button in this window). These are differentiated 
%                    using a lighter and darker hue, respectively). Pressing the 
%                    'Update Marks' or 'Reject' button (see 'reject' below)
%                    will then mark or reject all the colored data portions.
%                     {default: 0, show and act on new marks only}
%   reject     - 0 = Mark for rejection. Mark data portions by dragging the left mouse 
%                    button on the data windows (producing a background coloring indicating 
%                    the extent of the marked data portion).  Then press the screen button 
%                    'Update Marks' to store the data portions marked for rejection 
%                    (stretches of continuous data or whole data epochs). No 'Reject' button 
%                    is present, so data marked for rejection cannot be actually rejected 
%                    from this eegplot() window. 
%                1 = Reject marked trials. After inspecting/selecting data portions for
%                    rejection, press button 'Reject' to reject (remove) them from the EEG 
%                    dataset (i.e., those portions plottted on a colored background. 
%                    {default: 0, mark for rejection only}
% Outputs:
%   Modifications are applied to the current EEG dataset at the end of the
%   eegplot() call, when the user presses the 'Update Marks' or 'Reject' button.
%   NOTE: The modifications made are not saved into EEGLAB history.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), eegplot(), pop_rejepoch()

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
% Revision 1.23  2003/02/23 08:28:21  scott
% header edit -sm
%
% Revision 1.22  2003/02/18 22:17:52  arno
% typo in header
%
% Revision 1.21  2003/02/18 22:08:54  arno
% update header for GUI
%
% Revision 1.20  2003/02/03 18:36:53  scott
% header msg clarified -sm
%
% Revision 1.19  2003/01/10 01:11:49  arno
% removing default position
%
% Revision 1.18  2002/11/15 19:52:53  arno
% do not plot channel names for ICA continuous data
%
% Revision 1.17  2002/11/12 19:20:26  arno
% warning for rejecting continuous data
%
% Revision 1.16  2002/08/24 02:04:06  scott
% help msg
%
% Revision 1.15  2002/08/12 21:54:13  arno
% text
%
% Revision 1.14  2002/08/12 16:26:20  arno
% inputdlg2
%
% Revision 1.13  2002/08/08 01:33:56  arno
% adding history for continuous rejection
%
% Revision 1.12  2002/08/07 22:41:01  arno
% editing text
%
% Revision 1.11  2002/07/31 17:16:15  arno
% debugging
%
% Revision 1.10  2002/07/31 16:32:09  arno
% debugging
%
% Revision 1.9  2002/07/30 22:21:11  arno
% debugging
%
% Revision 1.8  2002/07/30 22:19:24  arno
% removing multieegplot.m
%
% Revision 1.7  2002/07/27 01:35:55  arno
% editing header
%
% Revision 1.6  2002/07/26 16:53:22  arno
% switching icacomp
%
% Revision 1.5  2002/07/08 22:03:02  arno
% nothing
%
% Revision 1.4  2002/06/25 22:23:56  scott
% *** empty log message ***
%
% Revision 1.3  2002/04/26 21:28:17  arno
% eeg_updatemenu
%
% Revision 1.2  2002/04/26 21:25:03  arno
% update call to eeg_store
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-27-02 added event latency recalculation for continuous data -ad

function [com] = pop_eegplot( EEG, icacomp, superpose, reject, topcommand);

com = '';
if nargin < 1
	help pop_eegplot;
	return;
end;	
if nargin < 2
	icacomp = 1;
end;	
if nargin < 3
	superpose = 0;
end;
if nargin < 4
	reject = 1;
end;
if icacomp == 0
	if isempty( EEG.icasphere )
		disp('Error: you must run ICA first'); return;
	end;
end;

if nargin < 3 & EEG.trials > 1

	% which set to save
	% -----------------
    promptstr    = { 'Add to previously marked rejections (yes/no)', ...
         	         'Reject marked trials (yes/no)', ...
						 };
	inistr       = { 'yes', 'no' };
	result       = inputdlg2( promptstr, fastif(icacomp==0, 'Manual component rejection -- pop_eegplot()', ...
								'Reject epochs by visual inspection -- pop_eegplot()'), 1,  inistr, 'pop_eegplot');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
   
   switch lower(result{1}), case 'yes', superpose=1; end;
   switch lower(result{2}), case 'no', reject=0; end;

end;

if EEG.trials > 1
    if icacomp == 1 macrorej  = 'EEG.reject.rejmanual';
        			macrorejE = 'EEG.reject.rejmanualE';
    else			macrorej  = 'EEG.reject.icarejmanual';
        			macrorejE = 'EEG.reject.icarejmanualE';
    end;
    elecrange = [1:EEG.nbchan];
	colrej = EEG.reject.rejmanualcol;
	rej  = eval(macrorej);
	rejE = eval(macrorejE);
	
	eeg_rejmacro; % script macro for generating command and old rejection arrays

else % case of a single trial (continuous data)
	     %if icacomp, 
    	 %    	command = ['if isempty(EEG.event) EEG.event = [eegplot2event(TMPREJ, -1)];' ...
         %         'else EEG.event = [EEG.event(find(EEG.event(:,1) ~= -1),:); eegplot2event(TMPREJ, -1, [], [0.8 1 0.8])];' ...
         %         'end;']; 
      	 %else, command = ['if isempty(EEG.event) EEG.event = [eegplot2event(TMPREJ, -1)];' ...
         %         'else EEG.event = [EEG.event(find(EEG.event(:,1) ~= -2),:); eegplot2event(TMPREJ, -1, [], [0.8 1 0.8])];' ...
         %         'end;']; 
      	 %end;
         %if reject
         %   command = ...
         %   [  command ...
         %      '[EEG.data EEG.xmax] = eegrej(EEG.data, EEG.event(find(EEG.event(:,1) < 0),3:end), EEG.xmax-EEG.xmin);' ...
         %      'EEG.xmax = EEG.xmax+EEG.xmin;' ...
         %   	'EEG.event = EEG.event(find(EEG.event(:,1) >= 0),:);' ...
         %      'EEG.icaact = [];' ...
         %      'EEG = eeg_checkset(EEG);' ];
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
    if reject == 0, command = [];
    else, 
        command = ...
            [  '[EEGTMP LASTCOM] = eeg_eegrej(EEG,eegplot2event(TMPREJ, -1));' ...
               'if ~isempty(LASTCOM),' ... 
               '  h(LASTCOM);' ...
               '  [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
               '  if ~isempty(LASTCOM),' ... 
               '     h(''[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);'');' ...
               '     eeglab(''redraw'');' ...
               '  end;' ...
               'end;' ...
               'clear EEGTMP;' ];
        res = questdlg2( strvcat('Mark stretches of continuous data for rejection', ...
                             'by dragging the left mouse button. Click on marked', ...
                             'stretches to unmark. When done,press "REJECT" to', ...
                             'excise marked stretches (Note: Leaves rejection', ...
                             'boundary markers in the event table).'), 'Warning', 'Cancel', 'Continue', 'Continue');
        if strcmpi(res, 'Cancel'), return; end;
    end; 
    eegplotoptions = { 'winlength', 5 };
    if ~isempty(EEG.chanlocs) & icacomp
        eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs };
    end;
end;

if icacomp == 1
	eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Scroll channel activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	%eeg_multieegplot( EEG.data, [], [], oldrej, oldrejE, 'title', 'Scroll channel activities -- eegplot()', 'srate', ...
	%	      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
else
	eeg_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
	    tmpdata = EEG.icaact;
	else
        tmpdata = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        tmpdata = reshape( tmpdata, size(tmpdata,1), EEG.pnts, EEG.trials);
    end;
	eegplot( tmpdata, 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	%eeg_multieegplot( tmpdata, [], [], oldrej, oldrejE, 'title', 'Scroll component activities -- eegplot()', 'srate', ...
	%	      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
end;

com = [ com sprintf('pop_eegplot( %s, %d, %d, %d);', inputname(1), icacomp, superpose, reject) ]; 
return;
