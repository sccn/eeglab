% pop_erpimage() - draw an ERP-image plot of a given EEG channel or independent
%                  component. Uses a pop-up window if less than three (or four 
%                  in one condition) input arguments are supplied. Calls erpimage(). 
%                  For futher details see >> help erpimage
% Usage:
%   >> pop_erpimage(EEG, typeplot);          % pop-up a data entry window
%   >> pop_erpimage(EEG, typeplot, lastcom); % pop-up a data entry window
%   >> pop_erpimage(EEG, typeplot, channel); % no pop-up window
%   >> pop_erpimage(EEG, typeplot, channel, projchan, title, ...
%                  smooth, decimate, sortingtype, sortingwin, ...
%                            sortingeventfield, renorm, options...);
% Graphic interface:
%   "Channel or Component" - [edit box] Enter channel number or component
%                 number to plot. erpimage() equivalent: 'channel' 
%   "Project to channel #" - [edit box] (for plotting independent components). 
%                 Allow reprojecting the component activity 
%                 to a given channel or group of channels. 
%                 erpimage() equivalent: 'projchan' 
%   "Smoothing" - [text box] Smoothing parameter in number of trials.
%                 erpimage() equivalent: 'smooth' 
%   "Downsampling" - [edit box] Decimate parameter. 
%                 erpimage() equivalent: 'decimate' 
%   "Time limits" - [edit box] Enter the time limits in milliseconds. 
%                 erpimage() equivalent: the 1st and 2nd parameters of the 'limit' array
%   "Figure title" - [edit box] Enter the figure title here.  If empty, a title
%                 is automatically generated. erpimage() equivalent: 'title' 
%   "Plot scalp map" - [checkbox] Setting this option plot a scalp map of the 
%                 channel location (or component topography) next to the 
%                 erpimage. erpimage() equivalent: 'topo' 
%   "plot ERP" - [checkbox] Setting this option plot the channel or component 
%                 ERP below the ERP image. erpimage() equivalent: 'erp' 
%   "Plot colorbar" - [checkbox] Plot the colorbar on the right of the erpimage. 
%                 erpimage() equivalent: 'cbar' 
%   "ERP limits" - [edit box] Set the minimum and maximum value for the ERP plot
%                 erpimage() equivalent: 3rd and 4th parameters of the 'limit' array 
%   "Color limits" - [edit box] Set the color limits for the ERP image. 
%                 erpimage() equivalent: 'caxis' 
%   "Epoch sorting field" - [button and edit box] Specify the event field which
%                 values will be used to sort the trials. For instance, if you
%                 select the 'latency' fields, trials will be sorted by the 
%                 latency of the selected events.
%                 erpimage() equivalent: 'sortingeventfield'
%   "Event type(s)" - [button and edit box] Specify which event subset to select,
%                 based on the event 'type' field values (to scan for event types, use 
%                 menu Edit > Events values and look at the values of the 'type'
%                 field). For instance, entering type 'rt' (if defined) and field 
%                 'latency' in option aboce will sort trials on reaction time latencies. 
%                 When several selected events are present in individual trials, 
%                 the first event values are used for sorting and a warning is issued. 
%                 erpimage() equivalent: 'sortingtype'
%   "Event time range" - [edit box] Specify which event subset to select based on 
%                 event latency values. As the option above, this further restrains 
%                 the selection of events. For example, entering [200 300] in this 
%                 box, 'rt' for the event type (above), and 'latency' for the 
%                 epoch sorting field will select trials with reaction-time latencies
%                 in between 200 and 300 ms. Trials with no such event will not be
%                 included in the ERP-image plot. erpimage() equivalent: 'sortingwin'  
%   "rescale" - [edit box] 'yes', 'no', or a Matlab formula. 
%                 erpimage() equivalent: 'renorm' 
%   "align" - [edit box] Set this to 'Inf' to re-align the individual trials 
%                 on the median latency of the selected events. Else, enter an epoch time 
%                 (in ms) to align the events to (Ex: 0). erpimage() equivalent: 'align' 
%   "Don't sort by value" - [checkbox] Check this box if you do not want to 
%                 sort the trials but do want to plot the selected event values. 
%                 erpimage() equivalent: 'nosort' 
%   "Don't plot value" - [checkbox] Check this box if you do not want to 
%                 plot the selected event values, but still want to sort 
%                 the data trials according to these values. 
%                 erpimage() eqivalent: 'noplot' 
%   "Sort by phase > Frequency" - [edit box] Specify the frequency or frequency 
%                 range for sorting trials by phase. erpimage() equivalent: 
%                 3rd and 4th inputs to 'phasesort' 
%   "Window center (ms)" - [edit box] erpimage() equivalent: 1st 'phasesort' input
%   "Percent low-amp. trials to ignore" - [edit box] erpimage() equivalent: 
%                 2nd 'phasesort' input 
%   "Wavelet cycle" - [text] Number of wavelet cycles used for spectral 
%                 decomposition at the specified latency. Cannot be edited.
%   "Inter-trial coherence options > Frequency" - [edit box] Frequency at which  
%                 to compute coherence. Constrained to be the same as the 
%                 "Sort by phase > Frequency" edit box. erpimage() equivalent: 'coher' 
%   "Signif. level" - [edit box] Coherence significance cutoff, as a proability
%                  (Ex: .05). erpimage() equivalent: 'signif' 
%   "Amplitude limit" - [edit box] Amplitude limits [min max] for the data power 
%                 plot at the selected frequency. erpimage() equivalent:
%                 5th and 6th inputs to 'limit' 
%   "Coher limits" - [edit box] Upper limit (<=1) for the coherence 
%                 plot. erpimage() equivalent: 7th and 8th inputs of 'limit' 
%   "Image amps" - [checkbox] Check this box for plotting the spectral amplitude
%                 image at the selected frequency (instead of plotting EEG potential). 
%                 erpimage() equivalent: 'plotamp'
%   "Plot spectrum" - [edit box] Plot the channel or component data spectrum in 
%                 the top right corner of the ERP image. erpimage() equivalent: 'spec' 
%   "Baseline ampl." - [edit box] Baseline amplitude for data power plot at the 
%                 selected frequency. erpimage() equivalent: 7th inputs of 'limit'
%   "Mark times" - [edit box] Time(s) in ms to plot vertical lines.
%                 erpimage() equivalent: 'vert' 
%   "More options" - [edit box] Enter 'key', 'value' sequences. Other erpimage()
%                 options not handled by this interface, including: 'erpstd' to 
%                 plot the ERP standard deviation; 'auxvar' to plot auxilary 
%                 variables; 'ampsort' to sort trials based on amplitude at 
%                 the selected frequency, etc.  For further information see  
%                 >> help erpimage() 
% Inputs:
%   EEG        - dataset structure
%   typeplot   - 1=channel, 0=component {default: 1}
%   lastcom    - string containing previous pop_erpimage command (from LASTCOM) 
%                or from the previous function call output.  The values from this 
%                function call are used as default in the graphic interface.
%
% Commandline options:
%   channel    - Index of channel or component(s) to plot {default: 1}
%   projchan   - Channel to back-project the selected component(s) to. 
%                If plotting channel activity, this argument is ignored. 
%                If [], the ICA component activation is plotted {default []}.
%   title      - ['string'] Plot title {default: []}
%   smooth     - Smoothing parameter (number of trials). {Default: 5} 
%                erpimage() equivalent: 'avewidth' 
%   decimate   - Decimate parameter (i.e. ratio of trials_in/trials_out).
%                erpaimge() equivalent: 'decimate' {Default: 0}
%   sortingtype - Sorting event type(s) ([int vector]; []=all). See Notes below.
%                Either a string or an integer.
%   sortingwin - Sorting event window [start, end] in seconds ([]=whole epoch)
%   sortingeventfield - Sorting field name. {default: none}. See Notes below.
%   options    - erpimage() options, separated by commas (Ex: 'erp', 'cbar'). 
%                {Default: none}. For further details see >> erpimage help   
%
% Outputs from pop-up: 
%   String containing the command used to evaluate this plotting function
%   (saved by eeglab() as LASTCOM). Enter it as the 'lastcom' input to restore
%   the previous parameters as defaults in a new pop_erpimage() pop-up window
%
% Outputs from commandline:
%   Same as erpimage(). Note: No outputs are returned when a window pops-up 
%   to ask for additional arguments
%   
% Notes:
%   1) A new figure is created only when the pop-up window is called, 
%      so you may call this command with >3 args to draw in sbplot() axes.
%   2) To sort epochs, first define the event field to be used with
%      the argument 'sortingeventfield' (for instance 'latency'). Then,
%      because there may be several events with different latencies in a
%      given epoch, it is possible to consider only a subsets of events
%      using the 'sortingtype' and 'sortingwin' arguments. The 
%      'sortingtype' argument selects events with definite types. The 
%      'sortingwin' argument helps to define a specific time window in the 
%      epoch to select events. For instance the epoch time range may be -1 
%      to 2 seconds but one may want to select events only in the range 0 
%      to 1 second. These three parameters are forwarded to the function
%      eeg_getepochevent(), whose help message contains more details.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), erpimage(), eeg_getepochevent()

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
% Revision 1.123  2004/09/21 17:46:42  hilit
% changed the if statement from if isempty(EEG.chanlocs) -> if ~
%
% isempty(EEG.chanlocs)
%
% Revision 1.122  2004/06/16 22:24:33  arno
% no scalp map if only channel labels
%
% Revision 1.121  2004/06/16 22:21:20  arno
% debug for numercial types
%
% Revision 1.120  2004/05/15 00:33:25  arno
% same
%
% Revision 1.119  2004/05/15 00:30:12  arno
% putting back microvolt
%
% Revision 1.118  2004/03/02 21:55:00  arno
% recover version 1.114
%
% Revision 1.114  2003/05/10 19:06:21  arno
% debug last
%
% Revision 1.113  2003/05/10 19:02:42  arno
% debug last
%
% Revision 1.112  2003/05/10 19:01:13  arno
% num2str -> vararg2str
%
% Revision 1.111  2003/05/10 18:49:18  arno
% same
%
% Revision 1.110  2003/05/10 18:47:23  arno
% condense output string
%
% Revision 1.109  2003/04/24 21:53:56  arno
% debuging renorm default option
%
% Revision 1.108  2003/04/23 23:14:53  arno
% cycle
%
% Revision 1.107  2003/04/23 01:09:34  arno
% change text for default wavelet cycle
%
% Revision 1.106  2003/03/12 20:11:59  scott
% header edits .m
% -sm
%
% Revision 1.105  2003/03/12 19:21:37  arno
% same
%
% Revision 1.104  2003/03/12 19:15:00  arno
% debug renorm
%
% Revision 1.103  2003/03/12 19:12:06  arno
% debuging history for rescale
%
% Revision 1.102  2003/03/12 03:13:12  arno
% help msg
%
% Revision 1.101  2003/03/12 02:43:51  arno
% two help buttons
%
% Revision 1.100  2003/02/25 00:59:31  scott
% edit header -sm
%
% Revision 1.99  2003/02/24 16:09:57  arno
% resolving ???
%
% Revision 1.98  2003/02/23 09:06:15  scott
% header edits -sm
%
% Revision 1.97  2003/02/18 19:21:40  scott
% header edits -sm
%
% Revision 1.96  2003/02/17 02:25:28  arno
% reformating text for new functionality in help2html
%
% Revision 1.95  2003/02/17 01:52:32  arno
% updating header
%
% Revision 1.94  2003/02/17 00:59:51  arno
% adding GUI info
%
% Revision 1.93  2003/01/02 16:52:55  arno
% removing "ERP image" from title, plotting component scalp map
%
% Revision 1.92  2002/11/19 22:55:02  arno
% debug for component plotting
%
% Revision 1.91  2002/11/18 18:21:21  arno
% nicer aspect ratio
%
% Revision 1.90  2002/11/13 18:30:04  scott
% same
%
% Revision 1.89  2002/11/13 18:28:16  scott
% freq range text
%
% Revision 1.88  2002/11/13 18:25:43  scott
% freq range legends
%
% Revision 1.87  2002/11/13 18:23:45  scott
% coher freq -> coher freq range
%
% Revision 1.86  2002/10/29 23:41:53  arno
% fixing colorbar problem and help button
%
% Revision 1.85  2002/10/29 23:33:27  arno
% same
%
% Revision 1.84  2002/10/29 23:32:44  arno
% debugging channel name when no channel location file
%
% Revision 1.83  2002/10/16 16:24:37  arno
% debug command line call
%
% Revision 1.82  2002/10/15 16:56:44  scott
% debug
%
% Revision 1.81  2002/10/15 16:56:07  scott
% debug
%
% Revision 1.80  2002/10/15 16:52:31  scott
% debug
%
% Revision 1.79  2002/10/15 16:42:27  scott
% debug
%
% Revision 1.78  2002/10/15 16:39:36  scott
% default title with channel name
%
% Revision 1.77  2002/10/15 15:08:17  scott
% *** empty log message ***
%
% Revision 1.76  2002/10/15 14:48:55  scott
% removed erpimopt, added 'Color limits' help
%
% Revision 1.75  2002/10/10 14:56:56  scott
% [] -> default in help msgs where applicable
%
% Revision 1.74  2002/09/06 18:24:39  arno
% reparing history
%
% Revision 1.73  2002/09/06 16:28:46  arno
% same
%
% Revision 1.72  2002/09/06 16:27:35  arno
% debug command
%
% Revision 1.71  2002/08/31 17:00:47  arno
% add yerplabel option
%
% Revision 1.70  2002/08/30 17:49:23  arno
% same
%
% Revision 1.69  2002/08/30 17:47:34  arno
% debug
%
% Revision 1.68  2002/08/29 01:20:55  arno
% nothing
%
% Revision 1.67  2002/08/23 20:17:20  scott
% help msg
%
% Revision 1.66  2002/08/23 18:31:26  arno
% new features, complete gui reprogrammation with tags, new components only options
%
% Revision 1.65  2002/08/20 23:55:27  arno
% updating text and tooltipstring
%
% Revision 1.64  2002/08/20 23:29:38  scott
% changed box labels again - but not the tooltipstring messages!
%
% Revision 1.63  2002/08/19 23:19:51  arno
% back to version 1.47
%
% Revision 1.47  2002/08/17 03:06:27  arno
% same
%
% Revision 1.46  2002/08/17 02:49:56  arno
% listdlg2
%
% Revision 1.45  2002/08/13 23:12:31  scott
% add erpaimge() to figure number
%
% Revision 1.44  2002/08/12 23:34:36  arno
% reordering fields
%
% Revision 1.43  2002/08/12 18:50:47  arno
% errordlg2
%
% Revision 1.42  2002/08/12 01:35:59  arno
% allamp debug
%
% Revision 1.41  2002/08/12 01:35:36  arno
% color
%
% Revision 1.40  2002/08/11 22:12:40  arno
% color
%
% Revision 1.39  2002/08/09 16:45:47  arno
% updating plotamp and phasesort
%
% Revision 1.38  2002/07/29 23:42:47  arno
% same
%
% Revision 1.37  2002/07/29 23:41:26  arno
% updating getdef
%
% Revision 1.36  2002/07/29 00:20:00  arno
% debugging
%
% Revision 1.35  2002/07/27 01:21:56  arno
% debugging last command
%
% Revision 1.34  2002/07/27 00:45:17  arno
% adding messages, changing output command
%
% Revision 1.33  2002/07/24 19:29:10  arno
% removing renorm
%
% Revision 1.32  2002/07/24 19:27:16  arno
% adding back renorm into text output
%
% Revision 1.31  2002/07/24 19:23:56  arno
% ebugging output message
%
% Revision 1.30  2002/07/24 19:22:05  arno
% end debugging
%
% Revision 1.29  2002/07/24 19:20:58  arno
% same
%
% Revision 1.28  2002/07/24 19:19:34  arno
% debugging
%
% Revision 1.27  2002/07/23 21:32:58  arno
% debugging not same size
%
% Revision 1.26  2002/05/01 23:43:05  arno
% additional test
%
% Revision 1.25  2002/04/29 20:22:22  arno
% typo
%
% Revision 1.24  2002/04/26 21:41:27  arno
% phase/coher consistency chack
%
% Revision 1.23  2002/04/25 18:18:24  arno
% debugging spec
%
% Revision 1.22  2002/04/25 18:13:11  arno
% updating command output
%
% Revision 1.21  2002/04/25 17:55:05  arno
% improving commandline message
% ,
%
% Revision 1.20  2002/04/25 17:52:03  arno
% modifying renorm
%
% Revision 1.19  2002/04/25 16:05:57  scott
% adding blank line above history entry -sm
%
% Revision 1.18  2002/04/25 16:03:12  scott
% trying to output com into history as well as popcom -sm
%
% Revision 1.17  2002/04/24 20:36:22  arno
% implementing buttons
%
% Revision 1.16  2002/04/24 19:50:28  arno
% updating parameters
%
% Revision 1.15  2002/04/24 15:22:28  scott
% [same] -sm
%
% Revision 1.14  2002/04/24 15:21:48  scott
% [same] -sm
%
% Revision 1.13  2002/04/24 15:20:47  scott
% [same] -sm
%
% Revision 1.12  2002/04/24 15:19:53  scott
% adding note re plotamps -sm
%
% Revision 1.11  2002/04/23 21:51:01  arno
% fixed vert
%
% Revision 1.10  2002/04/23 21:24:33  scott
% *** empty log message ***
%
% Revision 1.9  2002/04/23 21:20:15  scott
% edited help msg -sm
%
% Revision 1.8  2002/04/23 17:19:45  arno
% contextual help modif
%
% Revision 1.7  2002/04/23 02:01:32  arno
% New graphic interface
%
% Revision 1.6  2002/04/18 19:11:09  arno
% testing version control
%
% Revision 1.5  2002/04/18 15:56:33  scott
% editted msgs -sm
%
% Revision 1.4  2002/04/11 20:02:42  arno
% adding last command history
%
% Revision 1.3  2002/04/11 19:06:43  arno
% adding option to default parameters
%
% Revision 1.2  2002/04/05 23:05:50  arno
% Correct typo
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 02-12-02 added new event format compatibility -ad 
% 02-15-02 text interface editing -sm & ad 
% 03-07-02 add the eeglab computeica options -ad
% 02-15-02 modified the function accoring to the new event/epoch structure -ad
% 03-18-02 added title -ad & sm
% 04-04-02 added outputs -ad & sm

function varargout = pop_erpimage( EEG, typeplot, channel, projchan, titleplot, smooth, decimate, sortingtype, ...
            sortingwin, sortingeventfield, varargin)
 
varargout{1} = '';
if nargin < 1
   help pop_erpimage;
   return;
end;

if typeplot == 0 & isempty(EEG.icasphere)
   error('no ICA data for this set, first run ICA');
end;   
if EEG.trials == 1
   error('erpimage of one trial cannot be plotted');
end;   

if nargin < 2	
	typeplot = 1; %1=signal; 0=component
end;
	lastcom = [];
if nargin < 3
	popup = 1;
else
	popup = isstr(channel) | isempty(channel);
	if isstr(channel)
		lastcom = channel;
	end;
end;

if popup
	% get contextual help
	% -------------------
	[txt2 vars2] = gethelpvar('erpimage.m');
	txt  = { txt2{:}};
	vars = { vars2{:}};
	% [txt vars] = gethelpvar('erpimopt.m');
	% txt  = { txt{:} txt2{:}};
	% vars = { vars{:} vars2{:}};
	commandphase = [ 'if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''phase''),''string'')),' ...
					 '   if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''coher''),''string'')), ' ...
					 '      set(findobj(''parent'', gcbf, ''tag'', ''coher''), ''string'', get(findobj(''parent'', gcbf, ''tag'', ''phase''),''string''));' ...
					 'end; end;' ];
	commandcoher = [ 'if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''coher''),''string'')),' ...
					 '   if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''phase''),''string'')), ' ...
					 '      set(findobj(''parent'', gcbf, ''tag'', ''phase''), ''string'', get(findobj(''parent'', gcbf, ''tag'', ''coher''),''string''));' ...
					 'end; end;' ];
	
	commandfield = ['if isempty(EEG.event)' ...
				   '   errordlg2(''No events'');' ...
				   'else' ...
				   '   tmpfieldnames = fieldnames(EEG.event);' ...
				   '   [tmps,tmpv] = listdlg2(''PromptString'', ''Select fields'', ''SelectionMode'',''single'',''ListString'', tmpfieldnames);' ...
				   '   if tmpv' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''field''), ''string'', tmpfieldnames{tmps});' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmps tmpv tmpfieldnames;' ];
	commandtype = ['if ~isfield(EEG.event, ''type'')' ...
				   '   errordlg2(''No type field'');' ...
				   'else' ...
				   '   if isstr(EEG.event(1).type)' ...
				   '       tmpfieldnames = unique({EEG.event.type});' ...
				   '   else' ...
				   '       tmpstr = unique(cell2mat({EEG.event.type}));' ...
				   '       tmpfieldnames = cell(1, length(tmpstr));' ...
                   '       for tmps=1:length(tmpstr), tmpfieldnames{tmps} = num2str(tmpstr(tmps)); end;' ...
				   '   end;' ...
		           '   [tmps,tmpv, tmpstr] = listdlg2(''PromptString'',''Select fields'', ''ListString'', tmpfieldnames);' ...
				   '   if tmpv' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''type''), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmps tmpv tmpstr tmpfieldnames;' ];
	
	geometry = { [1 1 0.1 0.8 2.1] [1 1 1 1 1] [1 1 1 1 1] [1 1 1 1 1] [1] [1] [1 1 1 0.8 0.8 1.2] [1 1 1 0.8 0.8 1.2] [1] [1] ...
				 [1.6 1.7 1.2 1 .5] [1.6 1.7 1.2 1 .5] [1] [1] [1.5 1 1 1 1] [1.5 1 1 1 1] [1] [1] [1.5 1 1 2.2] [1.5 1 1 2.2]};
    uilist = { { 'Style', 'text', 'string', fastif(typeplot, 'Channel', 'Component(s)'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,3,[],'1'), 'tag', 'chan' } { } ...
			   { 'Style', 'text', 'string', 'Figure title', 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', '', 'tag', 'title'  } ...
			   ...
			   { 'Style', 'text', 'string', 'Smoothing', 'fontweight', 'bold', 'tooltipstring', context('avewidth',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 6, [], int2str(min(max(EEG.trials-5,0), 10))), 'tag', 'smooth' } ...
			   { 'Style', 'checkbox', 'string', 'Plot scalp map', 'tooltipstring', 'plot a 2-d head map (vector) at upper left', ...
				 'value', getkeyval(lastcom, 'topo', 'present', 1), 'tag', 'plotmap'  } { } { } ...
			   { 'Style', 'text', 'string', 'Downsampling', 'fontweight', 'bold', 'tooltipstring', context('decimate',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 7, [], '1'), 'tag', 'decimate' } ...
			   { 'Style', 'checkbox', 'string', 'Plot ERP', 'tooltipstring', context('erp',vars,txt), 'value', getkeyval(lastcom, '''erp''', 'present', 1), 'tag', 'erp' } ...
			   { 'Style', 'text', 'string', fastif(typeplot, 'ERP limits (uV)','ERP limits'), 'tooltipstring', [ 'Plotting limits for ERP trace [min_uV max_uV]' 10 '{Default: ERP data limits}']  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'limits', [3:4]), 'tag', 'limerp'  } ...
			   { 'Style', 'text', 'string', 'Time limits (ms)', 'fontweight', 'bold', 'tooltipstring',  'Select time subset in ms' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'limits', [1:2], num2str(1000*[EEG.xmin EEG.xmax])), 'tag', 'limtime' } ...
			   { 'Style', 'checkbox', 'string', 'Plot colorbar','tooltipstring', context('caxis',vars,txt), 'value', getkeyval(lastcom, 'cbar', 'present', 1), 'tag', 'cbar' } ...
			   { 'Style', 'text', 'string', 'Color limits (see Help)','tooltipstring', context('caxis',vars,txt) } ...
			   { 'Style', 'edit', 'string',  getkeyval(lastcom, 'caxis'), 'tag', 'caxis' } ...
			   {} ...
			   { 'Style', 'text', 'string', 'Sort/align trials by epoch event values', 'fontweight', 'bold'} ...
			   { 'Style', 'pushbutton', 'string', 'Epoch-sorting field', 'callback', commandfield, ...
				 'tooltipstring', 'Epoch-sorting event field name (Ex: latency; default: no sorting):' } ...
			   { 'Style', 'pushbutton', 'string', 'Event type(s)', 'callback', commandtype, 'tooltipstring', ['Event type(s) subset (default: all):' 10 ...
                '(See ''/Edit/Edit event values'' for event types)']  } ...
			   { 'Style', 'text', 'string', 'Event time range', 'tooltipstring', [ 'Sorting event window [start, end] in milliseconds (default: whole epoch):' 10 ...
												  'events are only selected within this time window (can be usefull if several' 10 ...
												  'events of the same type are in the same epoch, or for selecting trials with given response time)']} ...
			   { 'Style', 'text', 'string', 'Rescale', 'tooltipstring', 'Rescale sorting variable to plot window (yes|no|a*x+b)(Ex:3*x+2):' } ...
			   { 'Style', 'text', 'string', 'Align', 'tooltipstring',  context('align',vars,txt) } ...
			   { 'Style', 'checkbox', 'string', 'Don''t sort by value', 'tooltipstring', context('nosort',vars,txt), 'value', getkeyval(lastcom, 'nosort', 'present', 0), 'tag', 'nosort' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 10), 'tag', 'field' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 8), 'tag', 'type' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 9), 'tag', 'eventrange' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'renorm','', 'no'), 'tag', 'renorm' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'align'), 'tag', 'align' } ...
			   { 'Style', 'checkbox', 'string', 'Don''t plot values', 'tooltipstring', context('noplot',vars,txt), 'value', getkeyval(lastcom, 'noplot', 'present', 0), 'tag', 'noplot' } ...
			   {} ...
			   { 'Style', 'text', 'string', 'Sort trials by phase', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Frequency (Hz | minHz maxHz)', 'tooltipstring', ['sort by phase at maximum-power frequency' 10 ...
               'in the data within the range [minHz,maxHz]' 10 '(overrides frequency specified in ''coher'' flag)']  } ...
			   { 'Style', 'text', 'string', 'Percent low-amp. trials to ignore', 'tooltipstring', ['percent of trials to reject for low' ...
			    'amplitude. Else,' 10 'if prct is in the range [-100,0] -> percent to reject for high amplitude'] } ...
			   { 'Style', 'text', 'string', 'Window center (ms)', 'tooltipstring', 'Ending time fo the 3 cycle window'  } ...
			   { 'Style', 'text', 'string', 'Wavelet cycles', 'tooltipstring', 'cycles per wavelet window'  } {}...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'phasesort', [3:4]), 'tag', 'phase', 'callback', commandphase } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'phasesort', [2]), 'tag', 'phase2' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'phasesort', [1]), 'tag', 'phase3' } ...
			   { 'Style', 'text', 'string', '        3'  } {}...
			   {} ...
			   { 'Style', 'text', 'string', 'Inter-trial coherence options', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Frequency (Hz | minHz maxHz)', 'tooltipstring', [ '[freq] -> plot erp plus amp & coher at freq (Hz)' 10 ...
               '[minHz maxHz] -> find max in frequency range' 10 '(or at phase freq above, if specified)']} ...
			   { 'Style', 'text', 'string', 'Signif. level (<0.20)', 'tooltipstring', 'add coher. signif. level line at alpha (alpha range: (0,0.1])' } ...
			   { 'Style', 'text', 'string', 'Amplitude limits (dB)'  } ...
			   { 'Style', 'text', 'string', 'Coher limits (<=1)'  } ...
			   { 'Style', 'checkbox', 'string', 'Image amps', 'tooltipstring',  context('plotamps',vars,txt), 'value', getkeyval(lastcom, 'plotamps', 'present', 0), 'tag', 'plotamps' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'coher', [1:2]), 'tag', 'coher', 'callback', commandcoher } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'coher', [3]), 'tag', 'coher2'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'limits',[5:6]), 'tag', 'limamp' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'limits',[7:8]), 'tag', 'limcoher' } {'style', 'text', 'string', '   (Requires signif.)' } ... 
			   {} ...
			   { 'Style', 'text', 'string', 'Other options', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Plot spectrum (minHz maxHz)','tooltipstring',  context('spec',vars,txt)} ...
			   { 'Style', 'text', 'string', 'Baseline ampl. (dB)', 'tooltipstring', 'Use it to fix baseline amplitude' } ...
			   { 'Style', 'text', 'string', 'Mark times (ms)','tooltipstring',  context('vert',vars,txt)} ...
			   { 'Style', 'text', 'string', 'More options (see >> help erpimage)' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'spec'), 'tag', 'spec' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'limits',9), 'tag', 'limbaseamp' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'vert'), 'tag', 'vert' } ...
			   { 'Style', 'edit', 'string', '', 'tag', 'others' } ...
			};
	if typeplot == 0 % add extra param for components
		geometry = { [1 1 0.1 0.8 2.1] geometry{:} };
		uilist   = { { } { } { } { } { } uilist{:}};
		uilist{1} = uilist{6};
		uilist{2} = uilist{7};
		uilist{6} = { 'Style', 'text', 'string', 'Project to channel #', 'fontweight', 'bold','tooltipstring', ['Project component(s) to data channel' 10 ...
												  'This allow to plot component activity in microvolt'] };
		uilist{7} = { 'Style', 'edit', 'string', getkeyval(lastcom, 4), 'tag', 'projchan' };
	end;
    [oldres a b res] = inputgui( geometry, uilist, 'pophelp(''pop_erpimage'');', ...
							 fastif( typeplot, 'Channel ERP image -- pop_erpimage()', 'Component ERP image -- pop_erpimage()'));
	if isempty(oldres), return; end;

	% first rows
	% ---------
	limits(1:8)  = NaN;
	channel   	 = eval( [ '[' res.chan ']' ]);
	titleplot    = res.title;
	if isfield(res, 'projchan'), projchan = str2num(res.projchan); else, projchan = []; end;
    options = '';
	if ~typeplot & isempty(projchan)
		options = [options ',''yerplabel'',''''' ];
    else
   		options = [options ',''yerplabel'',''\muV''' ];
	end;
	if isempty(titleplot)
        if typeplot==1 % if channel plot
            if ~isempty(EEG.chanlocs) % if channel information exist
                  titleplot = [ EEG.chanlocs(channel).labels ];
            else, titleplot = [ int2str(channel) ];
            end
        else
            titleplot = [ 'Comp. ' int2str(channel) ];
            if ~isempty(projchan),
                if ~isempty(EEG.chanlocs) % if channel plot
                      titleplot = [ titleplot ' -> ' EEG.chanlocs(projchan).labels ];
                else, titleplot = [ titleplot ' -> Chan. ' int2str(projchan) ];
                end
            end;
        end
    end;
	smooth       = eval(res.smooth);
    if res.plotmap
		if isfield(EEG.chanlocs, 'theta')
            if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end;
			if typeplot == 0
				     options = [options ',''topo'', { EEG.icawinv(:,' int2str(channel) ') EEG.chanlocs EEG.chaninfo } '];
			else     options = [options ',''topo'', { ' int2str(channel) ' EEG.chanlocs EEG.chaninfo } '];
			end;	
		end;
	end;
	
	decimate     = eval( res.decimate );
	if res.erp
		options = [options ',''erp'''];
	end;
	
	% finding limits
	% --------------
	if ~isempty(res.limerp)
		limits(3:4) = eval( [ '[' res.limerp ']' ]); 
	end;
	if ~isempty(res.limtime) % time limits
		if ~strcmp(res.limtime, num2str(1000*[EEG.xmin EEG.xmax]))
			limits(1:2) = eval( [ '[' res.limtime ']' ]);
		end;
	end;
	if ~isempty(res.limamp)
		limits(5:6) = eval( [ '[' res.limamp ']' ]);
	end;
	if ~isempty(res.limcoher)
		limits(7:8) = eval( [ '[' res.limcoher ']' ]);
	end;
	if ~isempty(res.limbaseamp)
		limits(9) = eval( res.limbaseamp ); %bamp
	end;
	if ~all(isnan(limits))
		options = [ options ',''limits'',' vararg2str(limits) '' ];
	end;
	
	% color limits
	% --------------
	if res.cbar
		options = [options ',''cbar'''];
	end;
	if res.caxis
		options = [options ',''caxis'',  [' res.caxis ']' ];
	end;
	
	% event rows
	% ----------
	if res.nosort
		options = [options ',''nosort'''];
	end;
	try, sortingeventfield = eval( res.field ); catch, sortingeventfield = res.field; end;
	sortingtype  = parsetxt(res.type);
	sortingwin   = eval( [ '[' res.eventrange ']' ] );
	if ~isempty(res.field) & ~strcmp(res.renorm, 'no')
		options = [options ',''renorm'', ''' res.renorm '''' ];
	end;
	if ~isempty(res.align)
		options = [options ',''align'', ' res.align ];
	end;
	if res.noplot
		options = [options ',''noplot'''];
	end;

	% phase rows
	% ----------
	tmpphase = [];
	if ~isempty(res.phase)
		tmpphase = eval( [ '[ 0 0 ' res.phase ']' ]);
	end;
	if ~isempty(res.phase2)
		tmpphase(2) = eval( res.phase2 );
	end;
	if ~isempty(res.phase3)
		tmpphase(1) = eval( res.phase3 );
	end;
	if ~isempty(tmpphase)
		options = [ options ',''phasesort'',' vararg2str(tmpphase) ];
	end;
	
	% coher row
	% ----------
	tmpcoher = [];
	if res.plotamps
		options = [options ',''plotamps'''];
	end;
	if ~isempty(res.coher)
		tmpcoher = eval( [ '[' res.coher ']' ]);
	end;
	if ~isempty(res.coher2)
		if length(tmpcoher) == 1
			tmpcoher(2) = tmpcoher(1);
		end;
		tmpcoher(3) = eval( res.coher2 );
	end;
	if ~isempty(tmpcoher)
		options = [ options ',''coher'',' vararg2str(tmpcoher)  ];
	end;

	% options row
	% ------------
	if ~isempty(res.spec)
		options = [options ',''spec'', [' res.spec ']' ];
	end;
	if ~isempty(res.vert)
		options = [options ',''vert'', [' res.vert ']' ];
	end;
	if ~isempty(res.others)
		options = [ options ',' res.others ];
	end;
	figure;
else
	options = '';
	if nargin < 4
		projchan = [];
	end;
	if nargin < 5
		titleplot = ' ';
	end;
	if nargin < 6
		smooth = 5;
	end;
	if nargin < 7
		decimate = 0;
	end;
	if nargin < 8
		sortingtype = [];
	end;
	if nargin < 9
		sortingwin = [];
	end;
	if nargin < 10
		sortingeventfield = [];
	end;
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else  
		  if ~iscell( varargin{ i } )
		      options = [ options ', [' num2str(varargin{i}) ']' ];
		  else
		      options = [ options ', { [' num2str(varargin{ i }{1}') ']'' EEG.chanlocs }' ];
		  end;    
		end;
	end;	
end;
try, icadefs; set(gcf, 'color', BACKCOLOR,'Name',' erpimage()'); catch, end;

% testing inputs
% --------------
if typeplot == 0 & length(channel) > 1 & isempty(projchan)
	error('A channel must be selected for plotting several components');
end;

% find sorting latencies
% ---------------------
typetxt = '';
if ~isempty(sortingeventfield)
    %events = eeg_getepochevent( EEG, sortingtype, sortingwin, sortingeventfield);
	events = sprintf('eeg_getepochevent( EEG, %s)', vararg2str({sortingtype, sortingwin, sortingeventfield}));
	
    % generate text for the command
    % -----------------------------
    for index = 1:length(sortingtype)
        if isstr(sortingtype{index})
            typetxt = [typetxt ' ''' sortingtype{index} '''' ];
        else
            typetxt = [typetxt ' ' num2str(sortingtype{index}) ];
        end;
    end;    
% $$$ 	% renormalize latencies if necessary
% $$$ 	% ----------------------------------
% $$$ 	switch lower(renorm)
% $$$ 	    case 'yes',
% $$$ 	         disp('Pop_erpimage warning: *** sorting variable renormalized ***');
% $$$ 	         events = (events-min(events)) / (max(events) - min(events)) * ...
% $$$ 	                        0.5 * (EEG.xmax*1000 - EEG.xmin*1000) + EEG.xmin*1000 + 0.4*(EEG.xmax*1000 - EEG.xmin*1000);
% $$$ 	    case 'no',;
% $$$ 	    otherwise,
% $$$ 	        locx = findstr('x', lower(renorm))
% $$$ 	        if length(locx) ~= 1, error('Pop_erpimage error: unrecognize renormalazing formula'); end;
% $$$ 	        eval( [ 'events =' renorm(1:locx-1) 'events' renorm(locx+1:end) ';'] );
% $$$ 	end;
else
	events = 'ones(1, EEG.trials)*EEG.xmax*1000';
    %events = ones(1, EEG.trials)*EEG.xmax*1000;
    sortingeventfield = '';
end;           

if typeplot == 1
	tmpsig = ['EEG.data(' int2str(channel) ', :)'];
else
    % test if ICA was computed or if one has to compute on line
    % ---------------------------------------------------------
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
		tmpsig = ['EEG.icaact([' int2str(channel) '], :)'];
	else
		tmpsig = ['EEG.icaact([' int2str(channel) '], :)'];
        tmpsig = ['EEG.icaweights([' int2str(channel) '],:)*EEG.icasphere*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts)'];
    end;
	if ~isempty(projchan)
		tmpsig = [ 'EEG.icawinv(' int2str(projchan) ',[' int2str(channel) '])*' tmpsig ];
	end;
end;

% outputs
% -------
outstr = '';
if ~popup
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the datas and generate output command
% --------------------------------------------
if length( options ) < 2
    options = '';
end;

% varargout{1} = sprintf('figure; pop_erpimage(%s,%d,%d,''%s'',%d,%d,{%s},[%s],''%s'',''%s''%s);', inputname(1), typeplot, channel, titleplot, smooth, decimate, typetxt, int2str(sortingwin), sortingeventfield, renorm, options);
popcom = sprintf('figure; pop_erpimage(%s,%d, [%s],[%s],''%s'',%d,%d,{%s},[%s],''%s'' %s);', inputname(1), typeplot, int2str(channel), int2str(projchan), titleplot, smooth, decimate, typetxt, int2str(sortingwin), sortingeventfield, options);

com = sprintf('%s erpimage( %s, %s, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), ''%s'', %d, %d %s);', outstr, tmpsig, events, titleplot, smooth, decimate, options);
disp('Command executed by pop_erpimage:');
disp(' '); disp(com); disp(' ');
eval(com)

if popup
	varargout{1} = popcom; % [10 '% Call: ' com];
end;

return;

% get contextual help
% -------------------
function txt = context(var, allvars, alltext);
	loc = strmatch( var, allvars);
	if ~isempty(loc)
		txt= alltext{loc(1)};
	else
		disp([ 'warning: variable ''' var ''' not found']);
		txt = '';
	end;

