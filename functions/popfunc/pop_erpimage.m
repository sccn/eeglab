% pop_erpimage() - plot an erpimage of a given EEG channel or independent
%                  component. Uses a pop-up window if only two (or three 
%                  in a specific condition) input arguments are given.
%
% Usage:
%   >> pop_erpimage(EEG, typeplot);          % pop_up window
%   >> pop_erpimage(EEG, typeplot, lastcom); % pop_up window
%   >> pop_erpimage(EEG, typeplot, channel); % do not pop-up
%   >> pop_erpimage(EEG, typeplot, channel, projchan, title, ...
%                  smooth, decimate, sortingtype, sortingwin, ...
%                  sortingeventfield, renorm, options...);
%
% Graphical interface:
%   "Channel or Component" - [edit box] Enter channel number or component
%                 number to plot. Sets the 'channel' parameter in erpimage()).
%   "Project to channel #" - [edit box] (only present when plotting independent
%                 components). Allow reprojecting the component activity 
%                 to a given channel or group of channels. Sets the
%                 the 'projchan' parameter in erpimage().
%   "Smoothing" - [text box] Smoothing parameter in number of trials.
%                 Sets the 'smooth' parameter in erpimage(). 
%   "Downsampling" - [edit box] Decimate parameter. Sets the 
%                 'decimate' parameter in erpimage()).
%   "Time limits" - [edit box] Enter the time limits in milliseconds. 
%                 Sets the two first parameters of the optional 'limit' 
%                 array in erpimage()).
%   "Figure title" - [edit box] Enter the figure title here.  If empty, a 
%                 title is automatically generated. Sets the 'title' option 
%                 in erpimage().
%   "Plot scalp map" - [checkbox] Setting this option plot a scalp map of the 
%                 channel location (or component topography) next to the 
%                 erpimage. Sets the optional 'topo' input in erpimage().
%   "plot ERP" - [checkbox] Setting this option plot the channel or component 
%                 ERP below the ERP image. Sets the 'erp' option in erpimage().
%   "Plot colorbar" - [checkbox] Plot the colorbar on the right of the 
%                 erpimage. Sets the 'cbar' option in erpimage().
%   "ERP limits" - [edit box] Set the minimum and maximum value for the ERP 
%                 plot. Sets the 3rd and 4th parameters of the optional 
%                 'limit' array in erpimage().
%   "Color limits" - [edit box] Set the color limits for the ERP image. Sets
%                 the 'caxis' option in erpimage().
%   "Epoch sorting field" - [button and edit box] Specify the event field whose 
%                 values will be used to sort the trials. Sets the
%                 'sortingeventfield' parameter in erpimage(). ???
%   "Event type(s)" - [button and edit box] Specify which event to consider. 
%                 When several events are present in individual trials, 
%                 selects which event should be processed. After selecting 
%                 events with specified event types, if several events are 
%                 still present within some trials, the event field value for 
%                 the first event in the trial is used for sorting, and a 
%                 warning is issued. Sets the 'sortingtype' parameter in
%                 erpimage(). ???
%   "Event time range" - [edit box] Allow the user to specify a time range in ms 
%                 for selecting events. This further restrains the selection of 
%                 events (see the option above). Sets the 'sortingwin' option 
%                 in erpimage(). ???
%   "rescale" - [edit box] Can be 'yes', 'no', or a Matlab formula. Sets the
%                 'renorm' option in erpimage().
%   "align" - [edit box] Set this to 'Inf' to re-align the individual trials 
%                 on the selected event values. Sets the optional 'align' 
%                 parameter in erpimage().
%   "Don't sort by value" - [checkbox] Check this box if you do not want to 
%                 sort the trials but do want to plot the selected event 
%                 values. Sets the 'nosort' option in erpimage().
%   "Don't plot value" - [checkbox] Check this box if you do not want to 
%                 plot the selected event values, but still want to sort 
%                 the data trials according to these values. Sets the
%                 'noplot' option in erpimage().
%   "Sort by phase > Frequency" - [edit box] Specify the frequency or frequency 
%                 range for sorting trials by phase. Sets the 3rd and 
%                  4th inputs to 'phasesort' in erpimage().
%   "Window center (ms)" - [edit box] Sets the 1st optional 'phasesort'  
%                 input in erpimage().
%   "Percent low-amp. trials to ignore" - [edit box] Sets the 2nd 
%                 optional 'phasesort' input in erpimage().
%   "Wavelet cycle" - [text] Number of wavelet cycles used for spectral 
%                 decomposition at the specified latency. Cannot be edited.
%   "Inter-trial coherence options > Frequency" - [edit box] Frequency for 
%                 plotting coherence. Constrained to be the same as the 
%                 "Sort by phase > Frequency" edit box. Sets the 'coher' 
%                 option in erpimage().
%   "Signif. level" - [edit box] Significance level for coherence. Sets the
%                 the optional 'signif' parameter in erpimage().
%   "Amplitude limit" - [edit box] Amplitude limits [min max] for the data power 
%                 plot at the selected frequency (plotted automatically). Sets
%                 the 5th and 6th parameter of the optional 'limit' 
%                 parameter in erpimage(). 
%   "Coher limits" - [edit box] Upper limit (<=1) for the coherence 
%                 plot. Sets the 7th and 8th parameter of the optional 'limit' 
%                 parameter in erpimage().
%   "Image amps" - [checkbox] Check this box for plotting spectral amplitude 
%                 at the selected frequency instead of raw EEG potential. Sets
%                 the 'plotamp' parameter of erpimage(). ???
%   "Plot spectrum" - [edit box] Plot the channel or component data spectrum in 
%                 the right top corner of the ERP image plot. Sets the 'spec' 
%                 option in erpimage().
%   "Baseline ampl." - [edit box] Baseline amplitude for data power plot at the 
%                 selected frequency. Sets the 7th parameter of the optional 
%                 'limit' parameter of the erpimage() function.
%   "Mark times" - [edit box] Enter value(s) to plot vertical lines at specified 
%                 latencies. Sets the optional 'vert' parameter in erpimage().
%   "More options" - [edit box] Enter 'key', 'value' sequences. Other erpimage()
%                 option not handled by this interface: 'erpstd' for 
%                 plotting the ERP standard deviation; 'auxvar' for 
%                 plotting auxilary variables; 'ampsort' for sorting trials 
%                 based on their amplitude at the selected frequency. See 
%                 help erpimage() for more information.
%
% Inputs:
%   EEG        - dataset structure
%   typeplot   - 1=channel, 0=component {default: ?}
%   lastcom    - string containing previous pop_erpimage command 
%                (from LASTCOM) or from the previous function call output.
%                The values in this function call are used as default in the
%                graphical interface.
%
% Commandline options:
%   channel    - Index of channel or component(s) to plot {default: 1}
%   projchan   - Channel to back-project the selected component(s) to. 
%                If plotting channel activity, this argument is ignored. 
%                If [], the ICA component activation is plotted {default []}.
%   title      - ['string'] Plot title {default: []}
%   smooth     - Smoothing parameter (number of trials). Default is 5. 
%                Sets the 'avewidth' parameter in erpimage().
%   decimate   - Decimate parameter (i.e. ratio of trials_in/trials_out).
%                Sets the 'decimate' parameter in erpimage(). {Default: 0}
%   sortingtype - Sorting event type(s) ([int vector]; []=all). See notes.
%                Either a string or an integer.
%   sortingwin - Sorting event window [start, end] in seconds ([]=whole epoch)
%   sortingeventfield - Sorting field name. {default: none}
%   options    - erpimage() options. Default is none. Separate the options
%                by commas. Example 'erp', 'cbar'. See erpimage() help 
%                for further details. 
%
% Outputs from pop-up: 
%   String containing the command used to evaluate this plotting function
%   (saved by eeglab() as LASTCOM) put it into 'lastcom' input to restore
%   last input parameters as defaults in a new pop_erpimage() pop-up window
%
% Outputs from commandline:
%   Same as erpimage(). Note: No outputs are returned when a window pops-up 
%   to ask for additional arguments
%   
% Notes:
%   1) A new figure is created only when the pop_up window is called, 
%      so you may call this command to draw in e.g. sbplot() axes.
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
%      eeg_getepochevent() whose help message contains more details.
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
				   '       tmpfieldnames = num2str(unique(cell2mat({EEG.event.type}))'');' ...
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
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'renorm',[], 'no'), 'tag', 'renorm' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'align'), 'tag', 'align' } ...
			   { 'Style', 'checkbox', 'string', 'Don''t plot values', 'tooltipstring', context('noplot',vars,txt), 'value', getkeyval(lastcom, 'noplot', 'present', 0), 'tag', 'noplot' } ...
			   {} ...
			   { 'Style', 'text', 'string', 'Sort trials by phase', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Frequency (Hz | minHz maxHz)', 'tooltipstring', ['sort by phase at maximum-power frequency' 10 ...
               'in the data within the range [minHz,maxHz]' 10 '(overrides frequency specified in ''coher'' flag)']  } ...
			   { 'Style', 'text', 'string', 'Percent low-amp. trials to ignore', 'tooltipstring', ['percent of trials to reject for low' ...
			    'amplitude. Else,' 10 'if prct is in the range [-100,0] -> percent to reject for high amplitude'] } ...
			   { 'Style', 'text', 'string', 'Window center (ms)', 'tooltipstring', 'Ending time fo the 3 cycle window'  } ...
			   { 'Style', 'text', 'string', 'Wavelet cycles', 'tooltipstring', '3 cycles wavelet window'  } {}...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'phasesort', [3:4]), 'tag', 'phase', 'callback', commandphase } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'phasesort', [2]), 'tag', 'phase2' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom, 'phasesort', [1]), 'tag', 'phase3' } ...
			   { 'Style', 'text', 'string', '         3'  } {}...
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
    [oldres a b res] = inputgui( geometry, uilist, 'pophelp(''erpimage'');', ...
							 fastif( typeplot, 'Channel ERP image -- pop_erpimage()', 'Component ERP image -- pop_erpimage()'));
	if isempty(oldres), return; end;

	% first rows
	% ---------
	limits(1:8)  = NaN;
	channel   	 = eval( [ '[' res.chan ']' ]);
	titleplot    = res.title;
	if isfield(res, 'projchan'), projchan = str2num(res.projchan); else, projchan = []; end;
    options = '';
	if ~typeplot
		options = [options ',''yerplabel'',''''' ];
	end;
	if isempty(titleplot)
        if typeplot==1
            if ~isempty(EEG.chanlocs) % if channel plot
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
		if ~isempty(EEG.chanlocs)
			if typeplot == 0
				     options = [options ',''topo'', { EEG.icawinv(:,' int2str(channel) ') EEG.chanlocs } '];
			else     options = [options ',''topo'', { ' int2str(channel) ' EEG.chanlocs } '];
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
		options = [ options ',''limits'',[' num2str(limits) ']' ];
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
		options = [ options ',''phasesort'',[' num2str(tmpphase) ']' ];
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
		options = [ options ',''coher'',[' num2str(tmpcoher) ']' ];
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

