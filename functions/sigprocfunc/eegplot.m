% eegplot() - Scroll (horizontally and/or vertically) through multichannel data.
%             The current version (3) allows vertical scrolling through channels 
%             and manual marking/unmarking of data stretches or epochs for rejection.
% Usage: 
%           >> eegplot(data, 'key1', value1 ...); % use interface buttons, etc.
%      else
%           >> eegplot('noui', data, 'key1', value1 ...); % no user interface
%
% Menu items:
%    "Figure > print" - [menu] Print figure in portrait or landscape.
%    "Figure > Edit figure" - [menu] Remove menus and buttons and call up the standard
%                  Matlab figure menu. Select "Tools > Edit" to format the figure
%                  for publication. Command line equivalent: 'noui' 
%    "Figure > Accept and Close" - [menu] Same as the bottom-right "Reject" button. 
%    "Figure > Cancel and Close" - [menu] Cancel all editing, same as the "Cancel" button. 
%    "Display > Marking color" > [Hide|Show] marks" - [menu] Show or hide patches of 
%                  background color behind the data. Mark stretches of *continuous* 
%                  data (e.g., for rejection) by dragging the mouse horizontally 
%                  over the activity. With *epoched* data, click on the selected epochs.
%                  Clicked on a marked region to unmark it. Called from the
%                  command line, marked data stretches or epochs are returned in 
%                  the TMPREJ variable in the global workspace *if/when* the "Reject" 
%                  button is pressed (see Outputs); called from pop_eegplot() or 
%                  eeglab(), the marked data portions are removed from the current
%                  dataset, and the dataset is automatically updated.
%     "Display > Marking color > Choose color" - [menu] Change the background marking 
%                  color. The marking color(s) of previously marked trials are preserved. 
%                  Called from command line, subsequent functions eegplot2event() or 
%                  eegplot2trials() allow processing trials marked with different colors 
%                  in the TMPREJ output variable. Command line equivalent: 'wincolor'.
%     "Display > Grid > ..." - [menu] Toggle (on or off) time and/or channel axis grids 
%                  in the activity plot. Submenus allow modifications to grid aspects.
%                  Command line equivalents: 'xgrid' / 'ygrid' 
%     "Display > Show scale" - [menu] Show (or hide if shown) the scale on the bottom 
%                  right corner of the activity window. Command line equivalent: 'scale' 
%     "Display > Title" - [menu] Change the title of the figure. Command line equivalent:
%                  'title'
%     "Settings > Time range to display"  - [menu] For continuous EEG data, this item 
%                  pops up a query window for entering the number of seconds to display
%                  in the activity window. For epoched data, the query window asks
%                  for the number of epochs to display (this can be fractional). 
%                  Command line equivalent: 'winlength'
%     "Settings > Number of channels to display" - [menu] Number of channels to display
%                  in the activity window.  If not all channels are displayed, the 
%                  user may scroll through channels using the slider on the left 
%                  of the activity plot. Command line equivalent: 'dispchans'
%     "Settings > Channel labels > ..."  - [menu] Use numbers as channel labels or load
%                  a channel location file from disk. If called from the eeglab() menu or
%                  pop_eegplot(), the channel labels of the dataset will be used. 
%                  Command line equivalent: 'eloc_file'
%     "Settings > Zoom on/off" - [menu] Toggle Matlab figure zoom on or off for time and
%                  electrode axes. left-click to zoom (x2); right-click to reverse-zoom. 
%                  Else, draw a rectange in the activity window to zoom the display into 
%                  that region. NOTE: When zoom is on, data cannot be marked for rejection.
%     "Settings > Events" - [menu] Toggle event on or off (assuming events have been 
%                  given as input). Press "legend" to pop up a legend window for events.
% Display window interface:
%    "Activity plot" - [main window] This axis displays the channel activities.  For 
%                  continuous data, the time axis shows time in seconds. For epoched
%                  data, the axis label indicate time within each epoch.
%    "Cancel" - [button] Closes the window and cancels any data rejection marks.
%    "Event types" - [button] pop up a legend window for events.
%    "<<" - [button] Scroll backwards though time or epochs by one window length.
%    "<"  - [button] Scroll backwards though time or epochs by 0.2 window length.
%    "Navigation edit box" - [edit box] Enter a starting time or epoch to jump to.
%    ">"  - [button] Scroll forward though time or epochs by 0.2 window length.
%    ">>" - [button] Scroll forward though time or epochs by one window length.
%    "Chan/Time/Value" - [text] If the mouse is within the activity window, indicates
%                  which channel, time, and activity value the cursor is closest to.
%    "Scale edit box" - [edit box] Scales the displayed amplitude in activity units.
%                  Command line equivalent: 'spacing' 
%    "+ / -" - [buttons] Use these buttons to +/- the amplitude scale by 10%. 
%    "Reject" - [button] When pressed, save rejection marks and close the figure. 
%                  Optional input parameter 'command' is evaluated at that time. 
%                  NOTE: This button's label can be redefined from the command line
%                  (see 'butlabel' below). If no processing command is specified
%                  for the 'command' parameter (below), this button does not appear.
% Required input:
%    data        - Input data matrix, either continuous 2-D (channels,timepoints) or 
%                  epoched 3-D (channels,timepoints,epochs). If the data is preceded 
%                  by keyword 'noui', GUI control elements are omitted (useful for 
%                  plotting data for presentation). A set of power spectra at
%                  each channel may also be plotted (see 'freqlimits' below).
% Optional keywords:
%    'srate'      - Sampling rate in Hz {default|0: 256 Hz}
%    'spacing'    - Display range per channel (default|0: max(data)-min(data))
%    'eloc_file'  - Electrode filename (as in  >> topoplot example) to read
%                    ascii channel labels. Else,
%                   [vector of integers] -> Show specified channel numbers. Else,
%                   [] -> Do not show channel labels {default|0 -> Show [1:nchans]}
%    'limits'     - [start end] Time limits for data epochs in ms (for labelling 
%                   purpose only).
%    'freqlimits' - [start end] If plotting epoch spectra instead of data, frequency 
%                   limits of the display. (Data should contain spectral values).
%    'winlength'  - [value] Seconds (or epochs) of data to display in window {default: 5}
%    'dispchans'  - [integer] Number of channels to display in the activity window 
%                   {default: from data}.  If < total number of channels, a vertical  
%                   slider on the left side of the figure allows vertical data scrolling. 
%    'title'      - Figure title {default: none}
%    'xgrid'      - ['on'|'off'] Toggle display of the x-axis grid {default: 'off'}
%    'ygrid'      - ['on'|'off'] Toggle display of the y-axis grid {default: 'off'}
%
% Additional keywords:
%    'command'    - ['string'] Matlab command to evaluate when the 'REJECT' button is 
%                   clicked (see Outputs below). The 'REJECT' button is visible only 
%                   if this parameter is not empty.
%    'butlabel'   - Reject button label. {default: 'REJECT'}
%    'winrej'     - [start end R G B e1 e2 e3 ...] Matrix giving data periods to mark 
%                    for rejection, each row indicating a different period: 
%                      [start end] = period limits (in frames from beginning of data); 
%                      [R G B] = specifies the marking color; 
%                      [e1 e2 e3 ...] = a (1,nchans) logical [0|1] vector giving 
%                         channels (1) to mark and (0) to not mark for rejection.
%    'color'      - ['on'|'off'|cell array] Plot channels with different colors.
%                   Entering a RGB cell array {'r' 'b' 'g'}, channels will 
%                   be plotted using the cell-array color elements. {default: 'off'}. 
%    'wincolor'   - [color] Color mark data stretches or epochs  
%                   {default: [ 0.8345 1 0.956]}
%    'events'     - [struct] EEGLAB event structure.
%    'submean'    - ['on'|'off'] Remove channel means in each window {default: 'on'}
%    'position'   - [lowleft_x lowleft_y width height] Position of the figure in pixels.
%    'tag'        - [string] Matlab object tag to identify this eegplot() window (allows 
%                    keeping track of several simultaneous eegplot() windows). 
%    'children'   - [integer] Figure handle of a *dependent* eegplot() window. Scrolling
%                    horizontally in the master window will produce the same scroll in 
%                    the dependent window. Allows comparison of two concurrent datasets.
%    'scale'      - ['on'|'off'] Display the amplitude scale {default: 'on'}.
% Outputs:
%    TMPREJ       -  Matrix (same format as 'winrej' above) placed as a variable in
%                    the global workspace (only) when the REJECT button is clicked. 
%                    The command specified in the 'command' keyword argument can use 
%                    this variable. (See eegplot2trial() and eegplot2event()). 
%
% Author: Arnaud Delorme & Colin Humphries, CNL/Salk Institute, SCCN/UCSD , 1998-2001
%
% See also: eeg_multieegplot(), eegplot2event(), eegplot2trial(), eeglab()

% deprecated 
%    'colmodif'   - nested cell array of window colors that may be marked/unmarked. Default
%                   is current color only.

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme & Colin Humphries, Salk Institute, arno@salk.edu
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
% Revision 1.87  2003/07/22 17:29:54  arno
% adding event button
%
% Revision 1.86  2003/07/20 19:27:37  scott
% "electrodes" -> "channels"
%
% Revision 1.85  2003/07/20 19:25:10  scott
% fixed typo
%
% Revision 1.84  2003/07/14 23:03:57  arno
% handling numerical types
%
% Revision 1.83  2003/07/12 02:13:47  arno
% fixing color, implementing dash
%
% Revision 1.82  2003/07/12 01:23:46  arno
% now handling events
%
% Revision 1.81  2003/06/28 02:07:18  arno
% update header
%
% Revision 1.80  2003/06/28 02:04:36  arno
% debuging zoom and scale
%
% Revision 1.79  2003/05/23 16:05:47  arno
% debug last
%
% Revision 1.78  2003/05/23 15:12:11  arno
% same
%
% Revision 1.77  2003/05/23 15:11:13  arno
% allowing cell array of color as input
%
% Revision 1.76  2003/03/13 03:04:16  scott
% commenting error msg
%
% Revision 1.75  2003/03/12 03:15:48  arno
% adding help menu
%
% Revision 1.74  2003/02/20 20:56:38  scott
% header edit -sm
%
% Revision 1.73  2003/02/20 20:11:33  arno
% replaing ???
% /
%
% Revision 1.72  2003/02/20 19:42:01  scott
% edit header msg -sm
%
% Revision 1.71  2003/02/20 18:35:02  arno
% adding interface description and new menu option to edit figure
%
% Revision 1.70  2003/02/17 03:04:57  arno
% typo in header
%
% Revision 1.69  2003/01/10 01:17:20  arno
% default position once more
%
% Revision 1.68  2003/01/10 01:11:15  arno
% new default position
%
% Revision 1.67  2003/01/10 01:04:31  arno
% same, making the figure visible at the end of plot
%
% Revision 1.66  2003/01/10 01:01:29  arno
% updating plotting position
%
% Revision 1.65  2002/11/15 01:11:59  arno
% debugging incallback
%
% Revision 1.64  2002/11/14 17:03:41  arno
% debugging multiple window selection
%
% Revision 1.63  2002/11/13 00:47:45  arno
% debug multiple color display
%
% Revision 1.62  2002/11/12 23:12:49  arno
% compatibility if one extra channel
%
% Revision 1.61  2002/10/22 17:25:56  arno
% add max for selecting regions
%
% Revision 1.60  2002/10/22 17:12:38  arno
% debug 0 limit
%
% Revision 1.59  2002/10/20 00:13:46  arno
% still debuging very low freqs
%
% Revision 1.58  2002/10/19 23:18:57  arno
% debug lim at very low frequencies
%
% Revision 1.57  2002/10/19 23:10:16  arno
% debug last
%
% Revision 1.56  2002/10/19 23:07:24  arno
% implement more exact limits at very low freq.
%
% Revision 1.55  2002/10/17 18:45:09  arno
% implementing nan
%
% Revision 1.54  2002/10/17 00:42:13  arno
% debug last
%
% Revision 1.53  2002/10/16 23:06:34  arno
% default spacing 1
%
% Revision 1.52  2002/10/04 17:57:12  arno
% Y grid text
%
% Revision 1.51  2002/09/06 19:49:03  arno
% custum colors for channels
%
% Revision 1.50  2002/09/06 19:34:31  arno
% noui optimize
%
% Revision 1.49  2002/09/06 19:29:40  arno
% addscale scale argument to toggle on/off the scale at startup
%
% Revision 1.48  2002/09/05 15:03:43  arno
% debug scale
%
% Revision 1.47  2002/08/28 00:42:57  arno
% debugging if no arguments
%
% Revision 1.46  2002/08/27 23:34:00  arno
% debugging xgrid
%
% Revision 1.45  2002/08/19 19:47:19  arno
% debugging last
%
% Revision 1.44  2002/08/19 19:33:57  arno
% bebugging for mac
%
% Revision 1.43  2002/08/15 01:05:35  arno
% debug continuous
%
% Revision 1.42  2002/08/15 00:45:47  arno
% typo
%
% Revision 1.41  2002/08/14 01:20:41  arno
% close or cancel
%
% Revision 1.40  2002/08/13 00:06:29  arno
% [Aslidder color and scale title
%
% Revision 1.39  2002/08/12 18:58:58  arno
% inputdlg2
%
% Revision 1.38  2002/08/12 01:18:40  arno
% debug last
%
% Revision 1.37  2002/08/12 01:15:09  arno
% update color
%
% Revision 1.36  2002/08/12 00:27:16  arno
% color
%
% Revision 1.35  2002/08/11 23:10:04  arno
% text edit
%
% Revision 1.34  2002/08/11 23:07:48  arno
% done for title
%
% Revision 1.33  2002/08/11 22:47:11  arno
% update color
%
% Revision 1.32  2002/08/08 21:03:56  arno
% debugging continuous window select
%
% Revision 1.31  2002/08/08 00:24:43  arno
% header
%
% Revision 1.30  2002/08/08 00:20:59  arno
% adding colmodif option
%
% Revision 1.29  2002/08/07 17:38:26  arno
% debugging
%
% Revision 1.28  2002/07/31 16:53:36  arno
% changing button size
%
% Revision 1.27  2002/07/31 16:08:35  arno
% debugging epoch problem
%
% Revision 1.26  2002/07/31 01:08:46  arno
% resizing button
%
% Revision 1.25  2002/07/31 01:07:06  arno
% debugging but label
%
% Revision 1.24  2002/07/31 01:05:07  arno
% optional button label
%
% Revision 1.23  2002/07/30 19:37:49  arno
% debug wincolor
%
% Revision 1.22  2002/07/30 19:22:57  arno
% creating new argument for selecting colors
%
% Revision 1.21  2002/07/30 17:15:06  arno
% adding multiple window plotting
%
% Revision 1.20  2002/07/30 15:16:41  arno
% debugging frequency axis
%
% Revision 1.19  2002/07/26 21:58:03  arno
% same
%
% Revision 1.18  2002/07/26 21:57:21  arno
% debugging x axis
%
% Revision 1.17  2002/07/24 01:37:20  arno
% debugging submean
%
% Revision 1.16  2002/07/24 01:28:09  arno
% default spacing -> 3 stds
%
% Revision 1.15  2002/07/24 01:25:19  arno
% default amplitude to standard deviation
%
% Revision 1.14  2002/07/23 16:26:01  arno
% debugging slidder and X axis
%
% Revision 1.13  2002/07/22 23:14:40  arno
% adding one to the epoch number
%
% Revision 1.12  2002/06/27 01:48:02  arno
% removing frequency option
%
% Revision 1.11  2002/06/27 01:42:48  scott
% completed help message editing -sm & ad
%
% Revision 1.10  2002/06/26 22:12:30  arno
% debugging slider and editing
%
% Revision 1.9  2002/06/26 21:29:14  arno
% editing header
%
% Revision 1.8  2002/06/26 18:39:33  scott
% Edited help message -- see ??? for points needing more clarification. -sm
%
% Revision 1.7  2002/06/25 01:19:49  arno
% update tag position
%
% Revision 1.6  2002/06/25 01:05:45  arno
% new version with zoom
%
% Revision 1.3  2002/05/02 21:27:01  arno
% reject button debugging
%
% Revision 1.2  2002/04/30 18:16:33  arno
% conidtional REJECT button apparition
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 4/01 version 3 Arnaud Delorme, CNL / Salk Institute, La Jolla CA (arno@salk.edu)
% 4-28-01 change popup windows; regroup drawing functions; add the trials tags  -ad
% 4-29-01 colored timerange selection with mouse; g.time, electrode and value online display -ad 
% 4-30-01 reorganize menus and buttons, add new ones; electrode selection with mouse; output -ad 
% 5-01-01 multisignal display, add boudary contraints  -ad
% 9-11-01 debugging, inversing Y direction for everything -ad
% 9-16-01 normalization of the positions of button (PC and UNIX compatibility) -ad
% 9-26-01 using 'key', value, matlab calling convention -ad
% 10-10-01 adding 'xgrid', 'ygrid', 'color', 'freq' options  -ad
% 11-10-01 whole restructuration of userdata variable, regroupment in g -ad
% 12-10-01 add trial number (top) and trials limits (bottom) for EEG trials -ad   
% 01-25-02 reformated help & license -ad 
% 03-04-02 preserve data 3D structure (to avoid reserving new memory:worked) -ad 
% 03-08-02 debug cancel and reject button -ad 
% 03-15-02 added readlocs and the use of eloc input structure -ad 
% 03-15-02 re-program noui version -ad 
% 03-17-02 debugging and text -ad & sm
% 03-20-02 more complex names for the callback variables -ad & sm
% 03-22-02 change help message -ad & lf
% ---------------------------------------------------------------------- 
% from an original version by Colin Humphries, 5/98  
% CNL / Salk Institute, La Jolla CA (colin@salk.edu) 
% 5-14-98 v2.1 fixed bug for small-variance data -ch
% 1-31-00 v2.2 exchanged meaning of > and >>, < and << -sm
% 8-15-00 v2.3 turned on SPACING_EYE and added int vector input for eloc_file -sm
% 12-16-00 added undocumented figure position arg (if not 'noui') -sm

% internal variables structure
% All in g except for Eposition and Eg.spacingwhich are inside the boxes
%
% gcf
%    1 - winlength
%    2 - srate 
%    3 - children
% 'backeeg' axis
%    1 - trialtag
%    2 - g.winrej
%    3 - nested call flag
% 'eegaxis'
%    1 - data
%    2 - colorlist
%    3 - submean    % on or off, subtract the mean
%    4 - maxfreq    % empty [] if no gfrequency content
% 'buttons hold other informations' Eposition for instance hold the current postition

function [outvar1] = eegplot(data, varargin); % p1,p2,p3,p4,p5,p6,p7,p8,p9)

% Defaults (can be re-defined):

DEFAULT_PLOT_COLOR = { [0 0 1], [0.7 0.7 0.7]};         % EEG line color
try, icadefs;
	DEFAULT_FIG_COLOR = BACKCOLOR;
	BUTTON_COLOR = GUIBUTTONCOLOR;
catch
	DEFAULT_FIG_COLOR = [1 1 1];
	BUTTON_COLOR =[0.8 0.8 0.8];
end;
DEFAULT_AXIS_COLOR = 'k';         % X-axis, Y-axis Color, text Color
DEFAULT_GRID_SPACING = 1;         % Grid lines every n seconds
DEFAULT_GRID_STYLE = '-';         % Grid line style
YAXIS_NEG = 'off';                % 'off' = positive up 
DEFAULT_NOUI_PLOT_COLOR = 'k';    % EEG line color for noui option
                                  %   0 - 1st color in AxesColorOrder
SPACING_EYE = 'on';               % g.spacingI on/off
SPACING_UNITS_STRING = '';        % '\muV' for microvolt optional units for g.spacingI Ex. uV
DEFAULT_AXES_POSITION = [0.0964286 0.15 0.842 0.788095];
                                  % dimensions of main EEG axes
ORIGINAL_POSITION = [50 50 800 500];
                                  
if nargin < 1
   help eegplot
   return
end
				  
% %%%%%%%%%%%%%%%%%%%%%%%%
% Setup inputs
% %%%%%%%%%%%%%%%%%%%%%%%%

if ~isstr(data) % If NOT a 'noui' call or a callback from uicontrols

   try
       options = varargin;
       for index = 1:length(options)
           if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
       end;
       if ~isempty( varargin ), g=struct(options{:}); 
       else g= []; end;
   catch
   		disp('eegplot() error: calling convention {''key'', value, ... } error'); return;
   end;	
		 
   try, g.srate; 		    catch, g.srate		= 256; 	end;
   try, g.spacing; 			catch, g.spacing	= 0; 	end;
   try, g.eloc_file; 		catch, g.eloc_file	= 0; 	end; % 0 mean numbered
   try, g.winlength; 		catch, g.winlength	= 5; 	end; % Number of seconds of EEG displayed
   try, g.position; 	    catch, g.position	= ORIGINAL_POSITION; 	end;
   try, g.title; 		    catch, g.title		= ['Scroll activity -- eegplot()']; 	end;
   try, g.trialstag; 		catch, g.trialstag	= -1; 	end;
   try, g.winrej; 			catch, g.winrej		= []; 	end;
   try, g.command; 			catch, g.command	= ''; 	end;
   try, g.tag; 				catch, g.tag		= 'EEGPLOT'; end;
   try, g.xgrid;		    catch, g.xgrid		= 'off'; end;
   try, g.ygrid;		    catch, g.ygrid		= 'off'; end;
   try, g.color;		    catch, g.color		= 'off'; end;
   try, g.submean;			catch, g.submean	= 'on'; end;
   try, g.children;			catch, g.children	= 0; end;
   try, g.limits;		    catch, g.limits	    = [0 1000*(size(data,2)-1)/g.srate]; end;
   try, g.freqlimits;	    catch, g.freqlimits	= []; end;
   try, g.dispchans; 		catch, g.dispchans  = size(data,1); end;
   try, g.wincolor; 		catch, g.wincolor   = [ 0.8345 1 0.9560]; end;
   try, g.butlabel; 		catch, g.butlabel   = 'REJECT'; end;
   try, g.colmodif; 		catch, g.colmodif   = { g.wincolor }; end;
   try, g.scale; 		    catch, g.scale      = 'on'; end;
   try, g.events; 		    catch, g.events      = []; end;

   if ndims(data) > 2
   		g.trialstag = size(	data, 2);
   	end;	
      
   gfields = fieldnames(g);
   for index=1:length(gfields)
      switch gfields{index}
      case {'spacing', 'srate' 'eloc_file' 'winlength' 'position' 'title' ...
               'trialstag'  'winrej' 'command' 'tag' 'xgrid' 'ygrid' 'color' 'colmodif'...
               'freqlimits' 'submean' 'children' 'limits' 'dispchans' 'wincolor' 'butlabel' 'scale' 'events' },;
      otherwise, error(['eegplot: unrecognized option: ''' gfields{index} '''' ]);
      end;
   end;

   if length(g.srate) > 1
   		disp('Error: srate must be a single number'); return;
   end;	
   if length(g.spacing) > 1
   		disp('Error: ''spacing'' must be a single number'); return;
   end;	
   if length(g.winlength) > 1
   		disp('Error: winlength must be a single number'); return;
   end;	
   if isstr(g.title) > 1
   		disp('Error: title must be is a string'); return;
   end;	
   if isstr(g.command) > 1
   		disp('Error: command must be is a string'); return;
   end;	
   if isstr(g.tag) > 1
   		disp('Error: tag must be is a string'); return;
   end;	
   if length(g.position) ~= 4
   		disp('Error: position must be is a 4 elements array'); return;
   end;	
   switch lower(g.xgrid)
	   case { 'on', 'off' },; 
	   otherwise disp('Error: xgrid must be either ''on'' or ''off'''); return;
   end;	
   switch lower(g.ygrid)
	   case { 'on', 'off' },; 
	   otherwise disp('Error: ygrid must be either ''on'' or ''off'''); return;
   end;	
   switch lower(g.submean)
	   case { 'on' 'off' };
	   otherwise disp('Error: submean must be either ''on'' or ''off'''); return;
   end;	
   switch lower(g.scale)
	   case { 'on' 'off' };
	   otherwise disp('Error: scale must be either ''on'' or ''off'''); return;
   end;	
   
   if ~iscell(g.color)
	   switch lower(g.color)
		case 'on', g.color = { 'k', 'm', 'c', 'b', 'g' }; 
		case 'off', g.color = { [ 0 0 0.4] };  
		otherwise 
		 disp('Error: color must be either ''on'' or ''off'' or a cell array'); return;
	   end;	
   end;
   if length(g.dispchans) > size(data,1)
	   g.dispchans = size(data,1);
   end;
   if ~iscell(g.colmodif)
   		g.colmodif = { g.colmodif };
   end;
   if any(isnan(data(:))) & strcmpi(g.submean, 'on')
       g.submean = 'nan';
   end;
   
   % convert color to modify into array of float
   % -------------------------------------------
   for index = 1:length(g.colmodif)
	   tmpcolmodif(index) = g.colmodif{index}(1) + g.colmodif{index}(2)*10 + g.colmodif{index}(3)*100;
   end;
   g.colmodif = tmpcolmodif;
   
   [g.chans,g.frames, tmpnb] = size(data);
   g.frames = g.frames*tmpnb;
  
  if g.spacing == 0
    maxindex = min(1000, g.frames);  
	stds = std(data(:,1:maxindex),[],2);
	stds = sort(stds);
	if length(stds) > 2
		stds = mean(stds(2:end-1));
	else
		stds = mean(stds);
	end;	
    g.spacing = stds*3;  
    if g.spacing > 10
      g.spacing = round(g.spacing);
    end
    if g.spacing  == 0 | isnan(g.spacing)
        g.spacing = 1;
    end;
  end

  % set defaults
  % ------------ 
  g.incallback = 0;
  g.winstatus = 1;
  g.setelectrode  = 0;
  [g.chans,g.frames,tmpnb] = size(data);   
  g.frames = g.frames*tmpnb;
  g.nbdat = 1; % deprecated
  g.time  = 0;
  g.elecoffset = 0;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%
  % Prepare figure and axes
  % %%%%%%%%%%%%%%%%%%%%%%%%
  
  figh = figure('UserData', g,... % store the settings here
      'Color',DEFAULT_FIG_COLOR, 'name', g.title,...
      'MenuBar','none','tag', g.tag ,'Position',g.position, 'numbertitle', 'off', 'visible', 'off');

  pos = get(figh,'position'); % plot relative to current axes
  q = [pos(1) pos(2) 0 0];
  s = [pos(3) pos(4) pos(3) pos(4)]./100;
  clf;
      
  % Background axis
  % --------------- 
  ax0 = axes('tag','backeeg','parent',figh,...
      'Position',DEFAULT_AXES_POSITION,...
      'Box','off','xgrid','off', 'xaxislocation', 'top'); 

  % Drawing axis
  % --------------- 
  YLabels = num2str((1:g.chans)');  % Use numbers as default
  YLabels = flipud(str2mat(YLabels,' '));
  ax1 = axes('Position',DEFAULT_AXES_POSITION,...
             'userdata', data, ...% store the data here (when in g, slow down display)
   			'tag','eegaxis','parent',figh,...
      'Box','on','xgrid', g.xgrid,'ygrid', g.ygrid,...
      'gridlinestyle',DEFAULT_GRID_STYLE,...
      'Xlim',[0 g.winlength*g.srate],...
      'xtick',[0:g.srate*DEFAULT_GRID_SPACING:g.winlength*g.srate],...
      'Ylim',[0 (g.chans+1)*g.spacing],...
      'YTick',[0:g.spacing:g.chans*g.spacing],...
      'YTickLabel', YLabels,...
      'XTickLabel',num2str((0:DEFAULT_GRID_SPACING:g.winlength)'),...
      'TickLength',[.005 .005],...
      'Color','none',...
      'XColor',DEFAULT_AXIS_COLOR,...
      'YColor',DEFAULT_AXIS_COLOR);

  if isstr(g.eloc_file) | isstruct(g.eloc_file)  % Read in electrode names
      if isstruct(g.eloc_file) & length(g.eloc_file) > size(data,1)
          g.eloc_file(end) = []; % common reference channel location
      end;
      eegplot('setelect', g.eloc_file, ax1);
  end;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%
  % Set up uicontrols
  % %%%%%%%%%%%%%%%%%%%%%%%%%

% positions of buttons
  posbut(1,:) = [ 0.0464    0.0254    0.0385    0.0339 ]; % <<
  posbut(2,:) = [ 0.0924    0.0254    0.0288    0.0339 ]; % <
  posbut(3,:) = [ 0.1924    0.0254    0.0299    0.0339 ]; % >
  posbut(4,:) = [ 0.2297    0.0254    0.0385    0.0339 ]; % >>
  posbut(5,:) = [ 0.1287    0.0203    0.0561    0.0390 ]; % Eposition
  posbut(6,:) = [ 0.4744    0.0236    0.0582    0.0390 ]; % Espacing
  posbut(7,:) = [ 0.2762    0.01    0.0582    0.0390 ]; % elec
  posbut(8,:) = [ 0.3256    0.01    0.0707    0.0390 ]; % g.time
  posbut(9,:) = [ 0.4006    0.01    0.0582    0.0390 ]; % value
  posbut(14,:) = [ 0.2762    0.05    0.0582    0.0390 ]; % elec tag
  posbut(15,:) = [ 0.3256    0.05    0.0707    0.0390 ]; % g.time tag
  posbut(16,:) = [ 0.4006    0.05    0.0582    0.0390 ]; % value tag
  posbut(10,:) = [ 0.5437    0.0458    0.0275    0.0270 ]; % +
  posbut(11,:) = [ 0.5437    0.0134    0.0275    0.0270 ]; % -
  posbut(12,:) = [ 0.6    0.02    0.14    0.05 ]; % cancel
  posbut(13,:) = [-0.15   0.02    0.07    0.05 ]; % cancel
  posbut(17,:) = [-0.06    0.02    0.09    0.05 ]; % events types
  posbut(20,:) = [-0.17   0.15     0.015    0.8 ]; % slider
  posbut(:,1) = posbut(:,1)+0.2;

% Five move buttons: << < text > >> 

  u(1) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position', posbut(1,:), ...
	'Tag','Pushbutton1',...
	'string','<<',...
	'Callback','eegplot(''drawp'',1)');
  u(2) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position', posbut(2,:), ...
	'Tag','Pushbutton2',...
	'string','<',...
	'Callback','eegplot(''drawp'',2)');
  u(5) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',[1 1 1], ...
	'Position', posbut(5,:), ...
	'Style','edit', ...
	'Tag','EPosition',...
	'string', fastif(g.trialstag(1) == -1, '0', '1'),...
	'Callback', 'eegplot(''drawp'',0);' );
  u(3) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position',posbut(3,:), ...
	'Tag','Pushbutton3',...
	'string','>',...
	'Callback','eegplot(''drawp'',3)');
  u(4) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position',posbut(4,:), ...
	'Tag','Pushbutton4',...
	'string','>>',...
	'Callback','eegplot(''drawp'',4)');

% Text edit fields: ESpacing

  u(6) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',[1 1 1], ...
	'Position', posbut(6,:), ...
	'Style','edit', ...
	'Tag','ESpacing',...
	'string',num2str(g.spacing),...
	'Callback', 'eegplot(''draws'',0);' );

% Slider for vertical motion
  u(20) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position', posbut(20,:), ...
   'Style','slider', ...
   'visible', 'off', ...
   'sliderstep', [0.9 1], ...
   'Tag','eegslider', ...
   'callback', [ 'tmpg = get(gcbf, ''userdata'');' ... 
   				'tmpg.elecoffset = get(gcbo, ''value'')*(tmpg.chans-tmpg.dispchans);' ...
               'set(gcbf, ''userdata'', tmpg);' ...
               'eegplot(''drawp'',0);' ...
               'clear tmpg;' ], ...
   'value', 0);

% channels, position, value and tag

  u(9) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',DEFAULT_FIG_COLOR, ...
	'Position', posbut(7,:), ...
	'Style','text', ...
	'Tag','Eelec',...
	'string',' ');
  u(10) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',DEFAULT_FIG_COLOR, ...
	'Position', posbut(8,:), ...
	'Style','text', ...
	'Tag','Etime',...
	'string','0.00');
  u(11) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',DEFAULT_FIG_COLOR, ...
	'Position',posbut(9,:), ...
	'Style','text', ...
	'Tag','Evalue',...
	'string','0.00');

  u(14)= uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',DEFAULT_FIG_COLOR, ...
	'Position', posbut(14,:), ...
	'Style','text', ...
	'Tag','Eelecname',...
	'string','Chan.');
  u(15) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',DEFAULT_FIG_COLOR, ...
	'Position', posbut(15,:), ...
	'Style','text', ...
	'Tag','Etimename',...
	'string','Time');
  u(16) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'BackgroundColor',DEFAULT_FIG_COLOR, ...
	'Position',posbut(16,:), ...
	'Style','text', ...
	'Tag','Evaluename',...
	'string','Value');

% ESpacing buttons: + -
  u(7) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position',posbut(10,:), ...
	'Tag','Pushbutton5',...
	'string','+',...
	'FontSize',8,...
	'Callback','eegplot(''draws'',1)');
  u(8) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position',posbut(11,:), ...
	'Tag','Pushbutton6',...
	'string','-',...
	'FontSize',8,...
	'Callback','eegplot(''draws'',2)');

  if isempty(g.command) tmpcom = 'fprintf(''Rejections saved in variable TMPREJ\n'');';   
  else tmpcom = g.command;
  end;
  acceptcommand = [ 'g = get(gcbf, ''userdata'');' ... 
				    'TMPREJ = g.winrej;' ...
                    'if g.children, delete(g.children); end;' ...
                    'delete(gcbf);' ...
		  				  tmpcom ...
                    'clear g;']; % quitting expression
  if ~isempty(g.command)
	  u(12) = uicontrol('Parent',figh, ...
						'Units', 'normalized', ...
						'Position',posbut(12,:), ...
						'Tag','Accept',...
						'string',g.butlabel, 'callback', acceptcommand);
  end;
  u(13) = uicontrol('Parent',figh, ...
	'Units', 'normalized', ...
	'Position',posbut(13,:), ...
	'string',fastif(isempty(g.command),'CLOSE', 'CANCEL'), 'callback', ...
		[	'g = get(gcbf, ''userdata'');' ... 
            'if g.children, delete(g.children); end;' ...
			'close(gcbf);'] );

  if ~isempty(g.events)
      u(17) = uicontrol('Parent',figh, ...
                        'Units', 'normalized', ...
                        'Position',posbut(17,:), ...
                        'string', 'Event types', 'callback', 'eegplot(''drawlegend'', gcbf)');
  end;

  set(u,'Units','Normalized')
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set up uimenus
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Figure Menu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m(7) = uimenu('Parent',figh,'Label','Figure');
  m(8) = uimenu('Parent',m(7),'Label','Print');
  uimenu('Parent',m(7),'Label','Edit figure', 'Callback', 'eegplot(''noui'');');
  uimenu('Parent',m(7),'Label','Accept and close', 'Callback', acceptcommand );
  uimenu('Parent',m(7),'Label','Cancel and close', 'Callback','delete(gcbf)')
  
  % Portrait %%%%%%%%
  timestring = ['[OBJ1,FIG1] = gcbo;',...
	        'PANT1 = get(OBJ1,''parent'');',...
	        'OBJ2 = findobj(''tag'',''orient'',''parent'',PANT1);',...
		'set(OBJ2,''checked'',''off'');',...
		'set(OBJ1,''checked'',''on'');',...
		'set(FIG1,''PaperOrientation'',''portrait'');',...
		'clear OBJ1 FIG1 OBJ2 PANT1;'];
		
  uimenu('Parent',m(8),'Label','Portrait','checked',...
      'on','tag','orient','callback',timestring)
  
  % Landscape %%%%%%%
  timestring = ['[OBJ1,FIG1] = gcbo;',...
	        'PANT1 = get(OBJ1,''parent'');',...
	        'OBJ2 = findobj(''tag'',''orient'',''parent'',PANT1);',...
		'set(OBJ2,''checked'',''off'');',...
		'set(OBJ1,''checked'',''on'');',...
		'set(FIG1,''PaperOrientation'',''landscape'');',...
		'clear OBJ1 FIG1 OBJ2 PANT1;'];
  
  uimenu('Parent',m(8),'Label','Landscape','checked',...
      'off','tag','orient','callback',timestring)

  % Print command %%%%%%%
  uimenu('Parent',m(8),'Label','Print','tag','printcommand','callback',...
  		['RESULT = inputdlg2( { ''Enter command:'' }, ''Print'', 1,  { ''print -r72'' });' ...
		 'if size( RESULT,1 ) ~= 0' ... 
		 '  eval ( RESULT{1} );' ...
		 'end;' ...
		 'clear RESULT;' ]);
  
  % Display Menu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m(1) = uimenu('Parent',figh,...
      'Label','Display', 'tag', 'displaymenu');
  
  % window grid %%%%%%%%%%%%%
  % userdata = 4 cells : display yes/no, color, electrode yes/no, trial boundary adapt yes/no (1/0)  
  m(11) = uimenu('Parent',m(1),'Label','Marking color', 'tag', 'displaywin', ...
				 'userdata', { 1, [0.8 1 0.8], 0, fastif( g.trialstag(1) == -1, 0, 1)} );
  uimenu('Parent',m(11),'Label','Hide marks','Callback', ...
  	['g = get(gcbf, ''userdata'');' ...
  	 'if ~g.winstatus' ... 
  	 '  set(gcbo, ''label'', ''Hide marks'');' ...
  	 'else' ...
  	 '  set(gcbo, ''label'', ''Show marks'');' ...
  	 'end;' ...
  	 'g.winstatus = ~g.winstatus;' ...
  	 'set(gcbf, ''userdata'', g);' ...
  	 'eegplot(''drawb''); clear g;'] )

	% color
	uimenu('Parent',m(11),'Label','Choose color', 'Callback', ...
  	[ 'g = get(gcbf, ''userdata'');' ...
  	  'g.wincolor = uisetcolor(g.wincolor);' ...
      'set(gcbf, ''userdata'', g ); ' ...
      'clear g;'] )

	% set channels
	uimenu('Parent',m(11),'Label','Mark channels', 'enable', 'off', 'checked', 'off', 'Callback', ...
  	['g = get(gcbf, ''userdata'');' ...
  	 'g.setelectrode = ~g.setelectrode;' ...
  	 'set(gcbf, ''userdata'', g); ' ...
     'if ~g.setelectrode setgcbo, ''checked'', ''on''); else set(gcbo, ''checked'', ''off''); end;'...
     ' clear g;'] )

	% trials boundaries
	%uimenu('Parent',m(11),'Label','Trial boundaries', 'checked', fastif( g.trialstag(1) == -1, 'off', 'on'), 'Callback', ...
  	%['hh = findobj(''tag'',''displaywin'',''parent'', findobj(''tag'',''displaymenu'',''parent'', gcbf ));' ...
  	% 'hhdat = get(hh, ''userdata'');' ...
  	% 'set(hh, ''userdata'', { hhdat{1},  hhdat{2}, hhdat{3}, ~hhdat{4}} ); ' ...
    %'if ~hhdat{4} set(gcbo, ''checked'', ''on''); else set(gcbo, ''checked'', ''off''); end;' ... 
    %' clear hh hhdat;'] )

  % X grid %%%%%%%%%%%%
  m(3) = uimenu('Parent',m(1),'Label','Grid');
  timestring = ['FIGH = gcbf;',...
	            'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
	            'if size(get(AXESH,''xgrid''),2) == 2' ... %on
		        '  set(AXESH,''xgrid'',''off'');',...
		        '  set(gcbo,''label'',''X grid on'');',...
		        'else' ...
		        '  set(AXESH,''xgrid'',''on'');',...
		        '  set(gcbo,''label'',''X grid off'');',...
		        'end;' ...
		        'clear FIGH AXESH;' ];
  uimenu('Parent',m(3),'Label',fastif(strcmp(g.xgrid, 'off'), 'X grid on','X grid off'), 'Callback',timestring)
  
  % Y grid %%%%%%%%%%%%%
  timestring = ['FIGH = gcbf;',...
	            'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
	            'if size(get(AXESH,''ygrid''),2) == 2' ... %on
		        '  set(AXESH,''ygrid'',''off'');',...
		        '  set(gcbo,''label'',''Y grid on'');',...
		        'else' ...
		        '  set(AXESH,''ygrid'',''on'');',...
		        '  set(gcbo,''label'',''Y grid off'');',...
		        'end;' ...
		        'clear FIGH AXESH;' ];
  uimenu('Parent',m(3),'Label',fastif(strcmp(g.ygrid, 'off'), 'Y grid on','Y grid off'), 'Callback',timestring)

  % Grid Style %%%%%%%%%
  m(5) = uimenu('Parent',m(3),'Label','Grid Style');
  timestring = ['FIGH = gcbf;',...
	        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
		'set(AXESH,''gridlinestyle'',''--'');',...
		'clear FIGH AXESH;'];
  uimenu('Parent',m(5),'Label','- -','Callback',timestring)
  timestring = ['FIGH = gcbf;',...
	        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
		'set(AXESH,''gridlinestyle'',''-.'');',...
		'clear FIGH AXESH;'];
  uimenu('Parent',m(5),'Label','_ .','Callback',timestring)
  timestring = ['FIGH = gcbf;',...
	        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
		'set(AXESH,''gridlinestyle'','':'');',...
		'clear FIGH AXESH;'];
  uimenu('Parent',m(5),'Label','. .','Callback',timestring)
  timestring = ['FIGH = gcbf;',...
	        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
		'set(AXESH,''gridlinestyle'',''-'');',...
		'clear FIGH AXESH;'];
  uimenu('Parent',m(5),'Label','__','Callback',timestring)
  
  % Scale Eye %%%%%%%%%
  timestring = ['[OBJ1,FIG1] = gcbo;',...
	        'eegplot(''scaleeye'',OBJ1,FIG1);',...
		'clear OBJ1 FIG1;'];
  m(7) = uimenu('Parent',m(1),'Label','Show scale','Callback',timestring);
  
  % Title %%%%%%%%%%%%
  uimenu('Parent',m(1),'Label','Title','Callback','eegplot(''title'')')
  
  % Settings Menu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m(2) = uimenu('Parent',figh,...
      'Label','Settings'); 
  
  % Window %%%%%%%%%%%%
  uimenu('Parent',m(2),'Label','Time range to display',...
      'Callback','eegplot(''window'')')
  
  % Electrode window %%%%%%%%
  uimenu('Parent',m(2),'Label','Number of channels to display',...
      'Callback','eegplot(''winelec'')')
  
  % Electrodes %%%%%%%%
  m(6) = uimenu('Parent',m(2),'Label','Channel labels');
  
  timestring = ['FIGH = gcbf;',...
	        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
		'YTICK = get(AXESH,''YTick'');',...
		'YTICK = length(YTICK);',...
		'set(AXESH,''YTickLabel'',flipud(str2mat(num2str((1:YTICK-1)''),'' '')));',...
		'clear FIGH AXESH YTICK;'];
  uimenu('Parent',m(6),'Label','Show number','Callback',timestring)
  uimenu('Parent',m(6),'Label','Load .loc(s) file',...
      'Callback','eegplot(''loadelect'');')
  
  % Zooms %%%%%%%%
  zm = uimenu('Parent',m(2),'Label','Zoom off/on');
  commandzoom = [ 'tmpstr = get(gcbf, ''windowbuttondownfcn'');' ...
                  'set(gcbf, ''windowbuttondownfcn'', [ tmpstr ''; eegplot(''''zoom'''', gcbf, 1);'' ]);' ...
                  'tmpg = get(gcbf, ''userdata'');' ...
                  'set(gcbf, ''windowbuttonmotionfcn'', tmpg.commandselect{2}); clear tmpg tmpstr;'];
    
  %uimenu('Parent',zm,'Label','Zoom time', 'callback', ...
  %             [ 'zoom(gcbf, ''xon'');' commandzoom ]);
  %uimenu('Parent',zm,'Label','Zoom channels', 'callback', ...
  %             [ 'zoom(gcbf, ''yon'');' commandzoom ]);
  uimenu('Parent',zm,'Label','Zoom on', 'callback', ...
               [ 'zoom(gcbf, ''on'');' commandzoom ]);
  uimenu('Parent',zm,'Label','Zoom off', 'separator', 'on', 'callback', ...
     ['zoom(gcbf, ''off''); tmpg = get(gcbf, ''userdata'');' ...
     'set(gcbf, ''windowbuttondownfcn'', tmpg.commandselect{1});' ...
     'set(gcbf, ''windowbuttonmotionfcn'', tmpg.commandselect{2});' ...
     'set(gcbf, ''windowbuttonupfcn'', tmpg.commandselect{3});' ...
     'clear tmpg;' ]);
  uimenu('Parent',figh,'Label', 'Help', 'callback', 'pophelp(''eegplot'');');

  % Events %%%%%%%%
  zm = uimenu('Parent',m(2),'Label','Events');
  complotevent = [ 'tmpg = get(gcbf, ''userdata'');' ...
                  'tmpg.plotevent = ''on'';' ...                  
                  'set(gcbf, ''userdata'', tmpg); clear tmpg; eegplot(''drawp'', 0);'];
  comnoevent   = [ 'tmpg = get(gcbf, ''userdata'');' ...
                  'tmpg.plotevent = ''off'';' ...                  
                  'set(gcbf, ''userdata'', tmpg); clear tmpg; eegplot(''drawp'', 0);'];
  comeventleg  = [ 'eegplot(''drawlegend'', gcbf);'];
    
  uimenu('Parent',zm,'Label','Events on'    , 'callback', complotevent, 'enable', fastif(isempty(g.events), 'off', 'on'));
  uimenu('Parent',zm,'Label','Events off'   , 'callback', comnoevent  , 'enable', fastif(isempty(g.events), 'off', 'on'));
  uimenu('Parent',zm,'Label','Events'' legend', 'callback', comeventleg , 'enable', fastif(isempty(g.events), 'off', 'on'));
  
  % %%%%%%%%%%%%%%%%%
  % Set up autoselect
  % %%%%%%%%%%%%%%%%%

  % push button: create/remove window or select electrode
  commandpush = ['ax1 = findobj(''tag'',''backeeg'',''parent'',gcbf);' ... 
			 'tmppos = get(ax1, ''currentpoint'');' ...
  			 'g = get(gcbf,''UserData'');' ... % get data of backgroung image {g.trialstag g.winrej incallback}
             'if g.incallback ~= 1' ... % interception of nestest calls
 			 '   if g.trialstag ~= -1,' ...
			 '   	lowlim = round(g.time*g.trialstag+1);' ...
 			 '   else,' ...
			 '   	lowlim = round(g.time*g.srate+1);' ...
  			 '   end;' ...
			 '  if isempty(g.winrej) Allwin=0;' ...
			 '  else Allwin = (g.winrej(:,1) < lowlim+tmppos(1)) & (g.winrej(:,2) > lowlim+tmppos(1));' ...
			 '  end;' ...
			 '  if any(Allwin)' ... % remove the mark or select electrode if necessary
			 '    lowlim = find(Allwin==1);' ...
 			 '    if g.setelectrode' ...  % select electrode  
 			 '      ax2 = findobj(''tag'',''eegaxis'',''parent'',gcbf);' ...
			 '      tmppos = get(ax2, ''currentpoint'');' ...
    		 '      tmpelec = g.chans + 1 - round(tmppos(1,2) / g.spacing);' ...
    		 '      tmpelec = min(max(tmpelec, 1), g.chans);' ...
			 '      g.winrej(lowlim,tmpelec+5) = ~g.winrej(lowlim,tmpelec+5);' ... % set the electrode
			 '    else' ...  % remove mark
			 ... % '      for tmpi = lowlim''' ...
			 ... %'          tmpcolor = g.winrej(tmpi,3)+10*g.winrej(tmpi,4)+100*g.winrej(tmpi,5);' ...
		     ... %'          if any(tmpcolor == g.colmodif);' ...
			 ... %'              g.winrej(tmpi,:) = []; ' ...
			 ... %'          end;' ...
			 ... %'       end;' ... % THIS PART WAS TO DESELECT COLOR SELECTIVELLY 
             '       g.winrej(lowlim,:) = [];' ...
             '    end;' ...
			 '  else' ...
			 '    if g.trialstag ~= -1' ... % find nearest trials boundaries if epoched data
			 '      alltrialtag = [0:g.trialstag:g.frames];' ...
			 '      I1 = find(alltrialtag < (tmppos(1)+lowlim) );' ... 
			 '      if ~isempty(I1) & I1(end) ~= length(alltrialtag),' ...
			 '            g.winrej = [g.winrej'' [alltrialtag(I1(end)) alltrialtag(I1(end)+1) g.wincolor zeros(1,g.chans)]'']'';' ...
			 '      end;' ...
 			 '    else,' ...
			 '	    g.incallback = 1;' ...  % set this variable for callback for continuous data
			 '      g.winrej = [g.winrej'' [tmppos(1)+lowlim tmppos(1)+lowlim g.wincolor zeros(1,g.chans)]'']'';' ...
			 '    end;' ...
			 '  end;' ...
  			 '  set(gcbf,''UserData'', g);' ...
			 '  eegplot(''drawp'', 0);' ...  % redraw background
             'end;' ...
             'clear g hhdat hh tmpelec tmppos ax2 ESpacing lowlim Allwin Fs winlength EPosition ax1 I1 tmpi' ];
			 		
  % motion button: move windows or display current position (channel, g.time and activation)
  commandmove = ['ax1 = findobj(''tag'',''backeeg'',''parent'',gcbf);' ... 
			 'tmppos = get(ax1, ''currentpoint'');' ...
 			 'g = get(gcbf,''UserData'');' ...
    		 'if isstruct(g)' ...      %check if we are dealing with the right window
 			 '   if g.trialstag ~= -1,' ...
			 '   	lowlim = round(g.time*g.trialstag+1);' ...
 			 '   else,' ...
			 '   	lowlim = round(g.time*g.srate+1);' ...
  			 '   end;' ...
			 '   if g.incallback' ...
			 '      g.winrej = [g.winrej(1:end-1,:)'' [g.winrej(end,1) tmppos(1)+lowlim g.winrej(end,3:end)]'']'';' ...
 			 '      set(gcbf,''UserData'', g);' ...
 			 '      eegplot(''drawb'');' ...
 			 '   else' ...
 			 '     hh = findobj(''tag'',''Etime'',''parent'',gcbf);' ...
 			 '     if g.trialstag ~= -1,' ...
 			 '        set(hh, ''string'', num2str(mod(tmppos(1)+lowlim-1,g.trialstag)/g.trialstag*(g.limits(2)-g.limits(1)) + g.limits(1)));' ...
 			 '     else,' ...
  			 '     	  set(hh, ''string'', num2str((tmppos(1)+lowlim-1)/g.srate));' ... % put g.time in the box
  			 '     end;' ...
  			 '     ax1 = findobj(''tag'',''eegaxis'',''parent'',gcbf);' ... 
			 '     tmppos = get(ax1, ''currentpoint'');' ...
    		 '     tmpelec = round(tmppos(1,2) / g.spacing);' ...
    		 '     tmpelec = min(max(tmpelec, 1),g.chans);' ...
    		 '     labls = get(ax1, ''YtickLabel'');' ...
 			 '     hh = findobj(''tag'',''Eelec'',''parent'',gcbf);' ...  % put electrode in the box
 			 '     set(hh, ''string'', labls(tmpelec+1,:));' ...
 			 '     hh = findobj(''tag'',''Evalue'',''parent'',gcbf);' ...
             '     eegplotdata = get(ax1, ''userdata'');' ...
  			 '     set(hh, ''string'', num2str(eegplotdata(g.chans+1-tmpelec, min(g.frames,max(1,round(tmppos(1)+lowlim))))));' ...  % put value in the box
    		 '  end;' ...
			 'end;' ...
			 'clear g labls eegplotdata tmpelec nbdat ESpacing tmppos ax1 hh lowlim Fs winlength;' ];


  % release button: check window consistency, adpat to trial boundaries
  commandrelease = ['ax1 = findobj(''tag'',''backeeg'',''parent'',gcbf);' ... 
 			 'g = get(gcbf,''UserData'');' ...
 			 'g.incallback = 0;' ...
			 'set(gcbf,''UserData'', g); ' ... % early save in case of bug in the following
			 'if ~isempty(g.winrej)', ...
			 '	if g.winrej(end,1) == g.winrej(end,2)' ... % remove unitary windows
			 '		g.winrej = g.winrej(1:end-1,:);' ...
			 '  else' ...
             '    if g.winrej(end,1) > g.winrej(end,2)' ... % reverse values if necessary
             '       g.winrej(end, 1:2) = [g.winrej(end,2) g.winrej(end,1)];' ...
             '    end;' ...
			 '    g.winrej(end,1) = max(1, g.winrej(end,1));' ...
			 '    g.winrej(end,2) = min(g.frames, g.winrej(end,2));' ...
             '    if g.trialstag == -1' ... % find nearest trials boundaries if necessary
			 '      I1 = find((g.winrej(end,1) >= g.winrej(1:end-1,1)) & (g.winrej(end,1) <= g.winrej(1:end-1,2)) );' ... 
			 '      if ~isempty(I1)' ...
			 '		    g.winrej(I1,2) = max(g.winrej(I1,2), g.winrej(end,2));' ... % extend epoch
			 '		    g.winrej = g.winrej(1:end-1,:);' ... % remove if empty match
			 '      else,' ...
			 '          I2 = find((g.winrej(end,2) >= g.winrej(1:end-1,1)) & (g.winrej(end,2) <= g.winrej(1:end-1,2)) );' ... 
			 '          if ~isempty(I2)' ...
			 '		        g.winrej(I2,1) = min(g.winrej(I2,1), g.winrej(end,1));' ... % extend epoch
			 '		        g.winrej = g.winrej(1:end-1,:);' ... % remove if empty match
			 '		    else,' ...
			 '              I2 = find((g.winrej(end,1) <= g.winrej(1:end-1,1)) & (g.winrej(end,2) >= g.winrej(1:end-1,1)) );' ... 
			 '              if ~isempty(I2)' ...
			 '		            g.winrej(I2,:) = [];' ... % remove if empty match
			 '              end;' ...
			 '		    end;' ...
			 '		end;' ...
			 '    end;' ...
			 '  end;' ...
			 'end;' ...
          'set(gcbf,''UserData'', g);' ...
          'eegplot(''drawp'', 0);' ...
          'clear alltrialtag g tmptmp ax1 I1 I2 trialtag hhdat hh;'];

  set(figh, 'windowbuttondownfcn', commandpush);
  set(figh, 'windowbuttonmotionfcn', commandmove);
  set(figh, 'windowbuttonupfcn', commandrelease);
  set(figh, 'interruptible', 'off');
  set(figh, 'busyaction', 'cancel');
  
  % prepare event array if any
  % --------------------------
  if ~isempty(g.events)
      if ~isfield(g.events, 'type') | ~isfield(g.events, 'latency'), g.events = []; end;
  end;
      
  if ~isempty(g.events)
      if isstr(g.events(1).type)
           [g.eventtypes tmpind indexcolor] = unique({g.events.type});
      else [g.eventtypes tmpind indexcolor] = unique(cell2mat({g.events.type}));
      end;
      g.eventcolors     = { 'r', [0 0.8 0], 'm', 'c', 'k', 'b', [0 0.8 0] };  
      g.eventstyle      = { '-' '-' '-'  '-'  '-' '-' '-' '--' '--' '--'  '--' '--' '--' '--'};  
      g.eventtypecolors = g.eventcolors(mod([1:length(g.eventtypes)]-1 ,length(g.eventcolors))+1);
      g.eventcolors     = g.eventcolors(mod(indexcolor-1               ,length(g.eventcolors))+1);
      g.eventtypestyle  = g.eventstyle (mod([1:length(g.eventtypes)]-1 ,length(g.eventstyle))+1);
      g.eventstyle      = g.eventstyle (mod(indexcolor-1               ,length(g.eventstyle))+1);
      g.eventlatencies  = cell2mat({g.events.latency});
      g.plotevent       = 'on';
  end;
  if isempty(g.events)
      g.plotevent      = 'off';
  end;
  g.commandselect = { commandpush commandmove commandrelease };           
  set(figh, 'userdata', g);
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot EEG Data
  % %%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ax1)
  hold on
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot Spacing I
  % %%%%%%%%%%%%%%%%%%%%%%%%%%
  YLim = get(ax1,'Ylim');
  A = DEFAULT_AXES_POSITION;
  axes('Position',[A(1)+A(3) A(2) 1-A(1)-A(3) A(4)],'Visible','off','Ylim',YLim,'tag','eyeaxes')
  axis manual
  if strcmp(SPACING_EYE,'on'),  set(m(7),'checked','on')
  else set(m(7),'checked','off');
  end 
  eegplot('scaleeye', [], gcf);
  if strcmp(lower(g.scale), 'off')
	  eegplot('scaleeye', 'off', gcf);
  end;
  
  eegplot('drawp', 0);
  eegplot('drawp', 0);
  if g.dispchans ~= g.chans
  	   eegplot('zoom', gcf);
  end;  
  eegplot('scaleeye', [], gcf);
  
  h = findobj(gcf, 'style', 'pushbutton');
  set(h, 'backgroundcolor', BUTTON_COLOR);
  h = findobj(gcf, 'tag', 'eegslider');
  set(h, 'backgroundcolor', BUTTON_COLOR);
  set(figh, 'visible', 'on');
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Main Function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
  try, p1 = varargin{1}; p2 = varargin{2}; p3 = varargin{3}; catch, end;
  switch data
  case 'drawp' % Redraw EEG and change position

    % this test help to couple eegplot windows
    if exist('p3')
    	figh = p3;
    	figure(p3);
    else	
    	figh = gcf;                          % figure handle
	end;
	
    if strcmp(get(figh,'tag'),'dialog')
      figh = get(figh,'UserData');
    end
    ax0 = findobj('tag','backeeg','parent',figh); % axes handle
    ax1 = findobj('tag','eegaxis','parent',figh); % axes handle
    g = get(figh,'UserData');
    data = get(ax1,'UserData');
    ESpacing = findobj('tag','ESpacing','parent',figh);   % ui handle
    EPosition = findobj('tag','EPosition','parent',figh); % ui handle
    if g.trialstag(1) == -1
        g.time    = str2num(get(EPosition,'string'));  
    else 
        g.time    = str2num(get(EPosition,'string'));
        g.time    = g.time - 1;
    end; 
    g.spacing = str2num(get(ESpacing,'string'));
        
    if p1 == 1
      g.time = g.time-g.winlength;     % << subtract one window length
    elseif p1 == 2               
      g.time = g.time-fastif(g.winlength>=1, 1, g.winlength/5);             % < subtract one second
    elseif p1 == 3
      g.time = g.time+fastif(g.winlength>=1, 1, g.winlength/5);             % > add one second
    elseif p1 == 4
      g.time = g.time+g.winlength;     % >> add one window length
    end
    
	 if g.trialstag ~= -1 % time in second or in trials
		multiplier = g.trialstag;	
	 else
		multiplier = g.srate;
	 end;		
    
    % Update edit box
    % ---------------
    g.time = max(0,min(g.time,ceil((g.frames-1)/multiplier)-g.winlength));
    if g.trialstag(1) == -1
        set(EPosition,'string',num2str(g.time)); 
    else 
        set(EPosition,'string',num2str(g.time+1)); 
    end; 
    set(figh, 'userdata', g);

    lowlim = round(g.time*multiplier+1);
    highlim = round(min((g.time+g.winlength)*multiplier+2,g.frames));
    
    % Plot data and update axes
    % -------------------------
    switch lower(g.submean) % subtract the mean ?
     case 'on', meandata = mean(data(:,lowlim:highlim)');  
     case 'nan',meandata = nan_mean(data(:,lowlim:highlim)');
     otherwise, meandata = zeros(1,g.chans);
    end;
    axes(ax1)
    cla
     
    % plot data
    % ---------
    axes(ax1)
    hold on
    for i = 1:g.chans
        plot(data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing, ...
             'color', g.color{mod(i-1,length(g.color))+1}, 'clipping','on')
    end
     
    % draw selected channels
    % ------------------------
    if ~isempty(g.winrej)
    	for tpmi = 1:size(g.winrej,1) % scan rows
			if (g.winrej(tpmi,1) >= lowlim & g.winrej(tpmi,1) <= highlim) | ...
				(g.winrej(tpmi,2) >= lowlim & g.winrej(tpmi,2) <= highlim)
				abscmin = max(1,round(g.winrej(tpmi,1)-lowlim));	 
				abscmax = round(g.winrej(tpmi,2)-lowlim);
				maxXlim = get(gca, 'xlim');
				abscmax = min(abscmax, round(maxXlim(2)-1));	 
   				for i = 1:g.chans
   					if g.winrej(tpmi,g.chans-i+1+5)
   						plot(abscmin+1:abscmax+1,data(g.chans-i+1,abscmin+lowlim:abscmax+lowlim) ...
   							-meandata(g.chans-i+1)+i*g.spacing, 'color','r','clipping','on')
					end;
    			end
  			end;	
    	end;
    end;		
    
    set(ax1, 'Xlim',[1 g.winlength*multiplier+1],...
		     'XTick',[1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1]);
    set(ax1, 'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'))

    % ordinates: even if all elec are plotted, some may be hidden
    set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
    eegplot('drawb');
	
	 if g.children ~= 0
		if ~exist('p2')
			p2 =[];
		end;	
		eegplot( 'drawp', p1, p2, g.children);
		figure(figh);
	 end;	  
  
  case 'drawb' % Draw background ******************************************************
    % Redraw EEG and change position

    ax0 = findobj('tag','backeeg','parent',gcf); % axes handle
    ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
        
    g = get(gcf,'UserData');  % Data (Note: this could also be global)

    % Plot data and update axes
    axes(ax0);
	cla;
	hold on;
	% plot rejected windows
	if g.trialstag ~= -1
		multiplier = g.trialstag;	
	else
		multiplier = g.srate;
	end;		
   	lowlim = round(g.time*multiplier+1);
   	highlim = round(min((g.time+g.winlength)*multiplier+1));
  	displaymenu = findobj('tag','displaymenu','parent',gcf);
    if ~isempty(g.winrej) & g.winstatus
		if g.trialstag ~= -1 % epoched data
			indices = find((g.winrej(:,1)' >= lowlim & g.winrej(:,1)' <= highlim) | ...
						   (g.winrej(:,2)' >= lowlim & g.winrej(:,2)' <= highlim));
			if ~isempty(indices)
				tmpwins1 = g.winrej(indices,1)';
				tmpwins2 = g.winrej(indices,2)';
				tmpcols  = g.winrej(indices,3:5);
				try, eval('[cumul indicescount] = histc(tmpwins1, (min(tmpwins1)-1):g.trialstag:max(tmpwins2));');
				catch, [cumul indicescount] = myhistc(tmpwins1, (min(tmpwins1)-1):g.trialstag:max(tmpwins2));
				end;
				count = zeros(size(cumul));
				%if ~isempty(find(cumul > 1)), find(cumul > 1), end;
                for tmpi = 1:length(tmpwins1)
					poscumul = indicescount(tmpi);
					heightbeg = count(poscumul)/cumul(poscumul);
					heightend = heightbeg + 1/cumul(poscumul);
					count(poscumul) = count(poscumul)+1;
					h = patch([tmpwins1(tmpi)-lowlim tmpwins2(tmpi)-lowlim ...
							   tmpwins2(tmpi)-lowlim tmpwins1(tmpi)-lowlim], ...
							  [heightbeg heightbeg heightend heightend], ...
							  tmpcols(tmpi,:));  % this argument is color
					set(h, 'EdgeColor', get(h, 'facecolor')) 
				end;
			end;
		else
			for tpmi = 1:size(g.winrej,1) % scan rows
				if (g.winrej(tpmi,1) >= lowlim & g.winrej(tpmi,1) <= highlim) | ...
						(g.winrej(tpmi,2) >= lowlim & g.winrej(tpmi,2) <= highlim)	 
					h = patch([g.winrej(tpmi,1)-lowlim g.winrej(tpmi,2)-lowlim ...
							   g.winrej(tpmi,2)-lowlim g.winrej(tpmi,1)-lowlim], ...
							  [0 0 1 1], [g.winrej(tpmi,3) g.winrej(tpmi,4) g.winrej(tpmi,5)]);  
					set(h, 'EdgeColor', get(h, 'facecolor')) 
				end;	
			end;
		end;
	end;
    		
	% plot tags
	% ---------
    %if trialtag(1) ~= -1 & displaystatus % put tags at arbitrary places
    % 	for tmptag = trialtag
	%		if tmptag >= lowlim & tmptag <= highlim
	%			plot([tmptag-lowlim tmptag-lowlim], [0 1], 'b--');
   	%		end;	
    %	end;
    %end;

    % draw events if any
    % ------------------
    if strcmpi(g.plotevent, 'on')
        event2plot = find ( g.eventlatencies >=lowlim & g.eventlatencies <= highlim );
        for index = 1:length(event2plot)
            tmplat = g.eventlatencies(event2plot(index))-lowlim-1;
            tmph   = plot([ tmplat tmplat ], ylim, 'color', g.eventcolors{ event2plot(index) }, ...
                     'linestyle', g.eventstyle{ event2plot(index) } );
            %if g.trials == 1
            %    set(tmph, 'userdata', sprintf('Type: %s; Lat: %.4f s');
            %else
            %    set(tmph, 'userdata', sprintf('Type: %s; Lat: %.2f ms');
            %end;
        end;
    end;

    if g.trialstag(1) ~= -1
        
        % plot trial limits
        % -----------------
        tmptag = [lowlim:highlim];
       	tmpind = find(mod(tmptag-1, g.trialstag) == 0);
        for index = tmpind
            plot([tmptag(index)-lowlim-1 tmptag(index)-lowlim-1], [0 1], 'b--');
        end;
        alltag = tmptag(tmpind);

        % compute Xticks
        % --------------
        tagnum = (alltag-1)/g.trialstag+1;
     	set(ax0,'XTickLabel', tagnum,'YTickLabel', [],...
		'Xlim',[0 g.winlength*multiplier],...
		'XTick',alltag-lowlim+g.trialstag/2, 'YTick',[], 'tag','backeeg');
		
		axes(ax1);
		tagpos  = [];
		tagtext = [];
		if ~isempty(alltag)
			alltag = [alltag(1)-g.trialstag alltag alltag(end)+g.trialstag]; % add border trial limits
		else
			alltag = [ floor(lowlim/g.trialstag)*g.trialstag ceil(highlim/g.trialstag)*g.trialstag ]+1;
		end;
        
		nbdiv = 20/g.winlength; % approximative number of divisions
		divpossible = [ 100000./[1 2 4 5] 10000./[1 2 4 5] 1000./[1 2 4 5] 100./[1 2 4 5 10 20]]; % possible increments
		[tmp indexdiv] = min(abs(nbdiv*divpossible-(g.limits(2)-g.limits(1)))); % closest possible increment

		incrementpoint = divpossible(indexdiv)/1000*g.srate;
		tagzerooffset  = -g.limits(1)/1000*g.srate; 

		for i=1:length(alltag)-1
			if ~isempty(tagpos) & tagpos(end)-alltag(i)<2*incrementpoint/3
				tagpos  = tagpos(1:end-1);
			end;
			if ~isempty(g.freqlimits)
				tagpos  = [ tagpos linspace(alltag(i),alltag(i+1)-1, nbdiv) ];
			else
				if tagzerooffset ~= 0
					tmptagpos = [alltag(i)+tagzerooffset:-incrementpoint:alltag(i)];
				else
					tmptagpos = [];
				end;
				tagpos  = [ tagpos [tmptagpos(end:-1:2) alltag(i)+tagzerooffset:incrementpoint:(alltag(i+1)-1)]];
			end;
		end;
        % find corresponding epochs
        % -------------------------
        tagtext = eeg_point2lat(tagpos, floor((tagpos)/g.trialstag)+1, g.srate, g.limits, 1E-3);
     	set(ax1,'XTickLabel', tagtext,'XTick', tagpos-lowlim);
			
    else
     	set(ax0,'XTickLabel', [],'YTickLabel', [],...
		'Xlim',[0 g.winlength*multiplier],...
		'XTick',[], 'YTick',[], 'tag','backeeg');

		axes(ax1);
    	set(ax1,'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'),...
		'XTick',[1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1])
    end;
    		
    % ordinates: even if all elec are plotted, some may be hidden
    set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
    
    axes(ax1)	

  case 'draws'
    % Redraw EEG and change scale

    ax1 = findobj('tag','eegaxis','parent',gcf);         % axes handle
    g = get(gcf,'UserData');  
    data = get(ax1, 'userdata');
    ESpacing = findobj('tag','ESpacing','parent',gcf);   % ui handle
    EPosition = findobj('tag','EPosition','parent',gcf); % ui handle
    if g.trialstag(1) == -1
        g.time    = str2num(get(EPosition,'string'));  
    else 
        g.time    = str2num(get(EPosition,'string'))-1;   
    end;        
    g.spacing = str2num(get(ESpacing,'string'));  
    
    orgspacing= g.spacing;
    if p1 == 1
      	g.spacing= g.spacing+ 0.1*orgspacing; % increase g.spacing(5%)
	elseif p1 == 2
		g.spacing= max(0,g.spacing-0.1*orgspacing); % decrease g.spacing(5%)
    end
    if round(g.spacing*100) == 0
        maxindex = min(10000, g.frames);  
        g.spacing = 0.01*max(max(data(:,1:maxindex),[],2),[],1)-min(min(data(:,1:maxindex),[],2),[],1);  % Set g.spacingto max/min data
    end;

    set(ESpacing,'string',num2str(g.spacing,4))  % update edit box
    set(gcf, 'userdata', g);
	 eegplot('drawp', 0);
    set(ax1,'YLim',[0 (g.chans+1)*g.spacing],'YTick',[0:g.spacing:g.chans*g.spacing])
    set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
    
    % update scaling eye if it exists
    eyeaxes = findobj('tag','eyeaxes','parent',gcf);
    if ~isempty(eyeaxes)
      eyetext = findobj('type','text','parent',eyeaxes,'tag','thescalenum');
      set(eyetext,'string',num2str(g.spacing,4))
    end
    
    return;

  case 'window'  % change window size
    % get new window length with dialog box
    g = get(gcf,'UserData');
	result       = inputdlg2( { fastif(g.trialstag==-1,'Enter new window length(secs):', 'Enter number of epoch(s):') }, 'Change window length', 1,  { num2str(g.winlength) });
	if size(result,1) == 0 return; end;

	g.winlength = eval(result{1}); 
	set(gcf, 'UserData', g);
	eegplot('drawp',0);	
	return;
    
  case 'winelec'  % change channel window size
   % get new window length with dialog box
   fig = gcf;
   g = get(gcf,'UserData');
   result = inputdlg2( { 'Enter number of channels to show:' } , 'Change number of channels to show', 1,  { num2str(g.dispchans) });
   if size(result,1) == 0 return; end;

   g.dispchans = eval(result{1});
   if g.dispchans<0 | g.dispchans>g.chans
      g.dispchans =g.chans;
   end;
	set(gcf, 'UserData', g);
   eegplot('updateslidder', fig);
	eegplot('drawp',0);	
	return;
   
  case 'loadelect' % load channels
	[inputname,inputpath] = uigetfile('*','Channel locations file');
	if inputname == 0 return; end;
	if ~exist([ inputpath inputname ])
		error('no such file');
	end;

	AXH0 = findobj('tag','eegaxis','parent',gcf);
	eegplot('setelect',[ inputpath inputname ],AXH0);
	return;
  
  case 'setelect'
    % Set channels    
    eloc_file = p1;
    axeshand = p2;
    outvar1 = 1;
    if isempty(eloc_file)
      outvar1 = 0;
      return
    end
    
	[tmp YLabels] = readlocs(eloc_file);
    YLabels = strvcat(YLabels);
    
    YLabels = flipud(str2mat(YLabels,' '));
    set(axeshand,'YTickLabel',YLabels)
  
  case 'title'
    % Get new title
	h = findobj('tag', 'eegplottitle');
	
	if ~isempty(h)
		result       = inputdlg2( { 'Enter new title:' }, 'Change title', 1,  { get(h(1), 'string') });
		if ~isempty(result), set(h, 'string', result{1}); end;
	else 
		result       = inputdlg2( { 'Enter new title:' }, 'Change title', 1,  { '' });
		if ~isempty(result), h = textsc(result{1}, 'title'); set(h, 'tag', 'eegplottitle');end;
	end;
	
	return;

  case 'scaleeye'
    % Turn scale I on/off
    obj = p1;
    figh = p2;
	g = get(figh,'UserData');
    % figh = get(obj,'Parent');

    if ~isempty(obj)
		eyeaxes = findobj('tag','eyeaxes','parent',figh);
		children = get(eyeaxes,'children');
		if isstr(obj)
			if strcmp(obj, 'off')
				set(children, 'visible', 'off');
				set(eyeaxes, 'visible', 'off');
				return;
			else
				set(children, 'visible', 'on');
				set(eyeaxes, 'visible', 'on');
			end;
		else
			toggle = get(obj,'checked');
			if strcmp(toggle,'on')
				set(children, 'visible', 'off');
				set(eyeaxes, 'visible', 'off');
				set(obj,'checked','off');
				return;
			else
				set(children, 'visible', 'on');
				set(eyeaxes, 'visible', 'on');
				set(obj,'checked','on');
			end;
		end;
	end;
	
	eyeaxes = findobj('tag','eyeaxes','parent',figh);
    ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
	YLim = get(ax1, 'ylim');
    
	ESpacing = findobj('tag','ESpacing','parent',figh);
	g.spacing= str2num(get(ESpacing,'string'));
	
	axes(eyeaxes); cla; axis off;
    set(eyeaxes, 'ylim', YLim);
    
	Xl = [.35 .65; .5 .5; .35 .65];
    Yl = [ g.spacing g.spacing; g.spacing 0; 0 0] + YLim(1);
	plot(Xl(1,:),Yl(1,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline'); hold on;
	plot(Xl(2,:),Yl(2,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline');
	plot(Xl(3,:),Yl(3,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline');
	text(.5,(YLim(2)-YLim(1))/23+Yl(1),num2str(g.spacing,4),...
		 'HorizontalAlignment','center','FontSize',10,...
		 'tag','thescalenum')
    text(Xl(2)+.1,Yl(1),'+','HorizontalAlignment','left',...
         'verticalalignment','middle', 'tag', 'thescale')
    text(Xl(2)+.1,Yl(4),'-','HorizontalAlignment','left',...
         'verticalalignment','middle', 'tag', 'thescale')
	if ~isempty(SPACING_UNITS_STRING)
        text(.5,-YLim(2)/23+Yl(4),SPACING_UNITS_STRING,...
			 'HorizontalAlignment','center','FontSize',10, 'tag', 'thescale')
	end
	text(.5,(YLim(2)-YLim(1))/10+Yl(1),'Scale',...
		 'HorizontalAlignment','center','FontSize',10, 'tag', 'thescale')
    set(eyeaxes, 'tag', 'eyeaxes');
    
  case 'noui'
      if ~isempty(varargin)
          eegplot( varargin{:} ); fig = gcf;
      else 
          fig = findobj('tag', 'EEGPLOT');
      end;
      set(fig, 'menubar', 'figure');
      
      % find button and text
      obj = findobj(fig, 'style', 'pushbutton'); delete(obj);
      obj = findobj(fig, 'style', 'edit'); delete(obj);
      obj = findobj(fig, 'style', 'text'); 
      %objscale = findobj(obj, 'tag', 'thescale');
      %delete(setdiff(obj, objscale));
	  obj = findobj(fig, 'tag', 'Eelec');delete(obj);
	  obj = findobj(fig, 'tag', 'Etime');delete(obj);
	  obj = findobj(fig, 'tag', 'Evalue');delete(obj);
	  obj = findobj(fig, 'tag', 'Eelecname');delete(obj);
	  obj = findobj(fig, 'tag', 'Etimename');delete(obj);
	  obj = findobj(fig, 'tag', 'Evaluename');delete(obj);
	  obj = findobj(fig, 'type', 'uimenu');delete(obj);
 
   case 'zoom' % if zoom
      fig = varargin{1};
      ax1 = findobj('tag','eegaxis','parent',fig); 
      ax2 = findobj('tag','backeeg','parent',fig); 
      tmpxlim  = get(ax1, 'xlim');
      tmpylim  = get(ax1, 'ylim');
      tmpxlim2 = get(ax2, 'xlim');
      set(ax2, 'xlim', get(ax1, 'xlim'));
      g = get(fig,'UserData');
      
      % deal with abscicia
      % ------------------
      if g.trialstag ~= -1
          Eposition = str2num(get(findobj('tag','EPosition','parent',fig), 'string'));
          g.winlength = (tmpxlim(2) - tmpxlim(1))/g.trialstag;
          Eposition = Eposition + (tmpxlim(1) - tmpxlim2(1)-1)/g.trialstag;
          Eposition = round(Eposition*1000)/1000;
          set(findobj('tag','EPosition','parent',fig), 'string', num2str(Eposition));
      else
          Eposition = str2num(get(findobj('tag','EPosition','parent',fig), 'string'))-1;
          g.winlength = (tmpxlim(2) - tmpxlim(1))/g.srate;	
          Eposition = Eposition + (tmpxlim(1) - tmpxlim2(1)-1)/g.srate;
          Eposition = round(Eposition*1000)/1000;
          set(findobj('tag','EPosition','parent',fig), 'string', num2str(Eposition+1));
      end;  
      
      % deal with ordinate
      % ------------------
      g.elecoffset = tmpylim(1)/g.spacing;
      g.dispchans  = round(1000*(tmpylim(2)-tmpylim(1))/g.spacing)/1000;      
      
      set(fig,'UserData', g);
      eegplot('updateslidder', fig);
      eegplot('drawp', 0);
      eegplot('scaleeye', [], fig);

      % reactivate zoom if 3 arguments
      % ------------------------------
      if exist('p2') == 1
          zoom on;
          tmpstr = get(gcbf, 'windowbuttondownfcn');
          set(gcbf, 'windowbuttondownfcn', [ tmpstr '; eegplot(''zoom'', gcbf, 1);' ]);
          set(gcbf, 'windowbuttonmotionfcn', g.commandselect{2});
      end;

	case 'updateslidder' % if zoom
      fig = varargin{1};
      g = get(fig,'UserData');
      slidder = findobj('tag','eegslider','parent',fig);
      if g.elecoffset < 0
         g.elecoffset = 0;
      end;
      if g.dispchans >= g.chans
         g.dispchans = g.chans;
         g.elecoffset = 0;
         set(slidder, 'visible', 'off');
      else
         set(slidder, 'visible', 'on');         
		 set(slidder, 'value', g.elecoffset/g.chans, ...
					  'sliderstep', [1/(g.chans-g.dispchans) g.dispchans/(g.chans-g.dispchans)]);
         %'sliderstep', [1/(g.chans-1) g.dispchans/(g.chans-1)]);
      end;
      if g.elecoffset < 0
         g.elecoffset = 0;
      end;
      if g.elecoffset > g.chans-g.dispchans
         g.elecoffset = g.chans-g.dispchans;
      end;
      set(fig,'UserData', g);
	  eegplot('scaleeye', [], fig);
   
   case 'drawlegend'
      fig = varargin{1};
      g = get(fig,'UserData');
      
      if ~isempty(g.events)
          nleg = length(g.eventtypes);
          fig2 = figure('numbertitle', 'off', 'name', '', 'visible', 'off', 'menubar', 'none', 'color', DEFAULT_FIG_COLOR);
          pos = get(fig2, 'position');
          set(fig2, 'position', [ pos(1) pos(2) 130 10*nleg+20]);
          
          for index = 1:nleg
              plot([10 30], [(index-0.5) * 10 (index-0.5) * 10], 'color', g.eventtypecolors{index}, 'linestyle', ...
                          g.eventtypestyle{ index }); hold on;
              text(35, (index-0.5)*10, g.eventtypes{index});
          end;
          xlim([0 130]);
          ylim([0 nleg*10]);
          axis off;
          set(fig2, 'visible', 'on');
      end;
          
   otherwise
      error(['Error - invalid eegplot() parameter: ',data])
  end  
end

% function not supported under Mac
% --------------------------------
function [reshist, allbin] = myhistc(vals, intervals);

	reshist = zeros(1, length(intervals));
	allbin = zeros(1, length(vals));
	
	for index=1:length(vals)
		minvals = vals(index)-intervals;
		bintmp  = find(minvals >= 0);
		[mintmp indextmp] = min(minvals(bintmp));
		bintmp = bintmp(indextmp);
		
		allbin(index) = bintmp;
		reshist(bintmp) = reshist(bintmp)+1;
	end;
