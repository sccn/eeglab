% topoplot() - plot a topographic map of a scalp data field in a 2-D circular view 
%              (looking down at the top of the head) using interpolation on a fine 
%              cartesian grid. Can also show specified channnel location(s), or return 
%              an interpolated value at an arbitrary scalp location (see 'noplot').
%              Using option 'gridplot', the plot can be a rectangular imagesc() grid
%              or can add such a grid plot to the left or right of the head image.
%              By default, channel locations below head center (arc_length 0.5) are 
%              shown in a 'skirt' outside the cartoon head (see 'plotrad' and 'headrad' 
%              options below). Nose is at top of plot; left is left; right is right.
% Usage:
%        >>  topoplot(datavector, EEG.chanlocs);   % use a channel locations structure
%        >>  topoplot(datavector, 'my_chan.locs'); % read a channel locations file
%        >>  topoplot('example'); % Gives an example of an electrode location file
%        >>  [h grid_or_val [grid]]= topoplot(datavector, chan_locs, 'Param1','Value1', ...);
%
% Required Inputs:
%   datavector        - single vector of channel values. Else, if a vector of selected 
%                       channel numbers -> mark their location(s) using 'style' 'blank'.
%   chan_locs         - name of an EEG electrode position file (>> topoplot example).
%                       Else, an EEG.chanlocs structure (>> help pop_editset)
% Optional inputs:
%   'maplimits'       - 'absmax'   -> scale map colors to +/- the absolute-max (makes green 0); 
%                       'maxmin'   -> scale colors to the data range (makes green mid-range); 
%                       [lo.hi]    -> use user-definined lo/hi limits {default: 'absmax'}
%   'style'           - 'map'      -> plot colored map only
%                       'contour'  -> plot contour lines only
%                       'both'     -> plot both colored map and contour lines
%                       'fill'     -> plot constant color between contour lines
%                       'blank'    -> plot electrode locations only {default: 'both'}
%   'electrodes'      - 'on','off','labels','numbers','ptslabels','ptsnumbers' See Plot detail 
%                       options below. {default: 'on' -> mark electrode locations with points
%                       unless more than 64 channels, then 'off'}. 
%   'plotchans'       - vector of channel indices to interpolate in the head plot. Negative
%                       integers reverse the data polarity. Grid chans, if any, are not included 
%                       in plotchans (see 'gridplot' below). {default: [] -> plot all chans}
%   'plotrad'         - [0.15<=float<=1.0] plotting radius = max channel arc_length to plot.
%                       See >> topoplot example. If plotrad > 0.5, chans with arc_length > 0.5 
%                       (i.e. below ears-eyes) are plotted in a circular 'skirt' outside the
%                       cartoon head. See 'intrad' below. {default: max(max(chanlocs.radius),0.5)}. 
%   'headrad'         - [0.15<=float<=1.0] drawing radius (arc_length) for the cartoon head. 
%                       NOTE: Only headrad = 0.5 is anatomically correct! 0 -> don't draw head; 
%                       'rim' -> show cartoon head at outer edge of the plot {default: 0.5}
%   'intrad'          - [0.15<=float<=1.0] radius of the interpolation area (square or disk, see
%                       'intsquare' below). Interpolate electrodes in this area {default: channel max}
%   'shrink'          - ['on'|'off'|'force'|factor] Deprecated. 'on' -> If max channel arc_length 
%                       > 0.5, shrink electrode coordinates towards vertex to plot all channels
%                       by making max arc_length 0.5. 'force' -> Normalize arc_length 
%                       so the channel max is 0.5. factor -> Apply a specified shrink
%                       factor (range (0,1) = shrink fraction). {default: 'off'}
%   'intsquare'       - ['on'|'off'] 'on' -> Interpolate values at electrodes located in the whole 
%                       square containing the (radius plotrad) plotting disk. 'off' -> Interpolate
%                       values from electrodes shown in the plotting disk only. {default: 'on'}
%   'conv'            - ['on'|'off'] Show map interpolation only out to the convext hull of
%                       the electrode locations to minimize extrapolation.  {default: 'off'}
%   'noplot'          - ['on'|'off'|[rad theta]] do not plot (but return interpolated data).
%                       Else, if [rad theta] are coordinates of a (possibly missing) channel, 
%                       returns interpolated value for channel location.  For more info, 
%                       see >> topoplot 'example' {default: 'off'}
% Dipole plotting:
%   'dipole'          - [xi yi xe ye ze] plot dipole on the top of the scalp map
%                       from coordinate (xi,yi) to coordinates (xe,ye,ze) (dipole head 
%                       model has radius 1). If several rows, plot one dipole per row.
%                       Coordinates returned by dipplot() may be used. Can accept
%                       an EEG.dipfit.model structure (See >> help dipplot).
%                       Ex: ,'dipole',EEG.dipfit.model(17) % Plot dipole(s) for comp. 17.
%   'dipnorm'         - ['on'|'off'] normalize dipole length {default: 'on'}.
%   'diporient'       - [-1|1] invert dipole orientation {default: 1}.
%   'diplen'          - [real] scale dipole length {default: 1}.
%   'dipscale'        - [real] scale dipole size {default: 1}.
%   'dipsphere'       - [real] size of the dipole sphere. {default: 85 mm}.
%   'dipcolor'        - [color] dipole color as Matlab code code or [r g b] vector
%                       {default: 'k' = black}.
% Plot detail options:
%   'electcolor' {'k'}|'emarker' {'.'}|'emarkersize' {14}|'emarkersize1chan' {40}|'efontsize' {var}
%                       electrode marking details and their {defaults}. 
%   'shading'         - 'flat','interp'  {default: 'flat'}
%   'colormap'        -  (n,3) any size colormap {default: existing colormap}
%   'numcontour'      - number of contour lines {default: 6}
%   'ccolor'          - color of the contours {default: dark grey}
%   'hcolor'|'ecolor' - colors of the cartoon head and electrodes {default: black}
%   'gridscale'       - [int > 32] size (nrows) of interpolated data matrix {default: 67}
%   'circgrid'        - [int > 100] number of elements (angles) in head and border circles {201}
%   'verbose'         - ['on'|'off'] comment on operations on command line {default: 'on'}.
%
% Outputs:
%         h           - plot axes handle
%         grid_or_val - [matrix] the interpolated data image (with off-head points = NaN).  
%                       ELSE, single interpolated value at the specified 'noplot' arg channel 
%                       location ([rad theta]), if any.
%         [grid]      - IF interpolated value returned above, then the interpolated image grid
% Chan_locs format:
%    See >> topoplot 'example'
%
% Notes: - To change the plot map masking ring to a new figure background color,
%            >> set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
%        - Topoplots may be rotated. From the commandline >> view([deg 90]) {default: [0 90])
%
% Authors: Andy Spydell, Colin Humphries, Arnaud Delorme & Scott Makeig
%          CNL / Salk Institute, 8/1996-/10/2001; SCCN/INC/UCSD, Nov. 2001 -
%
% See also: timtopo(), envtopo()

% Unimplemented future options:
%   'plotgrid'        - [channels] or {[channels], position} where channels is a matrix of grid 
%                       channel numbers - in which 0s plot 0-values and negative integers, 
%                       polarity-reversed values - and char position (if either 'l' or 'r') 
%                       indicates which side of the head to plot the grid, or if 'o', plot the 
%                       grid only {the default grid location is to the left of the head}. The 
%                       channels matrix should show the topographic ordering of the channels. 
%                       {default: no grid plot} NOTE: Not yet implemented.
%                       Ex: >> topoplot(values,'chanlocs','plotgrid',{[11 12 0; 13 14 15],'o'});
%                       % Plot a 2x3 grid (only) of channels 11-15 data values plus one 0 value.

% Deprecated but still usable;
%   'interplimits'    - ['electrodes'|'head'] 'electrodes'-> interpolate the electrode grid; 
%                       'head'-> interpolate the whole disk {default: 'head'}.

% Copyright (C) Colin Humphries & Scott Makeig, CNL / Salk Institute, Aug, 1996
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
% Revision 1.221  2004/11/22 05:41:43  arno
% more debugging
%
% Revision 1.220  2004/11/22 05:39:05  arno
% function was crashing on regular topoplot, debuging
%
% Revision 1.219  2004/11/22 05:02:51  scott
% fixing topoplot([],EEG.chanlocs,'emarker','o')
%
% Revision 1.218  2004/11/22 04:58:21  scott
% fixed topoplot(32,EEG.chanlocs) and topoplot([],EEG.chanlocs,'emarker','o')
% to plot marked channel 32 in red disk
%
% Revision 1.217  2004/11/18 20:29:14  hilit
% enabled the 'example' option
%
% Revision 1.216  2004/11/18 19:22:22  scott
% made 3rd output, 'grid'. [] unless interpolated value asked for
%
% Revision 1.215  2004/11/09 19:25:08  arno
% move plotgrid help outsie of header since unimplemented
%
% Revision 1.214  2004/10/27 17:34:35  scott
% help msg adjust -sm
%
% Revision 1.213  2004/10/27 16:39:06  arno
% remove infinite and NaN values
%
% Revision 1.212  2004/10/09 22:26:18  scott
% iv interp. value output, then output grid too -sm
%
% Revision 1.211  2004/10/08 21:34:25  scott
% same -sm
%
% Revision 1.210  2004/10/08 21:32:09  scott
% help message clarification on outputs -sm
%
% Revision 1.209  2004/10/07 15:55:15  scott
% made Values==[] work with plotchans  -sm
%
% Revision 1.208  2004/09/29 15:44:46  scott
% added 'plotchans' option. upgraded 'plotgrid' (still unimplemented) -sm
%
% Revision 1.207  2004/09/29 01:04:22  scott
% created input 'plotgrid' - plotting not yet implemented -sm
%
% Revision 1.206  2004/09/10 00:53:08  hilit
% converted input arguments to text() to double
%
% Revision 1.205  2004/07/07 22:21:30  arno
% debug shrink
%
% Revision 1.204  2004/06/10 19:11:53  arno
% remove debug msg
%
% Revision 1.203  2004/05/14 23:41:09  arno
% allowing negative shrink
%
% Revision 1.202  2004/05/14 23:29:32  arno
% fixing toggle name/number
%
% Revision 1.201  2004/05/10 15:14:34  scott
% more flexible labels/numbers/points argument reading; defined ELECTRODE_HEIGHT
%
% Revision 1.200  2004/05/07 15:12:51  scott
% removed textax, instead plot3() electrode labels/pts/numbers above the rest with plot3() -sm
%
% Revision 1.199  2004/05/07 04:35:10  scott
% superimpose textax again - making both axes square
%
% Revision 1.198  2004/05/05 21:57:23  hilit
% removed text from the previous log message
%
% Revision 1.197  2004/05/05 20:56:20  hilit
% change the defult setting of dipnorm to 'on'.
%
% Revision 1.197  2004/05/05 13:55:00  hilit
% Set the defult option of dipnorm to 'on'
%
% Revision 1.196  2004/05/05 20:36:04  scott
% DIPOLE scaling
%
% Revision 1.195  2004/05/05 20:21:02  scott
% *** empty log message ***
%
% Revision 1.194  2004/04/29 18:58:48  scott
% removed new axes - scaling problem. Toggling pts|numbers doesnt work inside head cartoon.
%
% Revision 1.193  2004/04/29 18:36:46  scott
% test
%
% Revision 1.192  2004/04/29 18:23:03  scott
% make overplot axis limits the same as topoplot limits
%
% Revision 1.191  2004/04/28 18:19:06  scott
% put labels/numbers on another axes so that clicking numbers<->labels
% will work inside the head cartoon patch
%
% Revision 1.190  2004/04/28 17:00:42  scott
% no blanking ring when style is 'blank'
%
% Revision 1.189  2004/04/01 17:10:46  scott
% converted 'conv' interpolation to polar
%
% Revision 1.188  2004/03/31 18:23:15  scott
% debug 'conv' mode - plot ears and nose above map surface to avoid masking by 'conv'
%
% Revision 1.187  2004/03/31 18:06:53  scott
% adding 'conv' mode for plotting convex hull; corrected shrink in 'interp' mode
%
% Revision 1.186  2004/03/31 05:15:05  scott
% *** empty log message ***
%
% Revision 1.185  2004/03/31 05:06:27  scott
% implementing 'conv' (undocumented)
%
% Revision 1.184  2004/03/31 03:19:02  scott
% adjust ear lines
%
% Revision 1.183  2004/03/31 02:53:35  scott
% made blanking ring and head filled rings; made default electrodes 'off' iff chans>64; made contour color
% dark grey; adjusted nose and ear shapes
%
% Revision 1.182  2004/03/31 02:08:07  scott
% *** empty log message ***
%
% Revision 1.181  2004/03/30 18:48:21  scott
% same
%
% Revision 1.180  2004/03/30 18:29:08  scott
% testing fill ring
%
% Revision 1.179  2004/03/30 17:38:15  scott
% plot ring patch instead of blanking circle
%
% Revision 1.178  2004/03/25 22:30:13  arno
% same thing
%
% Revision 1.177  2004/03/25 22:26:45  arno
% same thing
%
% Revision 1.176  2004/03/25 22:24:41  arno
% fixing shrinkfactor bug
%
% Revision 1.175  2004/03/24 16:35:25  scott
% added 'cricgrid' plotting detail argument
%
% Revision 1.174  2004/03/23 19:19:34  scott
% made 'electrodes' default 'off'
%
% Revision 1.173  2004/03/23 19:18:32  scott
% default: plotrad >= 0.5
%
% Revision 1.172  2004/03/23 15:20:39  scott
% made only 2 outputs
%
% Revision 1.171  2004/03/23 00:40:06  scott
% clarifying handling of un-located channels
%
% Revision 1.170  2004/03/22 17:57:21  scott
% added arg 'intrad' - separated interpolation and plotting areas
% Now, by default, interpolates over all the (radius<=1) electrodes.
% Added 'intsquare' option - interpolated values in electrodes in the entire
% interpolation square, not just the (plotting) disk. Can give more accurate
% interpolation at edges of the plotting disk i.e. interpolation instead of
% extrapolation), if there are additional channel locations beyond the plotting area
%
% Revision 1.169  2004/03/22 03:25:41  scott
% re-implmenting shrink options
%
% Revision 1.168  2004/03/21 19:19:18  scott
% help message
%
% Revision 1.167  2004/03/21 18:02:08  scott
% debugged deprecated 'shrink' mode code
%
% Revision 1.166  2004/03/21 17:31:44  scott
% nothing
%
% Revision 1.165  2004/03/21 17:25:39  scott
% corrected dipole plotting
%
% Revision 1.164  2004/03/21 16:52:37  scott
% debugged plotrad, headrad plot size setting
%
% Revision 1.163  2004/03/20 18:20:14  scott
% created 'headrad' (removed 'forcehead'). Now uses only 'plotrad' and 'headrad'
% to set plotting scales. 'shrink' mode disabled temporarily
%
% Revision 1.162  2004/03/19 21:57:58  scott
% do not plot channels with arc_length > 1
%
% Revision 1.161  2004/03/19 19:47:13  arno
% remove str2num
%
% Revision 1.160  2004/03/19 19:05:26  scott
% read string plotrad from channel locations structure
%
% Revision 1.159  2004/03/19 17:46:19  scott
% added 'forcehead'; changed 'pointnumbers' and 'pointlabels' to 'ptnumbers', 'ptlabels'
% but kept backwards compatibility. Allowed marking of multiple channel locations
% without requiring an explicit 'style','blank'. Allowed [] data -> plot channel
% locations. Improved help message and 'example' text. Switched order of plotting
% of head border, electrodes and head cartoon. Made head cartoon not appear by
% default when plotrad<0.5 or 'shrink' is severe (but see 'forcehead'). -sm
%
% Revision 1.158  2004/03/19 02:33:40  scott
% plotting head, ears and/or skirt as appropriate from plotrad and shrink args
%
% Revision 1.157  2004/03/19 01:49:07  scott
% plotrad
%
% Revision 1.156  2004/03/19 00:30:08  scott
% plotrad minmax
%
% Revision 1.155  2004/03/18 17:05:20  arno
% fixed plotrad
%
% Revision 1.154  2004/03/18 16:36:53  arno
% debug shrink and plotrad
%
% Revision 1.153  2004/03/18 16:22:12  arno
% debug shrink
%
% Revision 1.152  2004/03/18 01:47:24  scott
% debug
%
% Revision 1.151  2004/03/18 01:44:28  scott
% 'plotrad' arg and help message re skirt
%
% Revision 1.150  2004/03/18 01:26:33  arno
% plotrad
%
% Revision 1.149  2004/03/18 00:29:07  arno
% debug skirt option
%
% Revision 1.148  2004/03/18 00:18:09  arno
% skirt option
%
% Revision 1.147  2004/02/25 15:29:39  scott
% dont plot border if shrinkfac < .01
%
% Revision 1.146  2004/02/25 15:25:07  scott
% adjust border of 'skirt'
%
% Revision 1.145  2004/02/25 15:19:38  scott
% not allowing shrink to be negative
%
% Revision 1.144  2004/02/23 16:55:51  scott
% don't let ears go outside axes if shrink is 'skirt' but shrink factor is 0 or small
%
% Revision 1.143  2004/02/19 15:56:28  scott
% plot dipole(s) last
%
% Revision 1.142  2004/02/19 15:49:58  scott
% plot dipoles inside head in 'skirt' mode
%
% Revision 1.141  2004/02/18 01:16:53  scott
% help message adjust
%
% Revision 1.140  2004/02/18 01:02:58  scott
% 'dipole' help message. Adaptive AXHEADFAC.
%
% Revision 1.139  2004/02/17 22:44:54  arno
% now processing DIPFIT structure and fixed normalization bug
%
% Revision 1.138  2004/02/17 18:16:35  scott
% adjust EMARKERSIZE
%
% Revision 1.137  2004/02/17 18:11:36  scott
% fixed 'skirt'&'fill' problem. Also, made heads bigger
%
% Revision 1.136  2004/02/17 16:58:24  scott
% change color of outer 'shrink' mode ring to almost white, to avoid print bug
%
% Revision 1.135  2004/02/17 03:14:44  scott
% expand skirt border radius
%
% Revision 1.134  2004/02/15 21:30:01  scott
% same
%
% Revision 1.133  2004/02/15 21:17:07  scott
% omit QUAD_SKIRT option - not ready !
%
% Revision 1.132  2004/02/15 21:02:13  scott
% same
% Revision 1.96  2004/02/15 19:41:48  scott
% skirt with wedges
%
% Revision 1.95  2004/02/15 17:35:49  scott
% added 'style','skirt'
%
% Revision 1.72  2004/02/15 15:58:33  scott
% formatting, try 'shrink','skirt' ...
%
% Revision 1.71  2004/01/20 04:25:05  scott
% help msg edit
% .,
%
% Revision 1.70  2003/12/17 15:49:45  arno
% debug chan with no coordinates
%
% Revision 1.69  2003/12/17 01:25:37  arno
% debug plot electrode subset
%
% Revision 1.68  2003/12/17 00:57:17  arno
% subset of electrodes
%
% Revision 1.67  2003/11/29 23:34:00  scott
% help msg
%
% Revision 1.66  2003/11/06 16:31:18  arno
% changing dipnorm
%
% Revision 1.65  2003/11/06 02:04:41  arno
% correct orientation
%
% Revision 1.64  2003/11/06 01:40:31  arno
% diporient
%
% Revision 1.63  2003/11/06 01:00:57  arno
% adjusting corrdinates
% for dipole
%
% Revision 1.62  2003/11/05 20:35:21  arno
% dipole options
%
% Revision 1.61  2003/11/05 19:44:32  arno
% header text
%
% Revision 1.60  2003/08/08 17:36:12  arno
% shrink factor overwrite problem fixed
%
% Revision 1.59  2003/08/08 17:34:41  arno
% -cos -> cos
%
% Revision 1.58  2003/08/07 20:49:12  arno
% option 'masksurf' to speed up display
%
% Revision 1.57  2003/08/07 16:02:35  scott
% typo
%
% Revision 1.56  2003/08/07 16:01:49  scott
% debug
%
% Revision 1.55  2003/08/07 15:56:54  scott
% debug
%
% Revision 1.54  2003/08/07 15:54:49  scott
% debug last
%
% Revision 1.53  2003/08/07 15:51:05  scott
% added 'noplot' option to return interpolated channel value
%
% Revision 1.52  2003/07/18 01:34:07  scott
% text placement
%
% Revision 1.51  2003/07/18 01:33:19  scott
% text placement
%
% Revision 1.50  2003/07/18 01:31:49  scott
% debug
%
% Revision 1.49  2003/07/18 01:27:17  scott
% debug
%
% Revision 1.48  2003/07/18 01:26:05  scott
% debug
%
% Revision 1.47  2003/07/18 01:18:12  scott
% debug last
%
% Revision 1.46  2003/07/18 01:17:34  scott
% formatting, debug axes size message
%
% Revision 1.45  2003/07/17 23:42:32  scott
% nothing
%
% Revision 1.44  2003/07/17 23:13:03  scott
% rm debug message
%
% Revision 1.43  2003/07/16 16:29:46  arno
% replacing with topoplottest - added image output, gridscale arg
%
% Revision 1.41  2003/07/15 23:55:40  arno
% retreiving version 1.28
%
% Revision 1.28  2003/06/27 18:53:04  arno
% header msg
%
% Revision 1.27  2003/05/12 22:27:44  arno
% debug verbose
%
% Revision 1.26  2003/05/12 22:23:38  arno
% adding verbose option
%
% Revision 1.25  2002/11/27 01:23:53  arno
% change warning message
%
% Revision 1.24  2002/11/12 23:06:48  arno
% still debugging last insert
%
% Revision 1.23  2002/11/12 22:19:01  arno
% typo
%
% Revision 1.22  2002/11/12 21:43:51  scott
% tmpelocs -> tmpeloc
%
% Revision 1.21  2002/11/12 19:33:24  arno
% remove last channel of eloc structure if necessary (common ref)
%
% Revision 1.20  2002/11/01 03:50:08  erik
% same
%
% Revision 1.19  2002/11/01 03:47:40  erik
% added test for locs_file string to readlocs call
%
% Revision 1.18  2002/10/31 22:51:25  luca
% now also plotting n < nchans single channels
%
% Revision 1.17  2002/10/30 18:50:37  arno
% debugging dipole
%
% Revision 1.16  2002/10/30 16:41:21  arno
% adding the dipole option
%
% Revision 1.15  2002/10/26 20:09:35  arno
% error typo
%
% Revision 1.14  2002/10/14 00:40:44  arno
% *** empty log message ***
%
% Revision 1.13  2002/09/23 18:09:11  arno
% fixing single channel plotting
%
% Revision 1.12  2002/08/13 17:45:58  arno
% undo last change
%
% Revision 1.11  2002/08/13 17:44:37  arno
% remove color setting
%
% Revision 1.10  2002/08/12 01:34:53  arno
% color
%
% Revision 1.9  2002/08/11 22:31:20  arno
% color
%
% Revision 1.8  2002/05/01 18:49:20  arno
% modifying default shrink
%
% Revision 1.7  2002/05/01 02:40:10  arno
% typo
%
% Revision 1.6  2002/04/24 17:30:47  arno
% auto shrink
%
% Revision 1.5  2002/04/24 17:07:28  arno
% debugging error message problem
%
% Revision 1.4  2002/04/17 18:40:23  arno
% display real electrode number
%
% Revision 1.3  2002/04/06 03:47:44  arno
% adding emarkersize1chan input
%
% Revision 1.2  2002/04/06 03:37:24  arno
% adding single channel vector input
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% Topoplot Version 2.1
% Early development history:
% Begun by Andy Spydell and Scott Makeig, NHRC,  7-23-96
% 8-96 Revised by Colin Humphries, CNL / Salk Institute, La Jolla CA
%   -changed surf command to imagesc (faster)
%   -can now handle arbitrary scaling of electrode distances
%   -can now handle non integer angles in chan_locs
% 4-4-97 Revised again by Colin Humphries, reformatted by SM
%   -added parameters
%   -changed chan_locs format
% 2-26-98 Revised by Colin
%   -changed image back to surface command
%   -added fill and blank styles
%   -removed extra background colormap entry (now use any colormap)
%   -added parameters for electrode colors and labels
%   -now each topoplot axes use the caxis command again.
%   -removed OUTPUT parameter
% 3-11-98 changed default emarkersize, improve help msg -sm
% 5-24-01 made default emarkersize vary with number of channels -sm
% 01-25-02 reformated help & license, added link -ad 
% 03-15-02 added readlocs and the use of eloc input structure -ad 
% 03-25-02 added 'labelpoint' options and allow Values=[] -ad &sm
% 03-25-02 added details to "Unknown parameter" warning -sm & ad

function [handle,Zi,grid] = topoplot2(Values,loc_file,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10)

%
%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
icadefs                 % read defaults MAXTOPOPLOTCHANS and DEFAULT_ELOC and BACKCOLOR
if ~exist('BACKCOLOR')  % if icadefs.m does not define BACKCOLOR
   BACKCOLOR = [.93 .96 1];  % EEGLAB standard
end
plotgrid = 'off';
plotchans = [];
gridpos = 'l';          % default grid position left of head
noplot  = 'off';
handle = [];
Zi = [];
chanval = NaN;
rmax = 0.5;             % head radius - don't change this!
INTERPLIMITS = 'head';  % head, electrodes
INTSQUARE = 'on';       % default, interpolate electrodes located though the whole square containing
                        % the plotting disk
MAPLIMITS = 'absmax';   % absmax, maxmin, [values]
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
AXHEADFAC = 1.3;        % head to axes scaling factor
CONTOURNUM = 6;         % number of contour levels to plot
STYLE = 'both';         % default 'style': both,straight,fill,contour,blank
HEADCOLOR = [0 0 0];    % default head color (black)
CCOLOR = [0.2 0.2 0.2]; % default contour color
ECOLOR = [0 0 0];       % default electrode color
ELECTRODES = [];        % default 'electrodes': on|off|label - set below
MAXDEFAULTSHOWLOCS = 64;% if more channels than this, don't show electrode locations by default
EMARKER = '.';          % mark electrode locations with small disks
EMARKERSIZE = [];       % default depends on number of electrodes, set in code
EMARKERSIZE1CHAN = 40;  % default selected channel location marker size
EMARKERCOLOR1CHAN = 'red'; % selected channel location marker color
EFSIZE = get(0,'DefaultAxesFontSize'); % use current default fontsize for electrode labels
HLINEWIDTH = 3;         % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring 
HEADRINGWIDTH    = .007;% width of the cartoon head ring
SHADING = 'flat';       % default 'shading': flat|interp
shrinkfactor = [];      % shrink mode (dprecated)
intrad       = [];      % default interpolation square is to outermost electrode (<=1.0)
plotrad      = [];      % plotting radius ([] = auto, based on outermost channel location)
headrad      = [];      % default plotting radius for cartoon head is 0.5
MINPLOTRAD = 0.15;      % can't make a topoplot with smaller plotrad (contours fail)
VERBOSE = 'off';
MASKSURF = 'off';
CONVHULL = 'off';       % dont mask outside the electrodes convex hull

%%%%%% Dipole defaults %%%%%%%%%%%%
DIPOLE  = [];           
DIPNORM   = 'on';
DIPSPHERE = 85;
DIPLEN    = 1;
DIPSCALE  = 1;
DIPORIENT  = 1;
DIPCOLOR  = [0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%% Handle arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin< 1
   help topoplot;
   return
end
nargs = nargin;
if nargs == 1
  if isstr(Values)
    if any(strcmp(lower(Values),{'example','demo'}))
      fprintf(['This is an example of an electrode location file,\n',...
               'an ascii file consisting of the following four columns:\n',...
               ' channel_number degrees arc_length channel_name\n\n',...
               'Example:\n',...
               ' 1               -18    .352       Fp1 \n',...
               ' 2                18    .352       Fp2 \n',...
               ' 5               -90    .181       C3  \n',...
               ' 6                90    .181       C4  \n',...
               ' 7               -90    .500       A1  \n',...
               ' 8                90    .500       A2  \n',...
               ' 9              -142    .231       P3  \n',...
               '10               142    .231       P4  \n',...
               '11                 0    .181       Fz  \n',...
               '12                 0    0          Cz  \n',...
               '13               180    .181       Pz  \n\n',...
                                                             ...
               'In topoplot() coordinates, 0 deg. points to the nose, positive\n',...
               'angles point to the right hemisphere, and negative to the left.\n',...
               'The model head sphere has a circumference of 2; the vertex\n',...
               '(Cz) has arc_length 0. Locations with arc_length > 0.5 are below\n',...
               'head center and are plotted outside the head cartoon.\n'....
               'Option plotrad controls how much of this lower-head "skirt" is shown.\n',...
               'Option headrad controls if and where the cartoon head will be drawn.\n',...
               'Option intrad controls how many channels will be included in the interpolation.\n',...
               ])
      return
    end
  end
end
if nargs < 2
  loc_file = DEFAULT_ELOC;
  if ~exist(loc_file)
      fprintf('default locations file "%s" not found - specify chan_locs in topoplot() call.\n',loc_file)
      error(' ')
  end
end
if isempty(loc_file)
  loc_file = 0;
end
if isnumeric(loc_file) & loc_file == 0
  loc_file = DEFAULT_ELOC;
end

if nargs > 2
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
  end
  for i = 3:2:nargs
    Param = eval(['p',int2str((i-3)/2 +1)]);
    Value = eval(['v',int2str((i-3)/2 +1)]);
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch lower(Param)
         case 'conv'
           CONVHULL = lower(Value);
           if ~strcmp(CONVHULL,'on') & ~strcmp(CONVHULL,'off')
             error('Value of ''conv'' must be ''on'' or ''off''.');
          end
	 case 'colormap'
	  if size(Value,2)~=3
          error('Colormap must be a n x 3 matrix')
	  end
	  colormap(Value)
	 case 'intsquare'
          INTSQUARE = lower(Value);
          if ~strcmp(INTSQUARE,'on') & ~strcmp(INTSQUARE,'off')
             error('Value of ''intsquare'' must be ''on'' or ''off''.');
          end
	 case {'interplimits','headlimits'}
	  if ~isstr(Value)
          error('''interplimits'' value must be a string')
	  end
	  Value = lower(Value);
	  if ~strcmp(Value,'electrodes') & ~strcmp(Value,'head')
          error('Incorrect value for interplimits')
	  end
	  INTERPLIMITS = Value;
	 case 'verbose'
	  VERBOSE = Value;
	 case 'maplimits'
	  MAPLIMITS = Value;
	 case 'masksurf'
	  MASKSURF = Value;
	 case 'circgrid'
	  CIRCGRID = Value;
          if isstr(CIRCGRID) | CIRCGRID<100
            error('''circgrid'' value must be an int > 100');
          end
	 case 'style'
	  STYLE = lower(Value);
	 case 'numcontour'
	  CONTOURNUM = Value;
	 case 'electrodes'
	  ELECTRODES = lower(Value);
         if strcmpi(ELECTRODES,'pointlabels') | strcmpi(ELECTRODES,'ptslabels') ...
              | strcmpi(ELECTRODES,'labelspts') | strcmpi(ELECTRODES,'ptlabels') ...
              | strcmpi(ELECTRODES,'labelpts') 
             ELECTRODES = 'labelpoint'; % backwards compatability
         end
         if strcmpi(ELECTRODES,'pointnumbers') | strcmpi(ELECTRODES,'ptsnumbers') ...
              | strcmpi(ELECTRODES,'numberspts') | strcmpi(ELECTRODES,'ptnumbers') ...
              | strcmpi(ELECTRODES,'numberpts')  | strcmpi(ELECTRODES,'ptsnums')  ...
              | strcmpi(ELECTRODES,'numspts') 
             ELECTRODES = 'numpoint'; % backwards compatability
         end
         if strcmpi(ELECTRODES,'nums') 
             ELECTRODES = 'numbers'; % backwards compatability
         end
         if strcmpi(ELECTRODES,'pts') 
             ELECTRODES = 'on'; % backwards compatability
         end
         if ~strcmpi(ELECTRODES,'labelpoint') ...
            & ~strcmpi(ELECTRODES,'numpoint') ...
            & ~strcmp(ELECTRODES,'on') ...
            & ~strcmp(ELECTRODES,'off') ...
            & ~strcmp(ELECTRODES,'labels') ...
            & ~strcmpi(ELECTRODES,'numbers') 
              error('Unknown value for keyword ''electrodes''');
         end
	 case 'dipole'
	  DIPOLE = Value;
	 case 'dipsphere'
	  DIPSPHERE = Value;
	 case 'dipnorm'
	  DIPNORM = Value;
	 case 'diplen'
	  DIPLEN = Value;
	 case 'dipscale'
	  DIPSCALE = Value;
	 case 'diporient'
	  DIPORIENT = Value;
	 case 'dipcolor'
	  DIPCOLOR = Value;
	 case 'emarker'
	  EMARKER = Value;
	 case 'shrink'
	  shrinkfactor = Value;
	 case 'intrad'
	  intrad = Value;
          if isstr(intrad) | (intrad < MINPLOTRAD | intrad > 1)
	     error('intrad argument should be a number between 0.15 and 1.0');
	  end
	 case 'plotrad'
	  plotrad = Value;
          if isstr(plotrad) | (plotrad < MINPLOTRAD | plotrad > 1)
	     error('plotrad argument should be a number between 0.15 and 1.0');
	  end
	case 'headrad'
	  headrad = Value;
	  if isstr(headrad) & ( strcmpi(headrad,'off') | strcmpi(headrad,'none') )
	    headrad = 0;       % undocumented 'no head' alternatives
	  end
	  if isempty(headrad) % [] -> none also
	    headrad = 0;
	  end
	  if ~isstr(headrad) 
	    if ~(headrad==0) & (headrad < MINPLOTRAD | headrad>1)
	      error('bad value for headrad');
	    end
	  elseif  ~strcmpi(headrad,'rim')
	    error('bad value for headrad');
	  end
	 case {'headcolor','hcolor'}
	  HEADCOLOR = Value;
	 case {'contourcolor','ccolor'}
	  CCOLOR = Value;
	 case {'electcolor','ecolor'}
	  ECOLOR = Value;
	 case {'emarkersize','emsize'}
	  EMARKERSIZE = Value;
	 case 'emarkersize1chan'
	  EMARKERSIZE1CHAN= Value;
	 case {'efontsize','efsize'}
	  EFSIZE = Value;
	 case 'shading'
	  SHADING = lower(Value);
	  if ~any(strcmp(SHADING,{'flat','interp'}))
	     error('Invalid shading parameter')
	  end
         case 'noplot'
          noplot = Value;
          if ~isstr(noplot)
            if length(noplot) ~= 2
              error('''noplot'' location should be [radius, angle]')
            else
              chanrad = noplot(1);
              chantheta = noplot(2);
              noplot = 'on';
            end
          end
         case 'gridscale'
          GRID_SCALE = Value;
          if isstr(GRID_SCALE) | GRID_SCALE ~= round(GRID_SCALE) | GRID_SCALE < 32
               error('''gridscale'' value must be integer > 32.');
          end
         case {'plotgrid','gridplot'}
           plotgrid = 'on';
           if iscell(Value)
               if length(Value) > 2
                 error('''plotgrid'' value must be a matrix or cell array of length 2');
               end
               gridchans = Value{1};
               gridpos   = Value{2};
               if ~ismember(gridpos,'lro')
                 error('''plotgrid'' position must be ''l'', ''r'', or ''o''');
               end
           else
               gridchans = Value;
           end
         case 'plotchans'
           plotchans = Value(:);
           if find(plotchans==0) 
               error('''plotchans'' values must all be > 0');
           end
           if max(abs(plotchans))>max(Values) | max(abs(plotchans))>length(Values)
               error('''plotchans'' values must be <= max channel index');
           end
	 otherwise
	  error(['Unknown input parameter ''' Param ''' ???'])
    end
  end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% test args for plotting an electrode grid %%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(plotgrid,'on')
   if abs(max(max(gridchans))) > length(Values)
        error('''plotgrid'' channel index larger than the number of input channel values');
   end
   gchans = sort(find(abs(gridchans(:))>0));
   if setdiff(gchans,unique(gchans))
        fprintf('topoplot() warning: ''plotgrid'' channel matrix has duplicate channels\n');
   end
   if plotchans
     if intersect(gchans,abs(plotchans))
        fprintf('topoplot() warning: ''plotgrid'' and ''plotchans'' have channels in common\n');
     end
   end
   error('''plotgrid'' option not yet implemented.'); % <============
end

if isempty(ELECTRODES)                     % if electrode labeling not specified
  if length(Values) > MAXDEFAULTSHOWLOCS   % if more channels than default max
    ELECTRODES = 'off';                    % don't show electrodes
  else                                     % else if fewer chans,
    ELECTRODES = 'on';                     % do
  end
end

if isempty(Values)
   STYLE = 'blank';
end
[r,c] = size(Values);
if r>1 & c>1,
  error('input data must be a single vector');
end
Values = Values(:); % make Values a column vector

if ~isempty(intrad) & ~isempty(plotrad) & intrad < plotrad
   error('intrad must be >= plotrad');
end

%
%%%%%%%%%%%%%%%%%%%% Read the channel location information %%%%%%%%%%%%%%%%%%%%%%%%
% 
if isstr(loc_file)
	[tmpeloc labels Th Rd indices] = readlocs(loc_file,'filetype','loc');
else % a locs struct
	[tmpeloc labels Th Rd indices] = readlocs(loc_file);
        % Note: Th and Rd correspond to indices channels-with-coordinates only
end

if length(tmpeloc) == length(Values) + 1 % remove last channel if necessary 
                                         % (common reference channel)
    tmpeloc(end) = [];
    labels(end) = [];
    Th(end) = [];
    Rd(end) = [];
    ni = length(indices);
    fixit = find(indices==ni);
    if ~isempty(ni)
        indices(fixit) = [];
    end
end;
Th = pi/180*Th;                              % convert degrees to radians

%
%%%%%%%%%% if channels-to-mark-only are given in Values vector %%%%%%%%%%%%%%%%%
%
if length(Values) < length(tmpeloc) 
  if isempty(plotchans)
    if Values ~= round(Values) % if not integer values
      error('plotting fewer channels than in chanlocs: needs channel indices in ''plotchans''');
    elseif strcmpi(VERBOSE, 'on')
        fprintf('topoplot(): max chan number (%d) in locs > channels in data (%d).\n',...
                                   max(indices),length(Values));
        fprintf('            Marking the locations of the %d indicated channels.\n', ...
                                    length(Values));
    end
    plotchans = Values;
    STYLE = 'blank'; % plot channels only, marking the indicated channel number
    if strcmpi(ELECTRODES,'off')
         ELECTRODES = 'on';
    end
  elseif length(plotchans) ~= length(Values)
    error('number of channel values must = number of ''plotchans''');
  end
end

%
%%%%%%%%%%%%%%%%%%% remove infinite and NaN values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if length(Values) > 1
    inds = union(find(isnan(Values)), find(isinf(Values)));
    indices = setdiff(indices, inds);
    tmpeloc(inds) = [];
    Th(inds) = [];
    Rd(inds) = [];
end;

if length(Values) > 1
   if max(indices)>length(Values)
      STYLE = 'blank';
   else
      Values = Values(indices);
   end
end;
labels = labels(indices); % remove labels for electrodes without locations
labels = strvcat(labels); % make a label string matrix

if ~isempty(Values) & ~strcmpi( STYLE, 'blank')
  plotchans = 1:length(Th);
end
%
%%%%%%%%%%%%%%%%%% Read plotting radius from chanlocs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(plotrad) & isfield(tmpeloc, 'plotrad'), 
    plotrad = tmpeloc(1).plotrad; 
    if isstr(plotrad)                        % plotrad shouldn't be a string
        plotrad = str2num(plotrad)           % just checking
    end
    if plotrad < MINPLOTRAD | plotrad > 1.0
       fprintf('Bad value (%g) for plotrad.\n',plotrad);
       error(' ');
    end
    if strcmpi(VERBOSE,'on') & ~isempty(plotrad)
       fprintf('Plotting radius plotrad (%g) set from EEG.chanlocs.\n',plotrad);
    end
end;
if isempty(plotrad) 
  plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
  plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
end                                           % don't plot channels with Rd > 1 (below head)
if isempty(intrad) 
  intrad = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location
end                                           % don't interpolate channels with Rd > 1 (below head)

if isstr(plotrad) | plotrad < MINPLOTRAD | plotrad > 1.0
   error('plotrad must be between 0.15 and 1.0');
end

% plotrad now set

%
%%%%%%%%%%%%%%%%%%%%%%% Set radius of head cartoon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(headrad)  % never set -> defaults
  if plotrad >= rmax
     headrad = rmax;  % (anatomically correct)
  else % if plotrad < rmax
     headrad = 0;    % don't plot head
     if strcmpi(VERBOSE, 'on')
       fprintf('topoplot(): not plotting cartoon head since plotrad (%5.4g) < 0.5\n',...
                                                                    plotrad);
     end
  end
elseif strcmpi(headrad,'rim') % force plotting at rim of map
  headrad = plotrad;
end

% headrad now set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shrink mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(shrinkfactor) | isfield(tmpeloc, 'shrink'), 
    
    if isempty(shrinkfactor) & isfield(tmpeloc, 'shrink'), 
        shrinkfactor = tmpeloc(1).shrink;
        if strcmpi(VERBOSE,'on')
            if isstr(shrinkfactor)
                fprintf('Automatically shrinking coordinates\n');
            else                
                fprintf('Automatically shrinking coordinates by 3.2f\n', shrinkfactor);
            end;
        end
    end;
    
    if isstr(shrinkfactor)
        if strcmpi(shrinkfactor, 'on') | strcmpi(shrinkfactor, 'force') | strcmpi(shrinkfactor, 'auto')  
            if strcmpi(VERBOSE,'on')
                fprintf('     Shrink flag -> plotting cartoon head at plotrad\n');
            end
            headrad = plotrad; % plot head around outer electrodes, no matter if 0.5 or not
        end
    else % apply shrinkfactor
        plotrad = rmax/(1-shrinkfactor);
        headrad = plotrad;  % make deprecated 'shrink' mode plot 
        if strcmpi(VERBOSE,'on')
            fprintf('    %g%% shrink  applied.');
            if abs(headrad-rmax) > 1e-2
                fprintf(' Cartoon head is not anatomically correct.\n');
            else
                fprintf('\n');
            end
        end
    end
end; % if shrink
      
%
%%%%%%%%%%%%%%%%% Issue warning if headrad ~= rmax  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

if headrad ~= 0.5 & strcmpi(VERBOSE, 'on')
   fprintf('     NB: Plotting map using ''plotrad'' %-4.3g,',plotrad);
   fprintf(    ' ''headrad'' %-4.3g\n',headrad);
   fprintf('Warning: The plotting radius of the cartoon head is NOT anatomically correct (0.5).\n')
end
%
%%%%%%%%%%%%%%%%%%%%% Find plotting channels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle

[x,y] = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates

if strcmpi(INTSQUARE,'on') &  ~strcmpi(STYLE,'blank') % interpolate channels in the radius intrad square
  intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square
else
  intchans = find(Rd <= intrad); % interpolate channels in the radius intrad circle only
end

if length(pltchans) < length(Rd) & strcmpi(VERBOSE, 'on')
        fprintf('Interpolating %d and plotting %d of the %d scalp electrodes.\n', ...
                   length(intchans),length(pltchans),length(Rd));    
end;	
%
%%%%%%%%%%%%%%%%%%%%% Eliminate channels not plotted  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

if strcmp(plotgrid,'on')
     plotchans(gchans) = [];   % remove grid chans from head plotchans
end

allx = x;
ally = y;

Values(find(plotchans<0)) = -Values(find(plotchans<0)); % reverse indicated channel polarities

intchans = intersect(intchans,abs(plotchans)); % plot only indicated 'plotchans' channels
pltchans = intersect(pltchans,abs(plotchans)); % plot only indicated 'plotchans' channels
% fprintf('topoplot(): plotting %d channels\n',length(pltchans));

if ~isempty(Values)
	if length(Values) == length(Th)  % if as many map Values as channel locs
		intValues = Values(intchans);
		Values = Values(pltchans);
	else if strcmp(STYLE,'blank')    % else if Values holds numbers of channels to mark
            tmpValues=[];
            cc=1;
            for kk=1:length(Values)
                tmpind = find(pltchans == Values(kk));
                if ~isempty(tmpind)
                    tmpValues(cc) = tmpind;
                    cc=cc+1;
                end;
            end
            Values=tmpValues;     % eliminate the channel indices outside plotting area
		end;
	end;	
end;   % now channel parameters and values all refer to plotting channels only

intTh = Th(intchans);           % eliminate channels outside the interpolation area
intRd = Rd(intchans);
intx  = x(intchans);
inty  = y(intchans);

Th    = Th(pltchans);              % eliminate channels outside the plotting area
Rd    = Rd(pltchans);
x     = x(pltchans);
y     = y(pltchans);

labels= labels(pltchans,:);
%
%%%%%%%%%%%%%%% Squeeze channel locations to <= rmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

squeezefac = rmax/plotrad;
intRd = intRd*squeezefac; % squeeze electrode arc_lengths towards the vertex
Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
                          % to plot all inside the head cartoon
intx = intx*squeezefac;   
inty = inty*squeezefac;  
x    = x*squeezefac;    
y    = y*squeezefac;   
allx    = allx*squeezefac;    
ally    = ally*squeezefac;   

% Note: Now outermost channel will be plotted just inside rmax

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~strcmpi(STYLE,'blank') % if draw interpolated scalp map
  %
  %%%%%%%%%%%%%%%% Find limits for interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if strcmpi(INTERPLIMITS,'head')
    xmin = min(-rmax,min(intx)); xmax = max(rmax,max(intx));
    ymin = min(-rmax,min(inty)); ymax = max(rmax,max(inty));

  else % INTERPLIMITS = rectangle containing electrodes   ??? SAME ???
    xmin = max(-rmax,min(intx)); xmax = min(rmax,max(intx));
    ymin = max(-rmax,min(inty)); ymax = min(rmax,max(inty));
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%% Interpolate scalp map data %%%%%%%%%%%%%%%%%%%%%%%%
  %
  xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
  yi = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

  [Xi,Yi,Zi] = griddata(inty,intx,intValues,yi',xi,'invdist'); % interpolate data
  %
  %%%%%%%%%%%%%%%%%%%%%%% Mask out data outside the head %%%%%%%%%%%%%%%%%%%%%
  %
  mask = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
  ii = find(mask == 0);
  Zi(ii) = NaN;                         % mask non-plotting voxels with NaNs
  grid = [];

  %
  %%%%%%%%%% Return interpolated value at designated scalp location %%%%%%%%%%
  %
  if exist('chanrad')   % optional first argument to 'noplot' 
      chantheta = (chantheta/360)*2*pi;
      chancoords = round(ceil(GRID_SCALE/2)+GRID_SCALE/2*2*chanrad*[cos(-chantheta),...
                                                      -sin(-chantheta)]);
      if chancoords(1)<1 ...
         | chancoords(1) > GRID_SCALE ...
            | chancoords(2)<1 ...
               | chancoords(2)>GRID_SCALE
          error('designated ''noplot'' channel out of bounds')
      else
        chanval = Zi(chancoords(1),chancoords(2));
        grid = Zi;
        Zi = chanval;  % return interpolated value instead of Zi
      end
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%% Return interpolated image only  %%%%%%%%%%%%%%%%%
  %
   if strcmpi(noplot, 'on') 
    if strcmpi(VERBOSE,'on')
       fprintf('topoplot(): no plot requested.\n')
    end
    return;
   end
  %
  %%%%%%%%%%%%%%%%%%%%%%% Calculate colormap limits %%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  m = size(colormap,1);
  if isstr(MAPLIMITS)
    if strcmp(MAPLIMITS,'absmax')
      amin = -max(max(abs(Zi)));
      amax = max(max(abs(Zi)));
    elseif strcmp(MAPLIMITS,'maxmin') | strcmp(MAPLIMITS,'minmax')
      amin = min(min(Zi));
      amax = max(max(Zi));
    else
      error('unknown ''maplimits'' value.');
    end
  else
    amin = MAPLIMITS(1);
    amax = MAPLIMITS(2);
  end
  delta = xi(2)-xi(1); % length of grid entry

  %
  %%%%%%%%%%%%%%%%%%%%%%%%%% Scale the axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  cla  % clear current axis
  hold on
  h = gca; % uses current axes

                          % instead of default larger AXHEADFAC 
  if squeezefac<0.92 & plotrad-headrad > 0.05  % (size of head in axes)
    AXHEADFAC = 1.05;     % do not leave room for external ears if head cartoon
                          % shrunk enough by the 'skirt' option
  end

  set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC);
                          % specify size of head axes in gca

  unsh = (GRID_SCALE+1)/GRID_SCALE; % un-shrink the effects of 'interp' SHADING
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Plot map contours only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if strcmp(STYLE,'contour')                     % plot surface contours only
    [cls chs] = contour(Xi,Yi,Zi,CONTOURNUM,'k'); 
    % for h=chs, set(h,'color',CCOLOR); end
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Else plot map and contours %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'both')  % plot interpolated surface and surface contours
    if strcmp(SHADING,'interp')
       tmph = surface(Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,...
               'EdgeColor','none','FaceColor',SHADING);                    
    else % SHADING == 'flat'
       tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,...
               'EdgeColor','none','FaceColor',SHADING);                    
    end
    if strcmpi(MASKSURF, 'on')
        set(tmph, 'visible', 'off');
        handle = tmph;
    end;
    [cls chs] = contour(Xi,Yi,Zi,CONTOURNUM,'k'); 
    for h=chs, set(h,'color',CCOLOR); end
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Else plot map only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'straight') | strcmp(STYLE,'map') % 'straight' was former arg

      if strcmp(SHADING,'interp') % 'interp' mode is shifted somehow... but how?
         tmph = surface(Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,'EdgeColor','none',...
                  'FaceColor',SHADING);
      else
         tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none',...
                 'FaceColor',SHADING);
      end
    if strcmpi(MASKSURF, 'on')
        set(tmph, 'visible', 'off');
        handle = tmph;
    end;
  %
  %%%%%%%%%%%%%%%%%% Else fill contours with uniform colors  %%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'fill')
    [cls chs] = contourf(Xi,Yi,Zi,CONTOURNUM,'k');

    % for h=chs, set(h,'color',CCOLOR); end 
    %     <- 'not line objects.' Why does 'both' work above???

  else
    error('Invalid style')
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set color axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  caxis([amin amax]) % set coloraxis

else % if STYLE 'blank'
%
%%%%%%%%%%%%%%%%%%%%%%% Draw blank head %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if strcmpi(noplot, 'on') 
   if strcmpi(VERBOSE,'on')
      fprintf('topoplot(): no plot requested.\n')
   end
   return;
  end
  cla
  hold on

  set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC)
   % pos = get(gca,'position');
   % fprintf('Current axes size %g,%g\n',pos(3),pos(4));

  if strcmp(ELECTRODES,'labelpoint') |  strcmp(ELECTRODES,'numpoint')
    text(-0.6,-0.6, ...
    [ int2str(length(Rd)) ' of ' int2str(length(tmpeloc)) ' electrode locations shown']);
    text(-0.6,-0.7, [ 'Click on electrodes to toggle name/number']);
    tl = title('Channel locations');
    set(tl, 'fontweight', 'bold');
  end;
end

if exist('handle') ~= 1
    handle = gca;
end;

%
%%%%%%%%%%%%%%%%%%% Plot filled ring to mask jagged grid boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
hwidth = HEADRINGWIDTH;                   % width of head ring 
hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius

if strcmp(SHADING,'interp')
  rwidth = BLANKINGRINGWIDTH*1.3;             % width of blanking outer ring
else
  rwidth = BLANKINGRINGWIDTH;         % width of blanking outer ring
end
rin    =  rmax*(1-rwidth/2);              % inner ring radius
if hin>rin
  rin = hin;                              % dont blank inside the head ring
end

if strcmp(CONVHULL,'on') %%%%%%%%% mask outside the convex hull of the electrodes %%%%%%%%%
  cnv = convhull(allx,ally);
  cnvfac = round(CIRCGRID/length(cnv)); % spline interpolate the convex hull
  if cnvfac < 1, cnvfac=1; end;
  CIRCGRID = cnvfac*length(cnv);

  startangle = atan2(allx(cnv(1)),ally(cnv(1)));
  circ = linspace(0+startangle,2*pi+startangle,CIRCGRID);
  rx = sin(circ); 
  ry = cos(circ); 

  allx = allx(:)';  % make x (elec locations; + to nose) a row vector
  ally = ally(:)';  % make y (elec locations, + to r? ear) a row vector
  erad = sqrt(allx(cnv).^2+ally(cnv).^2);  % convert to polar coordinates
  eang = atan2(allx(cnv),ally(cnv));
  eang = unwrap(eang);
  eradi =spline(linspace(0,1,3*length(cnv)), [erad erad erad], ...
                                      linspace(0,1,3*length(cnv)*cnvfac));
  eangi =spline(linspace(0,1,3*length(cnv)), [eang+2*pi eang eang-2*pi], ...
                                      linspace(0,1,3*length(cnv)*cnvfac));
  xx = eradi.*sin(eangi);           % convert back to rect coordinates
  yy = eradi.*cos(eangi);
  yy = yy(CIRCGRID+1:2*CIRCGRID);
  xx = xx(CIRCGRID+1:2*CIRCGRID);
  eangi = eangi(CIRCGRID+1:2*CIRCGRID);
  eradi = eradi(CIRCGRID+1:2*CIRCGRID);
  xx = xx*1.02; yy = yy*1.02;           % extend spline outside electrode marks

  splrad = sqrt(xx.^2+yy.^2);           % arc radius of spline points (yy,xx)
  oob = find(splrad >= rin);            %  enforce an upper bound on xx,yy
  xx(oob) = rin*xx(oob)./splrad(oob);   % max radius = rin
  yy(oob) = rin*yy(oob)./splrad(oob);   % max radius = rin

  splrad = sqrt(xx.^2+yy.^2);           % arc radius of spline points (yy,xx)
  oob = find(splrad < hin);             % don't let splrad be inside the head cartoon
  xx(oob) = hin*xx(oob)./splrad(oob);   % min radius = hin
  yy(oob) = hin*yy(oob)./splrad(oob);   % min radius = hin

  ringy = [[ry(:)' ry(1) ]*(rin+rwidth) yy yy(1)];
  ringx = [[rx(:)' rx(1) ]*(rin+rwidth) xx xx(1)];

  ringh2= patch(ringy,ringx,ones(size(ringy)),BACKCOLOR,'edgecolor','none'); hold on

  % plot(ry*rmax,rx*rmax,'b') % debugging line

else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mask the jagged border around rmax %%%%%%%%%%%%%%%5%%%%%%

  circ = linspace(0,2*pi,CIRCGRID);
  rx = sin(circ); 
  ry = cos(circ); 
  ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
  ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];

  if ~strcmpi(STYLE,'blank')
    ringh= patch(ringx,ringy,0.01*ones(size(ringx)),BACKCOLOR,'edgecolor','none'); hold on
  end
  % plot(ry*rmax,rx*rmax,'b') % debugging line
end

  %f1= fill(rin*[rx rX],rin*[ry rY],BACKCOLOR,'edgecolor',BACKCOLOR); hold on
  %f2= fill(rin*[rx rX*(1+rwidth)],rin*[ry rY*(1+rwidth)],BACKCOLOR,'edgecolor',BACKCOLOR);

% Former line-style border smoothing - width did not scale with plot
%  brdr=plot(1.015*cos(circ).*rmax,1.015*sin(circ).*rmax,...      % old line-based method
%      'color',HEADCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH);    % plot skirt outline
%  set(brdr,'color',BACKCOLOR,'linewidth',HLINEWIDTH + 4);        % hide the disk edge jaggies 

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot cartoon head, ears, nose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if headrad > 0                         % if cartoon head to be plotted
%
%%%%%%%%%%%%%%%%%%% Plot head outline %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];

ringh= patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR); hold on

% rx = sin(circ); rX = rx(end:-1:1);
% ry = cos(circ); rY = ry(end:-1:1);
% for k=2:2:CIRCGRID
%   rx(k) = rx(k)*(1+hwidth);
%   ry(k) = ry(k)*(1+hwidth);
% end
% f3= fill(hin*[rx rX],hin*[ry rY],HEADCOLOR,'edgecolor',HEADCOLOR); hold on
% f4= fill(hin*[rx rX*(1+hwidth)],hin*[ry rY*(1+hwidth)],HEADCOLOR,'edgecolor',HEADCOLOR);

% Former line-style head
%  plot(cos(circ).*squeezefac*headrad,sin(circ).*squeezefac*headrad,...
%      'color',HEADCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH);    % plot head outline

%
%%%%%%%%%%%%%%%%%%% Plot ears and nose %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  base  = rmax-.0046;
  basex = 0.18*rmax;                   % nose width
  tip   = 1.15*rmax; 
  tiphw = .04*rmax;                    % nose tip half width
  tipr  = .01*rmax;                    % nose tip rounding
  q = .04; % ear lengthening
  EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
  EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
  sf    = headrad/plotrad;                                          % squeeze the model ears and nose 
                                                                    % by this factor
  plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
         2*ones(size([basex;tiphw;0;-tiphw;-basex])),...
         'Color',HEADCOLOR,'LineWidth',HLINEWIDTH);                 % plot nose
  plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)    % plot left ear
  plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)   % plot right ear
end

%
% %%%%%%%%%%%%%%%%%%% Show electrode information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
 plotax = gca;
 axis square                                           % make plotax square
 axis off

 pos = get(gca,'position');
 xlm = get(gca,'xlim');
 ylm = get(gca,'ylim');
 % textax = axes('position',pos,'xlim',xlm,'ylim',ylm);  % make new axes so clicking numbers <-> labels 
                                                       % will work inside head cartoon patch
 % axes(textax);                   
 axis square                                           % make textax square

 pos = get(gca,'position');
 set(plotax,'position',pos);

 xlm = get(gca,'xlim');
 set(plotax,'xlim',xlm);

 ylm = get(gca,'ylim');
 set(plotax,'ylim',ylm);                               % copy position and axis limits again

%get(textax,'pos')    % test if equal!
%get(plotax,'pos')
%get(textax,'xlim')
%get(plotax,'xlim')
%get(textax,'ylim')
%get(plotax,'ylim')

 if isempty(EMARKERSIZE)
   EMARKERSIZE = 10;
   if length(y)>=32 
    EMARKERSIZE = 8;
   elseif length(y)>=48
    EMARKERSIZE = 6;
   elseif length(y)>=64
    EMARKERSIZE = 5;
   elseif length(y)>=80
    EMARKERSIZE = 4;
   elseif length(y)>=100
    EMARKERSIZE = 3;
   elseif length(y)>=128
    EMARKERSIZE = 2;
   elseif length(y)>=160
    EMARKERSIZE = 1;
   end
 end
%
%%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations only %%%%%%%%%%%%%%%%%%%%%%%%%%
%
ELECTRODE_HEIGHT = 2.1;  % z value for plotting electrode information (above the surf)

if strcmp(ELECTRODES,'on')   % plot electrodes as spots
  hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE);
%
%%%%%%%%%%%%%%%%%%%%%%%% Print electrode labels only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
    for i = 1:size(labels,1)
    text(double(y(i)),double(x(i)),...
        ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',ECOLOR,...
	'FontSize',EFSIZE)
  end
%
%%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus labels %%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'labelpoint') 
  hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE);
  for i = 1:size(labels,1)
    hh(i) = text(double(y(i)+0.01),double(x(i)),...
        ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','left',...
	'VerticalAlignment','middle','Color', ECOLOR,'userdata', num2str(pltchans(i)), ...
	'FontSize',EFSIZE, 'buttondownfcn', ...
	    ['tmpstr = get(gco, ''userdata'');'...
	     'set(gco, ''userdata'', get(gco, ''string''));' ...
	     'set(gco, ''string'', tmpstr); clear tmpstr;'] );
  end
%
%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus numbers %%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'numpoint') 
  hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE);
  for i = 1:size(labels,1)
    hh(i) = text(double(y(i)+0.01),double(x(i)),...
        ELECTRODE_HEIGHT,num2str(pltchans(i)),'HorizontalAlignment','left',...
	'VerticalAlignment','middle','Color', ECOLOR,'userdata', labels(i,:) , ...
	'FontSize',EFSIZE, 'buttondownfcn', ...
	    ['tmpstr = get(gco, ''userdata'');'...
	     'set(gco, ''userdata'', get(gco, ''string''));' ...
	     'set(gco, ''string'', tmpstr); clear tmpstr;'] );
  end
%
%%%%%%%%%%%%%%%%%%%%%% Print electrode numbers only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'numbers')
  for i = 1:size(labels,1)
    text(double(y(i)),double(x(i)),...
        ELECTRODE_HEIGHT,int2str(pltchans(i)),'HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',ECOLOR,...
	'FontSize',EFSIZE)
  end
end
%
%%%%%%%% Mark specified electrode locations with red filled disks  %%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(STYLE,'blank') % if mark-selected-channel-locations mode
  if strcmpi(ELECTRODES,'on') | strcmpi(ELECTRODES,'off')
   for kk = 1:length(plotchans)
     if strcmpi(EMARKER,'.')
        hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', EMARKERCOLOR1CHAN, ...
                                              'markersize', EMARKERSIZE1CHAN);
     else
        hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', EMARKERCOLOR1CHAN, ...
                                              'markersize', EMARKERSIZE);
     end
   end
   hold on
  end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot dipole(s) on the scalp map  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(DIPOLE)  
    hold on;
    tmp = DIPOLE;
    if isstruct(DIPOLE)
        if ~isfield(tmp,'posxyz')
           error('dipole structure is not an EEG.dipfit.model')
        end
        DIPOLE = [];  % Note: invert x and y from dipplot usage
        DIPOLE(:,1) = -tmp.posxyz(:,2)/DIPSPHERE; % -y -> x
        DIPOLE(:,2) =  tmp.posxyz(:,1)/DIPSPHERE; %  x -> y
        DIPOLE(:,3) = -tmp.momxyz(:,2);
        DIPOLE(:,4) =  tmp.momxyz(:,1);
    else
        DIPOLE(:,1) = -tmp(:,2);                    % same for vector input
        DIPOLE(:,2) =  tmp(:,1);
        DIPOLE(:,3) = -tmp(:,4);
        DIPOLE(:,4) =  tmp(:,3);
    end;
    for index = 1:size(DIPOLE,1)
        if ~any(DIPOLE(index,:))
             DIPOLE(index,:) = [];
        end
    end;
    DIPOLE(:,1:4)   = DIPOLE(:,1:4)*rmax*(rmax/plotrad); % scale radius from 1 -> rmax (0.5)
    DIPOLE(:,3:end) = (DIPOLE(:,3:end))*rmax/100000*(rmax/plotrad); 
    if strcmpi(DIPNORM, 'on')
        for index = 1:size(DIPOLE,1)
            DIPOLE(index,3:4) = DIPOLE(index,3:4)/norm(DIPOLE(index,3:end))*0.2;
        end;
    end;
    DIPOLE(:, 3:4) =  DIPORIENT*DIPOLE(:, 3:4)*DIPLEN;

    PLOT_DIPOLE=1;
    if sum(DIPOLE(1,3:4).^2) <= 0.00001  
      if strcmpi(VERBOSE,'on')
        fprintf('Note: dipole is length 0 - not plotted\n')
      end
      PLOT_DIPOLE = 0;
    end
    if 0 % sum(DIPOLE(1,1:2).^2) > plotrad
      if strcmpi(VERBOSE,'on')
        fprintf('Note: dipole is outside plotting area - not plotted\n')
      end
      PLOT_DIPOLE = 0;
    end
    if PLOT_DIPOLE
      for index = 1:size(DIPOLE,1)
        hh = plot( DIPOLE(index, 1), DIPOLE(index, 2), '.');
        set(hh, 'color', DIPCOLOR, 'markersize', DIPSCALE*30);
        hh = line( [DIPOLE(index, 1) DIPOLE(index, 1)+DIPOLE(index, 3)]', ...
                   [DIPOLE(index, 2) DIPOLE(index, 2)+DIPOLE(index, 4)]');
        set(hh, 'color', DIPCOLOR, 'linewidth', DIPSCALE*30/7);
      end;
    end;
end;

%
%%%%%%%%%%%%% Set EEGLAB background color to match head border %%%%%%%%%%%%%%%%%%%%%%%%
%
try, 
  icadefs; 
  set(gcf, 'color', BACKCOLOR); 
  catch, 
end; 

hold off
axis off
return
