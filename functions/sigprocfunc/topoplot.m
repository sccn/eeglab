% topoplot() - plot a topographic map of a scalp data field in a 2-D circular view 
%              (looking down at the top of the head) using interpolation on a fine 
%              cartesian grid. Can also show specified channnel location(s), or return 
%              an interpolated value at an arbitrary scalp location (see 'noplot').
%              By default, channel locations below head center (arc_length 0.5) are 
%              shown in a 'skirt' outside the cartoon head (see 'plotrad' and 'headrad' 
%              options below). Nose is at top of plot; left is left; right is right.
% Usage:
%        >>  topoplot(datavector, EEG.chanlocs); % use a channel locations structure
%        >>  topoplot(datavector, 'my_chan.locs'); % use a channel locations file
%        >>  [h val grid] = topoplot(datavector, chan_locs, 'Param1','Value1', ...);
%
% Required Inputs:
%   datavector        - single vector of channel values. Else, if a vector of selected 
%                       channel numbers -> mark their location(s) using 'style' 'blank'.
%   chan_locs         - name of an EEG electrode position file (>> topoplot example),
%                       else an EEG.chanlocs structure (>> help pop_editset)
% Optional Inputs:
%   'maplimits'       - 'absmax' scale to +/- the absolute-max; 'maxmin', scale to 
%                       the data range; [clim1,clim2], user-definined lo/hi limits
%                       {default: 'absmax'}
%   'style'           - 'straight' - plot colored map only
%                       'contour'  - plot contour lines only
%                       'both'     - plot both colored map and contour lines
%                       'fill'     - plot constant color between contour lines
%                       'blank'    - plot electrode locations only {default: 'both'}
%   'electrodes'      - 'on','off','labels','numbers','ptslabels','ptsnumbers'
%                       {default: 'on', marks electrode location points}
%   'electcolor'{'k'}|'emarker'{'.'}|'emarkersize'{14}|'emarkersize1chan'{40}|'efontsize'{var}
%                       electrode marking details and their {defaults}. 
%   'numcontour'      - number of contour lines {default: 6}
%   'shading'         - 'flat','interp'  {default: 'flat'}
%   'interplimits'    - ['electrodes'|'head'] 'electrodes'-> interpolate the electrode grid; 
%                       'head'-> interpolate the whole disk {default: 'head'}.
%   'plotrad'         - [0.15<=float<=1.0] plotting radius. (See >> topoplot example). If 
%                       plotrad > 0.5, plot channels with arc_length > 0.5 (below ears-eyes)
%                       in a circular 'skirt' outside the cartoon head outline. 
%                       {default: the larger of the max channel arc_length and 0.5}. 
%   'headrad'         - plotting radius (arc_length) of cartoon head. Value in range(0.15,1.0)
%                       -> specify arc_length. NB: 0.5 is anatomically correct; 0 -> none; 
%                       'rim' -> show cartoon head at outer edge of the plot {default: 0.5}
%   'colormap'        -  (n,3) any size colormap {default: existing colormap}
%   'verbose'         - ['on'|'off'] comment on operations on command line {default: 'on'}.
%   'noplot'          - ['on'|'off'|[rad theta]] do not plot (but return interpolated data).
%                       If [rad theta] are coordinates of a (possibly missing) channel, 
%                       returns interpolated value for channel location. For location 
%                       conventions, see >> topoplot 'example' {default: 'off'}
%   'ccolor'          - color of the contours {default: blue}
%   'hcolor'|'ecolor' - colors of the cartoon head and electrodes {default: black}
%   'gridscale'       - [int >> 1] - interpolated data matrix size (rows) (default: 67)
%   'shrink'          - ['on'|'off'|'force'|factor] 'on' -> If max channel arc_length > 0.5, 
%                       shrink electrode coordinates towards vertex to plot all channels
%                       by making max arc_length 0.5. 'force' -> Normalize arc_length 
%                       so the channel max is 0.5. factor -> Apply a specified shrink
%                       factor (range (0,1) = shrink fraction). {default: 'off'}
% Dipole plotting options:
%   'dipole'          - [xi yi xe ye ze] plot dipole on the top of the scalp map
%                       from coordinate (xi,yi) to coordinates (xe,ye,ze) (dipole head 
%                       model has radius 1). If several rows, plot one dipole per row.
%                       Coordinates returned by dipplot() may be used. Can accept
%                       an EEG.dipfit.model structure (See >> help dipplot).
%                       Ex: 'dipole',EEG.dipfit.model(17) % Plot dipole(s) for comp. 17.
%   'dipnorm'         - ['on'|'off'] normalize dipole length {default: 'off'}.
%   'diporient'       - [-1|1] invert dipole orientation {default: 1}.
%   'diplen'          - [real] scale dipole length {default: 1}.
%   'dipscale'        - [real] scale dipole size {default: 1}.
%   'dipsphere'       - [real] size of the dipole sphere. {default: 85 mm}.
%   'dipcolor'        - [color] dipole color as Matlab code code or [r g b] vector
%                       {default: 'k' -> black}.
% Outputs:
%         h           - plot axes handle
%         val         - interpolated value at given 'noplot' channel location, if any.
%         grid        - interpolated data image matrix (off-head points = NaN).
%
% Eloc_file format:
%    chan_number degrees arc_length reject_level amp_gain channel_name
%    (Angle-0 =Cz-to-Fz; C3-angle =-90; arc_length at edge of image = 0.5)
%    For more information and sample file: >> topoplot 'example'
%
% Authors: Andy Spydell, Colin Humphries, Scott Makeig & Arnaud Delorme 
%          CNL / Salk Institute, Aug, 1996 - SCCN/INC/UCSD, Nov. 2001 -
%
% See also: timtopo(), envtopo()

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

function [handle,chanval,Zi] = topoplot2(Values,loc_file,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10)

%
%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
icadefs                 % read defaults MAXTOPOPLOTCHANS and DEFAULT_ELOC and BACKCOLOR
if ~exist('BACKCOLOR')  % if icadefs.m does not define BACKCOLOR
   BACKCOLOR = [.93 .96 1];  % EEGLAB standard
end
noplot  = 'off';
handle = [];
Zi = [];
chanval = NaN;
rmax = 0.5;             % head radius - don't change this!
INTERPLIMITS = 'head';  % head, electrodes
MAPLIMITS = 'absmax';   % absmax, maxmin, [values]
GRID_SCALE = 67;        % plot map on a 67X67 grid
AXHEADFAC = 1.3;        % head to axes scaling factor
CONTOURNUM = 6;         % number of contour levels to plot
STYLE = 'both';         % default 'style': both,straight,fill,contour,blank
HCOLOR = [0 0 0];       % default head color
CCOLOR = [0 0 1];       % default contour color
ECOLOR = [0 0 0];       % default electrode color
ELECTRODES = 'on';      % default 'electrodes': on|off|label
EMARKER = '.';          % mark electrode locations with small disks
EMARKERSIZE = [];       % default depends on number of electrodes, set in code
EMARKERSIZE1CHAN = 40;  % default selected channel location marker size
EMARKERCOLOR1CHAN = 'red'; % selected channel location marker color
EFSIZE = get(0,'DefaultAxesFontSize'); % use current default fontsize for electrode labels
HLINEWIDTH = 2;         % default linewidth for head, nose, ears
SHADING = 'flat';       % default 'shading': flat|interp
shrinkfactor = 'off';   % shrinking mode
plotrad      = [];      % plotting radius ([] = auto, based on outermost channel location)
headrad      = [];      % default plotting radius for cartoon head is 0.5
MINPLOTRAD = 0.15;      % can't make a topoplot with smaller plotrad (contours fail)
VERBOSE = 'on';
MASKSURF = 'off';

%%%%%% Dipole defaults %%%%%%%%%%%%
DIPOLE  = [];           
DIPNORM   = 'off';
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
if nargs < 2
  loc_file = DEFAULT_ELOC;
  if ~exist(loc_file)
      fprintf('default locations file "%s" not found - specify chan_locs in topoplot() call.\n',loc_file)
      error(' ')
  end
end
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
               ])
      return
    end
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
    error('topoplot(): Odd number of inputs?')
  end
  for i = 3:2:nargs
    Param = eval(['p',int2str((i-3)/2 +1)]);
    Value = eval(['v',int2str((i-3)/2 +1)]);
    if ~isstr(Param)
      error('topoplot(): Parameter must be a string')
    end
    Param = lower(Param);
    switch lower(Param)
	 case 'colormap'
	  if size(Value,2)~=3
          error('topoplot(): Colormap must be a n x 3 matrix')
	  end
	  colormap(Value)
	 case {'interplimits','headlimits'}
	  if ~isstr(Value)
          error('topoplot(): interplimits value must be a string')
	  end
	  Value = lower(Value);
	  if ~strcmp(Value,'electrodes') & ~strcmp(Value,'head')
          error('topoplot(): Incorrect value for interplimits')
	  end
	  INTERPLIMITS = Value;
	 case 'verbose'
	  VERBOSE = Value;
	 case 'maplimits'
	  MAPLIMITS = Value;
	 case 'masksurf'
	  MASKSURF = Value;
	 case 'gridscale'
	  GRID_SCALE = Value;
	 case 'style'
	  STYLE = lower(Value);
	 case 'numcontour'
	  CONTOURNUM = Value;
	 case 'electrodes'
	  ELECTRODES = lower(Value);
         if strcmpi(ELECTRODES,'pointlabels') | strcmpi(ELECTRODES,'ptslabels')
             ELECTRODES = 'labelpoint'; % backwards compatability
         end
         if strcmpi(ELECTRODES,'pointnumbers') | strcmpi(ELECTRODES,'ptsnumbers')
             ELECTRODES = 'numpoint'; % backwards compatability
         end
         if ~strcmpi(ELECTRODES,'labelpoint') ...
            & ~strcmpi(ELECTRODES,'numpoint') ...
            & ~strcmp(ELECTRODES,'on') ...
            & ~strcmp(ELECTRODES,'off') ...
            & ~strcmp(ELECTRODES,'labels') ...
            & ~strcmpi(ELECTRODES,'numbers') 
                fprintf('topoplot(): Unknown value for keyword ''electrodes''.\n');
            return
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
	  SHRINKSET = 1;        % flag 'shrink' set
	  shrinkfactor = Value;
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
	  HCOLOR = Value;
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
          if GRID_SCALE ~= round(GRID_SCALE) | GRID_SCALE < 4
               fprintf('topoplot(): ''gridscale'' value must be integer > 4.\n');
               return
          end
	 otherwise
	  error(['topoplot(): Unknown input parameter ''' Param ''' ???'])
    end
  end
end

if isempty(Values)
   STYLE = 'blank';
end
[r,c] = size(Values);
if r>1 & c>1,
  error('topoplot(): input data must be a single vector');
elseif r==1 & c==1
  STYLE = 'blank'; % plot channels only, marking the indicated channel number
end

%
%%%%%%%%%%%%%%%%%%%% Read the channel location information %%%%%%%%%%%%%%%%%%%%%%%%
% 
if isstr(loc_file)
	[tmpeloc labels Th Rd indices] = readlocs(loc_file,'filetype','loc');
else % a locs struct
	[tmpeloc labels Th Rd indices] = readlocs(loc_file);
end
if length(tmpeloc) == length(Values) + 1 % remove last channel if necessary 
                                         % (common reference channel)
    tmpeloc(end) = [];
    labels(end) = [];
    Th(end) = [];
    Rd(end) = [];
end;
Th = pi/180*Th;                              % convert degrees to radians

if length(Values) > 1
   if max(indices)>length(Values)
      if strcmpi(VERBOSE, 'on')
        fprintf('topoplot(): max chan number (%d) in locs > channels in data (%d).\n',...
                                   max(indices),length(Values));
        fprintf('            Marking the locations of the %d indicated channels.\n', ...
                                    length(Values));
      end
      STYLE = 'blank';
   else
      Values     = Values(indices);
   end
end;
labels = labels(indices);
labels = strvcat(labels);

%
%%%%%%%%%%%%%%%%%% Read plotting radius from chanlocs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(plotrad) & isfield(tmpeloc, 'plotrad'), 
    plotrad = tmpeloc(1).plotrad; 
    if isstr(plotrad)                        % plotrad shouldn't be a string
        plotrad = str2num(plotrad)           % just checking
    end
    if plotrad < MINPLOTRAD | plotrad > 1.0
       fprintf('Bad value (%g) for plotrad.\n');
       error(' ');
    end
    if strcmpi(VERBOSE,'on') & ~isempty(plotrad)
       fprintf('Plotting radius plotrad (%g) set from EEG.chanlocs.\n',plotrad);
    end
end;
if isempty(plotrad) 
  plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
end                                           % don't plot channels with Rd > 1 (below head)

if isstr(plotrad) | plotrad < MINPLOTRAD | plotrad > 1.0
   error('plotrad must be between 0.15 and 1.0');
end

% plotrad now set

%
%%%%%%%%%%%%%%%%%%%%%%% Set radius of head cartoon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(headrad)  % never set -> defaults
  if plotrad >= 0.5
     headrad = 0.5;  % (anatomically correct)
  else % if plotrad < 0.5
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

if headrad ~= 0.5 & strcmpi(VERBOSE, 'on')
   fprintf('NB: Plotting map using: plotrad %-4.3g,',plotrad);
   fprintf(    ' headrad %-4.3g\n',headrad);
   fprintf('    The cartoon head is NOT anatomically correct (headrad 0.5).\n')
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shrink mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
if exist('SHRINKSET') % plotrad and headrad set above from those args and/or chanlocs
                      % and/or defaults - do not change them without careful thought!

% NOTE: This section should ONLY adjust plotrad and headrad based on 'shrink' values

%if strcmpi(shrinkfactor, 'off') & isfield(tmpeloc, 'shrink'), 
%   shrinkfactor = tmpeloc(1).shrink;         % read default shrinkfactor from chanlocs
%end;
%if ~isstr(shrinkfactor)
%    shrinkfactor = 0.5+0.5*shrinkfactor;     % convert input shrinkfactor value ????
%end;                                         
%
%if isstr(shrinkfactor) & ~strcmpi(shrinkfactor, 'off') & ~isempty(plotrad)
%    shrinkfactor = plotrad;                  % not specified 'shrink' mode
%end;
%
%if isstr(shrinkfactor)                       % if shrink arg is 'on' or 'force'
%  if (strcmp(lower(shrinkfactor), 'on') & max(Rd) >rmax) ...
%                   | strcmp(lower(shrinkfactor),'force') 
%    shrinkfactor = rmax/max(Rd);   % was ??? (2*max(r)-1)/(2*rmax);
%    if shrinkfactor > 1            % shrink, do not expand
%            shrinkfactor = 1;
%    elseif strcmpi(VERBOSE, 'on')
%            fprintf(...
%         'topoplot(): electrode arc_lengths shrunk by (1 -> %2.3g) to plot all\n', ...
%                shrinkfactor);
%    end;
%  end;
%
%else                                         % if numeric shrinkfactor read or set above
%    if strcmpi(VERBOSE, 'on')
%        fprintf('topoplot(): electrode arc_lengths shrunk by (1 -> %2.3g)\n',...
%                                                                  shrinkfactor);
%    end;
%end;
end % SHRINKSET

%
%%%%%%%%%%%%%%%%%%%%% Find plotting channels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
enum = find(Rd <= plotrad); % only interpolate and plot channels inside potting radius

if length(enum) < length(Rd) & strcmpi(VERBOSE, 'on')
        fprintf('Not using or showing the outermost %d of the %d electrodes.\n', ...
                   length(Rd)-length(enum),length(Rd));    
end;	
%
%%%%%%%%%%%%%%%%%%%%% Eliminate channels not plotted  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

if ~isempty(Values)
	if length(Values) == length(Th)  % if as many map Values as channel locs
		Values = Values(enum);
	else if strcmp(STYLE,'blank')    % else if Values holds numbers of channels to mark
            tmpValues=[];
            cc=1;
            for kk=1:length(Values)
                tmpind = find(enum == Values(kk));
                if ~isempty(tmpind)
                    tmpValues(cc) = tmpind;
                    cc=cc+1;
                end;
            end
            Values=tmpValues;     % eliminate the channel indices outside plotting area
		end;
	end;	
end;   % now channel parameters and values all refer to plotting channels only

Th = Th(enum);              % eliminate channels outside the ploting area
Rd = Rd(enum);
labels = labels(enum,:);

%
%%%%%%%%%%%%%%% Squeeze channel locations to <= rmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

squeezefac = rmax/plotrad;          
Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
                          % to plot all inside the head cartoon
% Note: Now outermost channel will be plotted just inside rmax

[x,y] = pol2cart(Th,Rd);  % transform from polar to cartesian coordinates

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~strcmpi(STYLE,'blank') % if draw interpolated scalp map
  %
  %%%%%%%%%%%%%%%% Find limits for interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if strcmp(INTERPLIMITS,'head')
    xmin = min(-rmax,min(x)); xmax = max(rmax,max(x));
    ymin = min(-rmax,min(y)); ymax = max(rmax,max(y));
  else % interplimits = rectangle containing electrodes
    xmin = max(-rmax,min(x)); xmax = min(rmax,max(x));
    ymin = max(-rmax,min(y)); ymax = min(rmax,max(y));
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%% Interpolate scalp map data %%%%%%%%%%%%%%%%%%%%%%%%
  %
  xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
  yi = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

  [Xi,Yi,Zi] = griddata(y,x,Values,yi',xi,'invdist'); % interpolate data
  %
  %%%%%%%%%%%%%%%%%%%%%%% Mask out data outside the head %%%%%%%%%%%%%%%%%%%%%
  %
  mask = (sqrt(Xi.^2+Yi.^2) <= rmax);
  ii = find(mask == 0);
  Zi(ii) = NaN;

  %
  %%%%%%%%%% Return interpolated value at designated scalp location %%%%%%%%%%
  %
  if exist('chanrad') % indicates optional 'noplot' (first argument)
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
      end
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%% Return interpolated image only  %%%%%%%%%%%%%%%%%
  %
   if strcmpi(noplot, 'on') 
       fprintf('topoplot(): no plot requested.\n')
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
      fprintf('topoplot(): unknown ''maplimits'' value.\n');
      return
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

  set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC)

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
    tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,...
               'EdgeColor','none','FaceColor',SHADING);                    
    if strcmpi(MASKSURF, 'on')
        set(tmph, 'visible', 'off');
        handle = tmph;
    end;
    [cls chs] = contour(Xi,Yi,Zi,CONTOURNUM,'k'); 
    for h=chs, set(h,'color',CCOLOR); end
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Else plot map only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'straight')
    tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none',...
	'FaceColor',SHADING);
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
    error('topoplot(): Invalid style')
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
      fprintf('topoplot(): no plot requested.\n')
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
%%%%%%%%%%%%%%%%%%%% Plot cartoon head, ears, nose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if headrad > 0                         % if cartoon head to be plotted
  circ  = 0:2*pi/100:2*pi;             % 101 circle vertices
  basex = 0.18*rmax;                   % nose width
  tip   = rmax*1.15; 
  base  = rmax-.004;
  EarX  = [.497  .510  .518  .5299 .5419  .54    .547   .532   .510   .489]; % rmax = 0.5
  EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];

  if headrad == plotrad                       % plot head outline at plot rim
     brdr=plot(1.015*cos(circ).*rmax,1.015*sin(circ).*rmax,...
      'color',HCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH);   % plot skirt outline
     set(brdr,'color',BACKCOLOR,'linewidth',HLINEWIDTH + 4);      % hide the disk edge jaggies 
  end

  plot(cos(circ).*squeezefac*headrad,sin(circ).*squeezefac*headrad,...
      'color',HCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH);    % plot head outline

  sf    = headrad/plotrad;                                       % squeeze the model ears and nose 
                                                                 % by this factor
  plot([basex;0;-basex]*sf,[base;tip;base]*sf,...
         'Color',HCOLOR,'LineWidth',HLINEWIDTH);                % plot nose
  plot(EarX*sf,EarY*sf,'color',HCOLOR,'LineWidth',HLINEWIDTH)   % plot left ear
  plot(-EarX*sf,EarY*sf,'color',HCOLOR,'LineWidth',HLINEWIDTH)  % plot right ear
end

%
% %%%%%%%%%%%%%%%%%%% Show electrode information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
if strcmp(ELECTRODES,'on')   % plot electrodes as spots
  hp2 = plot(y,x,EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE);
%
%%%%%%%%%%%%%%%%%%%%%%%% Print electrode labels only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
    for i = 1:size(labels,1)
    text(y(i),x(i),labels(i,:),'HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',ECOLOR,...
	'FontSize',EFSIZE)
  end
%
%%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus labels %%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'labelpoint') 
  hp2 = plot(y,x,EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE);
  for i = 1:size(labels,1)
    hh(i) = text(y(i)+0.01,x(i),labels(i,:),'HorizontalAlignment','left',...
	'VerticalAlignment','middle','Color', ECOLOR,'userdata', num2str(enum(i)), ...
	'FontSize',EFSIZE, 'buttondownfcn', ...
	    ['tmpstr = get(gco, ''userdata'');'...
	     'set(gco, ''userdata'', get(gco, ''string''));' ...
	     'set(gco, ''string'', tmpstr); clear tmpstr;'] );
  end
%
%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus numbers %%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'numpoint') 
  hp2 = plot(y,x,EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE);
  for i = 1:size(labels,1)
    hh(i) = text(y(i)+0.01,x(i),num2str(enum(i)),'HorizontalAlignment','left',...
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
    text(y(i),x(i),int2str(enum(i)),'HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',ECOLOR,...
	'FontSize',EFSIZE)
  end
end
%
%%%%%%%%%%%%%%%%%%%%%% Mark specified electrode locations %%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(STYLE,'blank') % if mark-selected-channel-locations mode
  if length(Values) < length(enum)  % mark selected electrodes
      for kk = 1:length(Values)
        hp2 = plot(y(Values(kk)),x(Values(kk)),'.','Color', EMARKERCOLOR1CHAN, ...
                                              'markersize', EMARKERSIZE1CHAN);
        hold on
      end
  end;
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
    DIPOLE(:,1:4)   = DIPOLE(:,1:4)*rmax^2/plotrad; % scale radius from 1 -> rmax 
    DIPOLE(:,3:end) = (DIPOLE(:,3:end)/10000)*rmax^2/plotrad; % NB: 10000 was 500 !?!?
    if strcmpi(DIPNORM, 'on')
        for index = 1:size(DIPOLE,1)
            DIPOLE(index,3:4) = DIPOLE(index,3:4)/norm(DIPOLE(index,3:end))*0.2;
        end;
    end;
    DIPOLE(:, 3:4) =  DIPORIENT*DIPOLE(:, 3:4)*DIPLEN;

    PLOT_DIPOLE=1;
    if sum(DIPOLE(1,3:4).^2) <= 0.0001  
      if strcmpi(VERBOSE,'on')
        fprintf('Note: dipole is length 0 - not plotted\n')
      end
      PLOT_DIPOLE = 0;
    end
    if sum(DIPOLE(1,1:2).^2) > plotrad
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

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Keep head round  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axis square; % keep head round!

hold off
axis off

return
