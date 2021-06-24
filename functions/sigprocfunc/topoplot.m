% topoplot() - plot a topographic map of a scalp data field in a 2-D circular view 
%              (looking down at the top of the head) using interpolation on a fine 
%              cartesian grid. Can also show specified channnel location(s), or return 
%              an interpolated value at an arbitrary scalp location (see 'noplot').
%              By default, channel locations below head center (arc_length 0.5) are 
%              shown in a 'skirt' outside the cartoon head (see 'plotrad' and 'headrad' 
%              options below). Nose is at top of plot; left is left; right is right.
%              Using option 'plotgrid', the plot may be one or more rectangular grids.
% Usage:
%        >>  topoplot(datavector, EEG.chanlocs);   % plot a map using an EEG chanlocs structure
%        >>  topoplot(datavector, 'my_chan.locs'); % read a channel locations file and plot a map
%        >>  topoplot('example');                  % give an example of an electrode location file
%        >>  [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
%                           topoplot(datavector, chan_locs, 'Input1','Value1', ...);
% Required Inputs:
%   datavector        - single vector of channel values. Else, if a vector of selected subset
%                       (int) channel numbers -> mark their location(s) using 'style' 'blank'.
%   chan_locs         - name of an EEG electrode position file (>> topoplot example).
%                       Else, an EEG.chanlocs structure (>> help readlocs or >> topoplot example)
% Optional inputs:
%   'maplimits'       - 'absmax'   -> scale map colors to +/- the absolute-max (makes green 0); 
%                       'maxmin'   -> scale colors to the data range (makes green mid-range); 
%                       [lo.hi]    -> use user-definined lo/hi limits
%                       {default: 'absmax'}
%   'style'           - 'map'      -> plot colored map only
%                       'contour'  -> plot contour lines only
%                       'both'     -> plot both colored map and contour lines
%                       'fill'     -> plot constant color between contour lines
%                       'blank'    -> plot electrode locations only {default: 'both'}
%   'electrodes'      - 'on','off','labels','numbers','ptslabels','ptsnumbers'. To set the 'pts' 
%                       marker,,see 'Plot detail options' below. {default: 'on' -> mark electrode 
%                       locations with points ('.') unless more than 64 channels, then 'off'}. 
%   'plotchans'       - [vector] channel numbers (indices) to use in making the head plot. 
%                       {default: [] -> plot all chans}
%   'plotgrid'        - [channels] Plot channel data in one or more rectangular grids, as 
%                       specified by [channels],  a position matrix of channel numbers defining 
%                       the topographic locations of the channels in the
%                       grid. Zero values are ignored (given the figure background color); 
%                       negative integers, the color of the polarity-reversed channel values.  
%                       Ex: >> figure; ...
%                             >> topoplot(values,'chanlocs','plotgrid',[11 12 0; 13 14 15]);
%                       % Plot a (2,3) grid of data values from channels 11-15 with one empty 
%                       grid cell (top right) {default: no grid plot} 
%   'nosedir'         - ['+X'|'-X'|'+Y'|'-Y'] direction of nose {default: '+X'}
%   'chaninfo'        - [struct] optional structure containing fields 'nosedir', 'plotrad'. 
%                       See these (separate) field definitions above, below.
%                       {default: nosedir +X, plotrad 0.5, all channels}
%   'plotrad'         - [0.15<=float<=1.0] plotting radius = max channel arc_length to plot.
%                       See >> topoplot example. If plotrad > 0.5, chans with arc_length > 0.5 
%                       (i.e. below ears-eyes) are plotted in a circular 'skirt' outside the
%                       cartoon head. See 'intrad' below. {default: max(max(chanlocs.radius),0.5);
%                       If the chanlocs structure includes a field chanlocs.plotrad, its value 
%                       is used by default}.
%   'headrad'         - [0.15<=float<=1.0] drawing radius (arc_length) for the cartoon head. 
%                       NOTE: Only headrad = 0.5 is anatomically correct! 0 -> don't draw head; 
%                       'rim' -> show cartoon head at outer edge of the plot {default: 0.5}
%   'intrad'          - [0.15<=float<=1.0] radius of the scalp map interpolation area (square or 
%                       disk, see 'intsquare' below). Interpolate electrodes in this area and use 
%                       this limit to define boundaries of the scalp map interpolated data matrix
%                       {default: max channel location radius}
%   'intsquare'       - ['on'|'off'] 'on' -> Interpolate values at electrodes located in the whole 
%                       square containing the (radius intrad) interpolation disk; 'off' -> Interpolate
%                       values from electrodes shown in the interpolation disk only {default: 'on'}.
%   'conv'            - ['on'|'off'] Show map interpolation only out to the convext hull of
%                       the electrode locations to minimize extrapolation. Use this option ['on'] when 
%                       plotting pvalues  {default: 'off'}. When plotting pvalues in totoplot, set 
%                       'conv' option to 'on' to minimize interpolation effects
%   'noplot'          - ['on'|'off'|[rad theta]] do not plot (but return interpolated data).
%                       Else, if [rad theta] are coordinates of a (possibly missing) channel, 
%                       returns interpolated value for channel location.  For more info, 
%                       see >> topoplot 'example' {default: 'off'}
%   'verbose'         - ['on'|'off'] comment on operations on command line {default: 'on'}.
%   'chantype'        - deprecated
%
% Plot detail options:
%   'drawaxis'        - ['on'|'off'] draw axis on the top left corner.
%   'emarker'         - Matlab marker char | {markerchar color size linewidth} char, else cell array 
%                       specifying the electrode 'pts' marker. Ex: {'s','r',32,1} -> 32-point solid 
%                       red square. {default: {'.','k',[],1} where marker size ([]) depends on the number 
%                       of channels plotted}.
%   'emarker2'        - {markchans}|{markchans marker color size linewidth} cell array specifying 
%                       an alternate marker for specified 'plotchans'. Ex: {[3 17],'s','g'} 
%                       {default: none, or if {markchans} only are specified, then {markchans,'o','r',10,1}}
%   'hcolor'          - color of the cartoon head. Use 'hcolor','none' to plot no head. {default: 'k' = black}
%   'shading'         - 'flat','interp'  {default: 'flat'}
%   'numcontour'      - number of contour lines {default: 6}. You may also enter a vector to set contours 
%                       at specified values.
%   'contourvals'     - values for contour {default: same as input values}
%   'pmask'           - values for masking topoplot. Array of zeros and 1 of the same size as the input 
%                       value array {default: []}
%   'color'           - color of the contours {default: dark grey}
%   'whitebk '        -  ('on'|'off') make the background color white (e.g., to print empty plotgrid channels) 
%                       {default: 'off'}
%   'gridscale'       - [int > 32] size (nrows) of interpolated scalp map data matrix {default: 67}
%   'colormap'        -  (n,3) any size colormap {default: existing colormap}
%   'circgrid'        - [int > 100] number of elements (angles) in head and border circles {201}
%   'emarkercolor'    - cell array of colors for 'blank' option.
%   'plotdisk'        - ['on'|'off'] plot disk instead of dots for electrodefor 'blank' option. Size of disk
%                       is controled by input values at each electrode. If an imaginary value is provided, 
%                       plot partial circle with red for the real value and blue for the imaginary one.
%
% Dipole plotting options:
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
% Outputs:
%              handle - handle of the colored surface.If
%                       contour only is plotted, then is the handle of
%                       the countourgroup. (If no surface or contour is plotted,
%                       return "gca", the handle of the current plot)
%         grid_or_val - [matrix] the interpolated data image (with off-head points = NaN).  
%                       Else, single interpolated value at the specified 'noplot' arg channel 
%                       location ([rad theta]), if any.
%     plotrad_or_grid - IF grid image returned above, then the 'plotrad' radius of the grid.
%                       Else, the grid image
%     xmesh, ymesh    - x and y values of the returned grid (above)
%
% Chan_locs format:
%    See >> topoplot 'example'
%
% Examples:
%
%    To plot channel locations only:
%    >> figure; topoplot([],EEG.chanlocs,'style','blank','electrodes','labelpoint','chaninfo',EEG.chaninfo);
%    
% Notes: - To change the plot map masking ring to a new figure background color,
%            >> set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
%        - Topoplots may be rotated. From the commandline >> view([deg 90]) {default: [0 90])
%        - When plotting pvalues make sure to use the option 'conv' to minimize extrapolation effects 
%
% Authors: Andy Spydell, Colin Humphries, Arnaud Delorme & Scott Makeig
%          CNL / Salk Institute, 8/1996-/10/2001; SCCN/INC/UCSD, Nov. 2001 -
%
% See also: timtopo(), envtopo()

% Deprecated options: 
%           'shrink' - ['on'|'off'|'force'|factor] Deprecated. 'on' -> If max channel arc_length 
%                       > 0.5, shrink electrode coordinates towards vertex to plot all channels
%                       by making max arc_length 0.5. 'force' -> Normalize arc_length 
%                       so the channel max is 0.5. factor -> Apply a specified shrink
%                       factor (range (0,1) = shrink fraction). {default: 'off'}
%   'electcolor' {'k'}  ... electrode marking details and their {defaults}. 
%   'emarker' {'.'}|'emarkersize' {14}|'emarkersizemark' {40}|'efontsize' {var} -
%                       electrode marking details and their {defaults}. 
%   'ecolor'          - color of the electrode markers {default: 'k' = black}
%   'interplimits'    - ['electrodes'|'head'] 'electrodes'-> interpolate the electrode grid; 
%                       'head'-> interpolate the whole disk {default: 'head'}.

% Unimplemented future options:

% Copyright (C) Colin Humphries & Scott Makeig, CNL / Salk Institute, Aug, 1996
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

function [handle,Zi,grid,Xi,Yi] = topoplot(Values,loc_file,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
icadefs                 % read defaults MAXTOPOPLOTCHANS and DEFAULT_ELOC and BACKCOLOR
if ~exist('BACKCOLOR')  % if icadefs.m does not define BACKCOLOR
   BACKCOLOR = [.93 .96 1];  % EEGLAB standard
end
whitebk = 'off';  % by default, make gridplot background color = EEGLAB screen background color

persistent warningInterp;

plotgrid = 'off';
plotchans = [];
noplot  = 'off';
handle = [];
Zi = [];
chanval = NaN;
rmax = 0.5;             % actual head radius - Don't change this!
INTERPLIMITS = 'head';  % head, electrodes
INTSQUARE = 'on';       % default, interpolate electrodes located though the whole square containing
                        % the plotting disk
default_intrad = 1;     % indicator for (no) specified intrad
MAPLIMITS = 'absmax';   % absmax, maxmin, [values]
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
AXHEADFAC = 1.3;        % head to axes scaling factor
CONTOURNUM = 6;         % number of contour levels to plot
STYLE = 'both';         % default 'style': both,straight,fill,contour,blank
HEADCOLOR = [0 0 0];    % default head color (black)
CCOLOR = [0.2 0.2 0.2]; % default contour color
ELECTRODES = [];        % default 'electrodes': on|off|label - set below
MAXDEFAULTSHOWLOCS = 64;% if more channels than this, don't show electrode locations by default
EMARKER = '.';          % mark electrode locations with small disks
ECOLOR = [0 0 0];       % default electrode color = black
EMARKERSIZE = [];       % default depends on number of electrodes, set in code
EMARKERLINEWIDTH = 1;   % default edge linewidth for emarkers
EMARKERSIZE1CHAN = 20;  % default selected channel location marker size
EMARKERCOLOR1CHAN = 'red'; % selected channel location marker color
EMARKER2CHANS = [];      % mark subset of electrode locations with small disks
EMARKER2 = 'o';          % mark subset of electrode locations with small disks
EMARKER2COLOR = 'r';     % mark subset of electrode locations with small disks
EMARKERSIZE2 = 10;      % default selected channel location marker size
EMARKER2LINEWIDTH = 1;
EFSIZE = get(0,'DefaultAxesFontSize'); % use current default fontsize for electrode labels
HLINEWIDTH = 2;         % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring 
HEADRINGWIDTH    = .007;% width of the cartoon head ring
SHADING = 'flat';       % default 'shading': flat|interp
shrinkfactor = [];      % shrink mode (dprecated)
intrad       = [];      % default interpolation square is to outermost electrode (<=1.0)
plotrad      = [];      % plotting radius ([] = auto, based on outermost channel location)
headrad      = [];      % default plotting radius for cartoon head is 0.5
squeezefac = 1.0;
MINPLOTRAD = 0.15;      % can't make a topoplot with smaller plotrad (contours fail)
VERBOSE = 'off';
MASKSURF = 'off';
CONVHULL = 'off';       % dont mask outside the electrodes convex hull
DRAWAXIS = 'off';
PLOTDISK = 'off';
ContourVals = Values;
PMASKFLAG   = 0;
COLORARRAY  = { [1 0 0] [0.5 0 0] [0 0 0] };
%COLORARRAY2 = { [1 0 0] [0.5 0 0] [0 0 0] };
gb = [0 0];
COLORARRAY2 = { [gb 0] [gb 1/4] [gb 2/4] [gb 3/4] [gb 1] };

%%%%%% Dipole defaults %%%%%%%%%%%%
DIPOLE  = [];           
DIPNORM   = 'on';
DIPNORMMAX = 'off';
DIPSPHERE = 85;
DIPLEN    = 1;
DIPSCALE  = 1;
DIPORIENT  = 1;
DIPCOLOR  = [0 0 0];
NOSEDIR   = '+X';
CHANINFO  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%% Handle arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin< 1
   help topoplot;
   return
end

% calling topoplot from Fieldtrip
% -------------------------------
fieldtrip = 0;
if nargin < 2, loc_file = []; end
if isstruct(Values) || ~isstruct(loc_file), fieldtrip == 1; end
if ischar(loc_file), if exist(loc_file) ~= 2, fieldtrip == 1; end; end
if fieldtrip
    error('Wrong calling format, are you trying to use the topoplot Fieldtrip function?');
end

nargs = nargin;
if nargs == 1
  if ischar(Values)
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
               'head center and are plotted outside the head cartoon.\n',...
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
if isnumeric(loc_file) && loc_file == 0
  loc_file = DEFAULT_ELOC;
end

if nargs > 2
    if ~(round(nargs/2) == nargs/2)
        error('Odd number of input arguments??')
    end
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
            error('Flag arguments must be strings')
        end
        Param = lower(Param);
        switch Param
            case 'conv'
                CONVHULL = lower(Value);
                if ~strcmp(CONVHULL,'on') && ~strcmp(CONVHULL,'off')
                    error('Value of ''conv'' must be ''on'' or ''off''.');
                end
            case 'colormap'
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
                colormap(Value)
            case 'gridscale'
                GRID_SCALE = Value;
            case 'plotdisk'
                PLOTDISK = lower(Value);
                if ~strcmp(PLOTDISK,'on') && ~strcmp(PLOTDISK,'off')
                    error('Value of ''plotdisk'' must be ''on'' or ''off''.');
                end
            case 'intsquare'
                INTSQUARE = lower(Value);
                if ~strcmp(INTSQUARE,'on') && ~strcmp(INTSQUARE,'off')
                    error('Value of ''intsquare'' must be ''on'' or ''off''.');
                end
            case 'emarkercolors'
                COLORARRAY = Value;
            case {'interplimits','headlimits'}
                if ~ischar(Value)
                    error('''interplimits'' value must be a string')
                end
                Value = lower(Value);
                if ~strcmp(Value,'electrodes') && ~strcmp(Value,'head')
                    error('Incorrect value for interplimits')
                end
                INTERPLIMITS = Value;
            case 'verbose'
                VERBOSE = Value;
            case 'nosedir'
                NOSEDIR = Value;
                if isempty(strmatch(lower(NOSEDIR), { '+x', '-x', '+y', '-y' }))
                    error('Invalid nose direction');
                end
            case 'chaninfo'
                CHANINFO = Value;
                if isfield(CHANINFO, 'nosedir'), NOSEDIR      = CHANINFO.nosedir; end
                if isfield(CHANINFO, 'shrink' ), shrinkfactor = CHANINFO.shrink;  end
                if isfield(CHANINFO, 'plotrad') && isempty(plotrad), plotrad = CHANINFO.plotrad; end
            case 'chantype'
            case 'drawaxis'
                DRAWAXIS = Value;
            case 'maplimits'
                MAPLIMITS = Value;
            case 'masksurf'
                MASKSURF = Value;
            case 'circgrid'
                CIRCGRID = Value;
                if ischar(CIRCGRID) || CIRCGRID<100
                    error('''circgrid'' value must be an int > 100');
                end
            case 'style'
                STYLE = lower(Value);
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
                if strcmpi(ELECTRODES,'pointlabels') || strcmpi(ELECTRODES,'ptslabels') ...
                        | strcmpi(ELECTRODES,'labelspts') | strcmpi(ELECTRODES,'ptlabels') ...
                        | strcmpi(ELECTRODES,'labelpts')
                    ELECTRODES = 'labelpoint'; % backwards compatability
                elseif strcmpi(ELECTRODES,'pointnumbers') || strcmpi(ELECTRODES,'ptsnumbers') ...
                        | strcmpi(ELECTRODES,'numberspts') | strcmpi(ELECTRODES,'ptnumbers') ...
                        | strcmpi(ELECTRODES,'numberpts')  | strcmpi(ELECTRODES,'ptsnums')  ...
                        | strcmpi(ELECTRODES,'numspts')
                    ELECTRODES = 'numpoint'; % backwards compatability
                elseif strcmpi(ELECTRODES,'nums')
                    ELECTRODES = 'numbers'; % backwards compatability
                elseif strcmpi(ELECTRODES,'pts')
                    ELECTRODES = 'on'; % backwards compatability
                elseif ~strcmp(ELECTRODES,'off') ...
                        & ~strcmpi(ELECTRODES,'on') ...
                        & ~strcmp(ELECTRODES,'labels') ...
                        & ~strcmpi(ELECTRODES,'numbers') ...
                        & ~strcmpi(ELECTRODES,'labelpoint') ...
                        & ~strcmpi(ELECTRODES,'numpoint')
                    error('Unknown value for keyword ''electrodes''');
                end
            case 'dipole'
                DIPOLE = Value;
            case 'dipsphere'
                DIPSPHERE = Value;
            case {'dipnorm', 'dipnormmax'}
                if strcmp(Param,'dipnorm')
                    DIPNORM = Value;
                    if strcmpi(Value,'on')
                        DIPNORMMAX = 'off';
                    end
                else
                    DIPNORMMAX = Value;
                    if strcmpi(Value,'on')
                        DIPNORM = 'off';
                    end
                end
                
            case 'diplen'
                DIPLEN = Value;
            case 'dipscale'
                DIPSCALE = Value;
            case 'contourvals'
                ContourVals = Value;
            case 'pmask'
                ContourVals = Value;
                PMASKFLAG   = 1;
            case 'diporient'
                DIPORIENT = Value;
            case 'dipcolor'
                DIPCOLOR = Value;
            case 'emarker'
                if ischar(Value)
                    EMARKER = Value;
                elseif ~iscell(Value) || length(Value) > 4
                    error('''emarker'' argument must be a cell array {marker color size linewidth}')
                else
                    EMARKER = Value{1};
                end
                if length(Value) > 1
                    ECOLOR = Value{2};
                end
                if length(Value) > 2
                    EMARKERSIZE = Value{3};
                end
                if length(Value) > 3
                    EMARKERLINEWIDTH = Value{4};
                end
            case 'emarker2'
                if ~iscell(Value) || length(Value) > 5
                    error('''emarker2'' argument must be a cell array {chans marker color size linewidth}')
                end
                EMARKER2CHANS = abs(Value{1}); % ignore channels < 0
                if length(Value) > 1
                    EMARKER2 = Value{2};
                end
                if length(Value) > 2
                    EMARKER2COLOR = Value{3};
                end
                if length(Value) > 3
                    EMARKERSIZE2 = Value{4};
                end
                if length(Value) > 4
                    EMARKER2LINEWIDTH = Value{5};
                end
            case 'shrink'
                shrinkfactor = Value;
            case 'intrad'
                intrad = Value;
                if ischar(intrad) || (intrad < MINPLOTRAD || intrad > 1)
                    error('intrad argument should be a number between 0.15 and 1.0');
                end
            case 'plotrad'
                plotrad = Value;
                if ~isempty(plotrad) && (ischar(plotrad) || (plotrad < MINPLOTRAD || plotrad > 1))
                    error('plotrad argument should be a number between 0.15 and 1.0');
                end
            case 'headrad'
                headrad = Value;
                if ischar(headrad) && ( strcmpi(headrad,'off') || strcmpi(headrad,'none') )
                    headrad = 0;       % undocumented 'no head' alternatives
                end
                if isempty(headrad) % [] -> none also
                    headrad = 0;
                end
                if ~ischar(headrad)
                    if ~(headrad==0) && (headrad < MINPLOTRAD || headrad>1)
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
            case {'emarkersize1chan','emarkersizemark'}
                EMARKERSIZE1CHAN= Value;
            case {'efontsize','efsize'}
                EFSIZE = Value;
            case 'shading'
                SHADING = lower(Value);
                if ~any(strcmp(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
                if strcmpi(SHADING,'interp') && isempty(warningInterp)
                    warning('Using interpolated shading in scalp topographies prevent to export them as vectorized figures');
                    warningInterp = 1;
                end
            case 'noplot'
                noplot = Value;
                if ~ischar(noplot)
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
                if ischar(GRID_SCALE) || GRID_SCALE ~= round(GRID_SCALE) || GRID_SCALE < 32
                    error('''gridscale'' value must be integer > 32.');
                end
            case {'plotgrid','gridplot'}
                plotgrid = 'on';
                gridchans = Value;
            case 'plotchans'
                plotchans = Value(:);
                if find(plotchans<=0)
                    error('''plotchans'' values must be > 0');
                end
                % if max(abs(plotchans))>max(Values) | max(abs(plotchans))>length(Values) -sm ???
            case {'whitebk','whiteback','forprint'}
                whitebk = Value;
            case {'iclabel'} % list of options to ignore
            otherwise
                error(['Unknown input parameter ''' Param ''' ???'])
        end
    end
end

if strcmpi(whitebk, 'on')
    BACKCOLOR = [ 1 1 1 ];
end

if isempty(find(strcmp(varargin,'colormap')))
    if exist('DEFAULT_COLORMAP','var')
        cmap = colormap(DEFAULT_COLORMAP);
    else
        cmap = parula;
    end
else
    cmap = colormap;
end
if strcmp(noplot,'on'), close(gcf); end
cmaplen = size(cmap,1);

if strcmp(STYLE,'blank')    % else if Values holds numbers of channels to mark
    if length(Values) < length(loc_file)
        ContourVals = zeros(1,length(loc_file));
        ContourVals(Values) = 1;
        Values = ContourVals;
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% test args for plotting an electrode grid %%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(plotgrid,'on')
   STYLE = 'grid';
   gchans = sort(find(abs(gridchans(:))>0));

   % if setdiff(gchans,unique(gchans))
   %      fprintf('topoplot() warning: ''plotgrid'' channel matrix has duplicate channels\n');
   % end

   if ~isempty(plotchans)
     if intersect(gchans,abs(plotchans))
        fprintf('topoplot() warning: ''plotgrid'' and ''plotchans'' have channels in common\n');
     end
   end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% misc arg tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
if r>1 && c>1,
  error('input data must be a single vector');
end
Values = Values(:); % make Values a column vector
ContourVals = ContourVals(:); % values for contour

if ~isempty(intrad) && ~isempty(plotrad) && intrad < plotrad
   error('intrad must be >= plotrad');
end

if ~strcmpi(STYLE,'grid')                     % if not plot grid only

%
%%%%%%%%%%%%%%%%%%%% Read the channel location information %%%%%%%%%%%%%%%%%%%%%%%%
% 
  if ischar(loc_file)
      [tmpeloc labels Th Rd indices] = readlocs( loc_file);
  elseif isstruct(loc_file) % a locs struct
      [tmpeloc labels Th Rd indices] = readlocs( loc_file );
      % Note: Th and Rd correspond to indices channels-with-coordinates only
  else
       error('loc_file must be a EEG.locs struct or locs filename');
  end
  Th = pi/180*Th;                              % convert degrees to radians
  allchansind = 1:length(Th);

  
  if ~isempty(plotchans)
      if max(plotchans) > length(Th)
          error('''plotchans'' values must be <= max channel index');
      end
  end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% channels to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(plotchans)
    plotchans = intersect_bc(plotchans, indices);
end
if ~isempty(Values) && ~strcmpi( STYLE, 'blank') && isempty(plotchans)
    plotchans = indices;
end
if isempty(plotchans) && strcmpi( STYLE, 'blank')
    plotchans = indices;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% filter channels used for components %%%%%%%%%%%%%%%%%%%%% 
%
if isfield(CHANINFO, 'icachansind') && ~isempty(Values) && length(Values) ~= length(tmpeloc)

    % test if ICA component
    % ---------------------
    if length(CHANINFO.icachansind) == length(Values)
        
        % if only a subset of channels are to be plotted
        % and ICA components also use a subject of channel
        % we must find the new indices for these channels
        
        plotchans = intersect_bc(CHANINFO.icachansind, plotchans);
        tmpvals   = zeros(1, length(tmpeloc));
        tmpvals(CHANINFO.icachansind) = Values;
        Values    = tmpvals;
        tmpvals   = zeros(1, length(tmpeloc));
        tmpvals(CHANINFO.icachansind) = ContourVals;
        ContourVals = tmpvals;
        
    end
end

%
%%%%%%%%%%%%%%%%%%% last channel is reference? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if length(tmpeloc) == length(Values) + 1 % remove last channel if necessary 
                                         % (common reference channel)
    if plotchans(end) == length(tmpeloc)
        plotchans(end) = [];
    end

end

%
%%%%%%%%%%%%%%%%%%% remove infinite and NaN values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if length(Values) > 1
    inds          = union_bc(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
    plotchans     = setdiff_bc(plotchans, inds);
end
if strcmp(plotgrid,'on')
    plotchans = setxor(plotchans,gchans);   % remove grid chans from head plotchans   
end

[x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
Th        = Th(plotchans);
Rd        = Rd(plotchans);
x         = x(plotchans);
y         = y(plotchans);
labels    = labels(plotchans); % remove labels for electrodes without locations
labels    = strvcat(labels); % make a label string matrix
if ~isempty(Values) && length(Values) > 1
    Values      = Values(plotchans);
    ContourVals = ContourVals(plotchans);
end

%
%%%%%%%%%%%%%%%%%% Read plotting radius from chanlocs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(plotrad) && isfield(tmpeloc, 'plotrad'), 
    plotrad = tmpeloc(1).plotrad; 
    if ischar(plotrad)                        % plotrad shouldn't be a string
        plotrad = str2num(plotrad)           % just checking
    end
    if plotrad < MINPLOTRAD || plotrad > 1.0
       fprintf('Bad value (%g) for plotrad.\n',plotrad);
       error(' ');
    end
    if strcmpi(VERBOSE,'on') && ~isempty(plotrad)
       fprintf('Plotting radius plotrad (%g) set from EEG.chanlocs.\n',plotrad);
    end
end
if isempty(plotrad) 
  plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
  plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
end                                           % don't plot channels with Rd > 1 (below head)

if isempty(intrad) 
  default_intrad = 1;     % indicator for (no) specified intrad
  intrad = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location
else
  default_intrad = 0;                         % indicator for (no) specified intrad
  if plotrad > intrad
     plotrad = intrad;
  end
end                                           % don't interpolate channels with Rd > 1 (below head)
if ischar(plotrad) || plotrad < MINPLOTRAD || plotrad > 1.0
   error('plotrad must be between 0.15 and 1.0');
end

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

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shrink mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(shrinkfactor) || isfield(tmpeloc, 'shrink'), 
    if isempty(shrinkfactor) && isfield(tmpeloc, 'shrink'), 
        shrinkfactor = tmpeloc(1).shrink;
        if strcmpi(VERBOSE,'on')
            if ischar(shrinkfactor)
                fprintf('Automatically shrinking coordinates to lie above the head perimter.\n');
            else                
                fprintf('Automatically shrinking coordinates by %3.2f\n', shrinkfactor);
            end
        end
    end
    
    if ischar(shrinkfactor)
        if strcmpi(shrinkfactor, 'on') || strcmpi(shrinkfactor, 'force') || strcmpi(shrinkfactor, 'auto')  
            if abs(headrad-rmax) > 1e-2
             fprintf('     NOTE -> the head cartoon will NOT accurately indicate the actual electrode locations\n');
            end
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
                fprintf(' Warning: With this "shrink" setting, the cartoon head will NOT be anatomically correct.\n');
            else
                fprintf('\n');
            end
        end
    end
end; % if shrink
      
%
%%%%%%%%%%%%%%%%% Issue warning if headrad ~= rmax  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

if headrad ~= 0.5 && strcmpi(VERBOSE, 'on')
   fprintf('     NB: Plotting map using ''plotrad'' %-4.3g,',plotrad);
   fprintf(    ' ''headrad'' %-4.3g\n',headrad);
   fprintf('Warning: The plotting radius of the cartoon head is NOT anatomically correct (0.5).\n')
end
%
%%%%%%%%%%%%%%%%%%%%% Find plotting channels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

pltchans = find(Rd <= plotrad); % plot channels inside plotting circle

if strcmpi(INTSQUARE,'on') % interpolate channels in the radius intrad square
  intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square
else
  intchans = find(Rd <= intrad); % interpolate channels in the radius intrad circle only
end

%
%%%%%%%%%%%%%%%%%%%%% Eliminate channels not plotted  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

allx      = x;
ally      = y;
intchans; % interpolate using only the 'intchans' channels
pltchans; % plot using only indicated 'plotchans' channels

if length(pltchans) < length(Rd) && strcmpi(VERBOSE, 'on')
        fprintf('Interpolating %d and plotting %d of the %d scalp electrodes.\n', ...
                   length(intchans),length(pltchans),length(Rd));    
end;	


% fprintf('topoplot(): plotting %d channels\n',length(pltchans));
if ~isempty(EMARKER2CHANS)
    if strcmpi(STYLE,'blank')
       error('emarker2 not defined for style ''blank'' - use marking channel numbers in place of data');
    else % mark1chans and mark2chans are subsets of pltchans for markers 1 and 2
       [tmp1, mark1chans, tmp2] = setxor(pltchans,EMARKER2CHANS);
       [tmp3, tmp4, mark2chans] = intersect_bc(EMARKER2CHANS,pltchans);
    end
end

if ~isempty(Values)
	if length(Values) == length(Th)  % if as many map Values as channel locs
		intValues      = Values(intchans);
		intContourVals = ContourVals(intchans);
        Values         = Values(pltchans);
		ContourVals    = ContourVals(pltchans);
	end;	
end;   % now channel parameters and values all refer to plotting channels only

allchansind = allchansind(pltchans);
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

else % if strcmpi(STYLE,'grid')
   intx = rmax; inty=rmax;
end % if ~strcmpi(STYLE,'grid')

%
%%%%%%%%%%%%%%%% rotate channels based on chaninfo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(lower(NOSEDIR), '+x')
     rotate = 0;
else
    if strcmpi(lower(NOSEDIR), '+y')
        rotate = 3*pi/2;
    elseif strcmpi(lower(NOSEDIR), '-x')
        rotate = pi;
    else rotate = pi/2;
    end
    allcoords = (inty + intx*sqrt(-1))*exp(sqrt(-1)*rotate);
    intx = imag(allcoords);
    inty = real(allcoords);
    allcoords = (ally + allx*sqrt(-1))*exp(sqrt(-1)*rotate);
    allx = imag(allcoords);
    ally = real(allcoords);
    allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*rotate);
    x = imag(allcoords);
    y = real(allcoords);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~strcmpi(STYLE,'blank') % if draw interpolated scalp map
 if ~strcmpi(STYLE,'grid') %  not a rectangular channel grid
  %
  %%%%%%%%%%%%%%%% Find limits for interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if default_intrad % if no specified intrad
   if strcmpi(INTERPLIMITS,'head') % intrad is 'head'
    xmin = min(-rmax,min(intx)); xmax = max(rmax,max(intx));
    ymin = min(-rmax,min(inty)); ymax = max(rmax,max(inty));

   else % INTERPLIMITS = rectangle containing electrodes -- DEPRECATED OPTION!
    xmin = max(-rmax,min(intx)); xmax = min(rmax,max(intx));
    ymin = max(-rmax,min(inty)); ymax = min(rmax,max(inty));
   end
  else % some other intrad specified
    xmin = -intrad*squeezefac; xmax = intrad*squeezefac;   % use the specified intrad value 
    ymin = -intrad*squeezefac; ymax = intrad*squeezefac;
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%% Interpolate scalp map data %%%%%%%%%%%%%%%%%%%%%%%%
  %
  xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
  yi = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

  try
      [Xi,Yi,Zi] = griddata(inty,intx,double(intValues),yi',xi,'v4'); % interpolate data
      [Xi,Yi,ZiC] = griddata(inty,intx,double(intContourVals),yi',xi,'v4'); % interpolate data
  catch,
      [Xi,Yi] = meshgrid(yi',xi);
      Zi  = gdatav4(inty,intx,double(intValues), Xi, Yi);
      ZiC = gdatav4(inty,intx,double(intContourVals), Xi, Yi);
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%% Mask out data outside the head %%%%%%%%%%%%%%%%%%%%%
  %
  mask = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
  ii = find(mask == 0);
  Zi(ii)  = NaN;                         % mask non-plotting voxels with NaNs  
  ZiC(ii) = NaN;                         % mask non-plotting voxels with NaNs
  grid = plotrad;                       % unless 'noplot', then 3rd output arg is plotrad
  %
  %%%%%%%%%% Return interpolated value at designated scalp location %%%%%%%%%%
  %
  if exist('chanrad')   % optional first argument to 'noplot' 
      chantheta = (chantheta/360)*2*pi;
      chancoords = round(ceil(GRID_SCALE/2)+GRID_SCALE/2*2*chanrad*[cos(-chantheta),...
                                                      -sin(-chantheta)]);
      if chancoords(1)<1 ...
         || chancoords(1) > GRID_SCALE ...
            || chancoords(2)<1 ...
               || chancoords(2)>GRID_SCALE
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
  if ischar(MAPLIMITS)
    if strcmp(MAPLIMITS,'absmax')
      amax = max(max(abs(Zi)));
      amin = -amax;
    elseif strcmp(MAPLIMITS,'maxmin') || strcmp(MAPLIMITS,'minmax')
      amin = min(min(Zi));
      amax = max(max(Zi));
    else
      error('unknown ''maplimits'' value.');
    end
  elseif length(MAPLIMITS) == 2
    amin = MAPLIMITS(1);
    amax = MAPLIMITS(2);
  else
    error('unknown ''maplimits'' value');
  end
  delta = xi(2)-xi(1); % length of grid entry

 end % if ~strcmpi(STYLE,'grid')
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%% Scale the axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %cla  % clear current axis
  hold on
  h = gca; % uses current axes

                          % instead of default larger AXHEADFAC 
  if squeezefac<0.92 && plotrad-headrad > 0.05  % (size of head in axes)
    AXHEADFAC = 1.05;     % do not leave room for external ears if head cartoon
                          % shrunk enough by the 'skirt' option
  end

  set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC);
                          % specify size of head axes in gca

  unsh = (GRID_SCALE+1)/GRID_SCALE; % un-shrink the effects of 'interp' SHADING

  %
  %%%%%%%%%%%%%%%%%%%%%%%% Plot grid only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if strcmpi(STYLE,'grid')                     % plot grid only

    %
    % The goal below is to make the grid cells square - not yet achieved in all cases? -sm
    %
    g1 = size(gridchans,1); 
    g2 = size(gridchans,2); 
    gmax = max([g1 g2]);
    Xi = linspace(-rmax*g2/gmax,rmax*g2/gmax,g1+1);
    Xi = Xi+rmax/g1; Xi = Xi(1:end-1);
    Yi = linspace(-rmax*g1/gmax,rmax*g1/gmax,g2+1);
    Yi = Yi+rmax/g2; Yi = Yi(1:end-1); Yi = Yi(end:-1:1); % by trial and error!
    %
    %%%%%%%%%%% collect the gridchans values %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    gridvalues = zeros(size(gridchans));
    for j=1:size(gridchans,1)
      for k=1:size(gridchans,2)
         gc = gridchans(j,k);
         if gc > 0
              gridvalues(j,k) = Values(gc);
         elseif gc < 0
              gridvalues(j,k) = -Values(abs(gc));
         else 
              gridvalues(j,k) = nan; % not-a-number = no value
         end
      end
    end
    %
    %%%%%%%%%%% reset color limits for grid plot %%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ischar(MAPLIMITS) 
      if strcmp(MAPLIMITS,'maxmin') || strcmp(MAPLIMITS,'minmax')
        amin = min(min(gridvalues(~isnan(gridvalues))));
        amax = max(max(gridvalues(~isnan(gridvalues))));
      elseif strcmp(MAPLIMITS,'absmax')
        % 11/21/2005 Toby edit
        % This should now work as specified. Before it only crashed (using
        % "plotgrid" and "maplimits>absmax" options).
        amax = max(max(abs(gridvalues(~isnan(gridvalues)))));
        amin = -amax;
        %amin = -max(max(abs([amin amax])));
        %amax = max(max(abs([amin amax])));
      else
        error('unknown ''maplimits'' value');
      end
    elseif length(MAPLIMITS) == 2
      amin = MAPLIMITS(1);
      amax = MAPLIMITS(2);
    else
      error('unknown ''maplimits'' value');
    end
    %
    %%%%%%%%%% explicitly compute grid colors, allowing BACKCOLOR  %%%%%%
    %
    gridvalues = 1+floor(cmaplen*(gridvalues-amin)/(amax-amin));
    gridvalues(find(gridvalues == cmaplen+1)) = cmaplen;
    gridcolors = zeros([size(gridvalues),3]);
    for j=1:size(gridchans,1)
      for k=1:size(gridchans,2)
         if ~isnan(gridvalues(j,k))
             gridcolors(j,k,:) = cmap(gridvalues(j,k),:);
         else
            if strcmpi(whitebk,'off')
                gridcolors(j,k,:) = BACKCOLOR; % gridchans == 0 -> background color
                % This allows the plot to show 'space' between separate sub-grids or strips
            else % 'on'
                gridcolors(j,k,:) = [1 1 1]; BACKCOLOR; % gridchans == 0 -> white for printing
            end
         end
      end
    end

    %
    %%%%%%%%%% draw the gridplot image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    handle=imagesc(Xi,Yi,gridcolors); % plot grid with explicit colors
    axis square
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Plot map contours only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'contour')                     % plot surface contours only
    [cls chs] = contour(Xi,Yi,ZiC,CONTOURNUM,'k'); 
    handle = chs;                                   % handle to a contourgroup object
    % for h=chs, set(h,'color',CCOLOR); end
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Else plot map and contours %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'both')  % plot interpolated surface and surface contours
      if strcmp(SHADING,'interp')
       tmph = surface(Xi*unsh,Yi*unsh,zeros(size(Zi))-0.1,Zi,...
               'EdgeColor','none','FaceColor',SHADING);                    
    else % SHADING == 'flat'
       tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi))-0.1,Zi,...
               'EdgeColor','none','FaceColor',SHADING);                    
    end
    if strcmpi(MASKSURF, 'on')
        set(tmph, 'visible', 'off');
        handle = tmph;
    end
    
    warning off;
    if ~PMASKFLAG
        [cls chs] = contour(Xi,Yi,ZiC,CONTOURNUM,'k'); 
    else
        ZiC(find(ZiC > 0.5 )) = NaN;
        [cls chs] = contourf(Xi,Yi,ZiC,0,'k');
        subh = get(chs, 'children');
        for indsubh = 1:length(subh)
            numfaces = size(get(subh(indsubh), 'XData'),1); 
            set(subh(indsubh), 'FaceVertexCData', ones(numfaces,3), 'Cdatamapping', 'direct', 'facealpha', 0.5, 'linewidth', 2);
        end
    end
    handle = tmph;                                   % surface handle
    try, for h=chs, set(h,'color',CCOLOR); end, catch, end % the try clause is for Octave
    warning on;
  %
  %%%%%%%%%%%%%%%%%%%%%%%% Else plot map only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'straight') || strcmp(STYLE,'map') % 'straight' was former arg

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
    end
    handle = tmph;                                   % surface handle
  %
  %%%%%%%%%%%%%%%%%% Else fill contours with uniform colors  %%%%%%%%%%%%%%%%%%
  %
  elseif strcmp(STYLE,'fill')
    [cls chs] = contourf(Xi,Yi,Zi,CONTOURNUM,'k');
    
    handle = chs;                                   % handle to a contourgroup object

    % for h=chs, set(h,'color',CCOLOR); end 
    %     <- 'not line objects.' Why does 'both' work above???

  else
    error('Invalid style')
  end
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set color axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
%   caxis([amin amax]); % set coloraxis

% 7/30/2014 Ramon: +-5% for the color limits were added
cax_sgn = sign([amin amax]);                                                  % getting sign
caxis([amin+cax_sgn(1)*(0.05*abs(amin)) amax+cax_sgn(2)*(0.05*abs(amax))]);   % Adding 5% to the color limits

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
  %cla
  hold on

  set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC)
   % pos = get(gca,'position');
   % fprintf('Current axes size %g,%g\n',pos(3),pos(4));

  if strcmp(ELECTRODES,'labelpoint') ||  strcmp(ELECTRODES,'numpoint')
    text(-0.6,-0.6, ...
    [ int2str(length(Rd)) ' of ' int2str(length(tmpeloc)) ' electrode locations shown']); 
    text(-0.6,-0.7, [ 'Click on electrodes to toggle name/number']);
    tl = title('Channel locations');
    set(tl, 'fontweight', 'bold');
  end
end % STYLE 'blank'

if exist('handle') ~= 1
    handle = gca;
end

if ~strcmpi(STYLE,'grid')                     % if not plot grid only

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
  if cnvfac < 1, cnvfac=1; end
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

if ~ischar(HEADCOLOR) || ~strcmpi(HEADCOLOR,'none')
   %ringh= patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'linewidth', HLINEWIDTH); hold on
   headx = [rx(:)' rx(1)]*hin;
   heady = [ry(:)' ry(1)]*hin;
   ringh= plot(headx,heady);
   set(ringh, 'color',HEADCOLOR,'linewidth', HLINEWIDTH); hold on
end

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
  if ~ischar(HEADCOLOR) || ~strcmpi(HEADCOLOR,'none')
    plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
         2*ones(size([basex;tiphw;0;-tiphw;-basex])),...
         'Color',HEADCOLOR,'LineWidth',HLINEWIDTH);                 % plot nose
    plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)    % plot left ear
    plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)   % plot right ear
  end
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

axis equal;
lim = [-0.525 0.525];
%lim = [-0.56 0.56];
set(gca, 'xlim', lim); set(plotax, 'xlim', lim);
set(gca, 'ylim', lim); set(plotax, 'ylim', lim);
set(gca, 'xlim', lim); set(plotax, 'xlim', lim);
set(gca, 'ylim', lim); set(plotax, 'ylim', lim);
 
%get(textax,'pos')    % test if equal!
%get(plotax,'pos')
%get(textax,'xlim')
%get(plotax,'xlim')
%get(textax,'ylim')
%get(plotax,'ylim')

 if isempty(EMARKERSIZE)
   EMARKERSIZE = 10;
   if length(y)>=160
    EMARKERSIZE = 3;
   elseif length(y)>=128
    EMARKERSIZE = 3;
   elseif length(y)>=100
    EMARKERSIZE = 3;
   elseif length(y)>=80
    EMARKERSIZE = 4;
   elseif length(y)>=64
    EMARKERSIZE = 5;
   elseif length(y)>=48
    EMARKERSIZE = 6;
   elseif length(y)>=32 
    EMARKERSIZE = 8;
   end
 end
%
%%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations only %%%%%%%%%%%%%%%%%%%%%%%%%%
%
ELECTRODE_HEIGHT = 2.1;  % z value for plotting electrode information (above the surf)

if strcmp(ELECTRODES,'on')   % plot electrodes as spots
  if isempty(EMARKER2CHANS)
    hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
  else % plot markers for normal chans and EMARKER2CHANS separately
    hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
        EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
  end
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
  if isempty(EMARKER2CHANS)
    hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
  else
    hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
        EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
  end
  for i = 1:size(labels,1)
    hh(i) = text(double(y(i)+0.01),double(x(i)),...
        ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','left',...
	'VerticalAlignment','middle','Color', ECOLOR,'userdata', num2str(allchansind(i)), ...
	'FontSize',EFSIZE, 'buttondownfcn', ...
	    ['tmpstr = get(gco, ''userdata'');'...
	     'set(gco, ''userdata'', get(gco, ''string''));' ...
	     'set(gco, ''string'', tmpstr); clear tmpstr;'] );
  end
%
%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus numbers %%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'numpoint') 
  if isempty(EMARKER2CHANS)
    hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
  else
    hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
        EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
        EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
  end
  for i = 1:size(labels,1)
    hh(i) = text(double(y(i)+0.01),double(x(i)),...
        ELECTRODE_HEIGHT,num2str(allchansind(i)),'HorizontalAlignment','left',...
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
        ELECTRODE_HEIGHT,int2str(allchansind(i)),'HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',ECOLOR,...
	'FontSize',EFSIZE)
  end
%
%%%%%%%%%%%%%%%%%%%%%% Mark emarker2 electrodes only  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif strcmp(ELECTRODES,'off') && ~isempty(EMARKER2CHANS)
    hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
        EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
end
%
%%%%%%%% Mark specified electrode locations with red filled disks  %%%%%%%%%%%%%%%%%%%%%%
%
try,
    if strcmpi(STYLE,'blank') % if mark-selected-channel-locations mode
        for kk = 1:length(1:length(x))
            if abs(Values(kk))
                if strcmpi(PLOTDISK, 'off')
                    angleRatio = real(Values(kk))/(real(Values(kk))+imag(Values(kk)))*360;
                    radius     = real(Values(kk))+imag(Values(kk));
                    allradius  = [0.02 0.03 0.037 0.044 0.05];
                    radius     = allradius(radius);
                    hp2 = disk(y(kk),x(kk),radius, [1 0 0], 0 , angleRatio, 16);
                    if angleRatio ~= 360
                        hp2 = disk(y(kk),x(kk),radius, [0 0 1], angleRatio, 360, 16);
                    end
                else
                    tmpcolor = COLORARRAY{max(1,min(Values(kk), length(COLORARRAY)))};
                    hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', tmpcolor, 'markersize', EMARKERSIZE1CHAN);
                    hp2 = disk(y(kk),x(kk),real(Values(kk))+imag(Values(kk)), tmpcolor, 0, 360, 10);
                end
            end
        end
    end
catch, end
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
        DIPOLE(:,5) =  tmp.momxyz(:,3);
    else
        DIPOLE(:,1) = -tmp(:,2);                    % same for vector input
        DIPOLE(:,2) =  tmp(:,1);
        DIPOLE(:,3) = -tmp(:,4);
        DIPOLE(:,4) =  tmp(:,3);
    end
    for index = 1:size(DIPOLE,1)
        if ~any(DIPOLE(index,:))
             DIPOLE(index,:) = [];
        end
    end
    DIPOLE(:,1:4)   = DIPOLE(:,1:4)*rmax*(rmax/plotrad); % scale radius from 1 -> rmax (0.5)
    DIPOLE(:,3:end) = (DIPOLE(:,3:end))*rmax/100000*(rmax/plotrad); 
    if strcmpi(DIPNORM, 'on')
        for index = 1:size(DIPOLE,1)
            DIPOLE(index,3:4) = DIPOLE(index,3:4)/norm(DIPOLE(index,3:end))*0.2;
        end
    elseif strcmpi(DIPNORMMAX, 'on')     
        for inorm =  1: size(DIPOLE,1)
            normtmp(inorm) = norm(DIPOLE(inorm,3:5)); % Max norm of projection on XY
        end
        [maxnorm,maxnormindx] = max(normtmp);
        for index = 1:size(DIPOLE,1)
            DIPOLE(index,3:4) = DIPOLE(index,3:4)/norm(DIPOLE(index,3:4))*0.2*normtmp(index)/normtmp(maxnormindx);
        end; 
    end
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
                   [DIPOLE(index, 2) DIPOLE(index, 2)+DIPOLE(index, 4)]',[10 10]);
        set(hh, 'color', DIPCOLOR, 'linewidth', DIPSCALE*30/7);
      end
    end
end

end % if ~ 'gridplot'

%
%%%%%%%%%%%%% Plot axis orientation %%%%%%%%%%%%%%%%%%%%
%
if strcmpi(DRAWAXIS, 'on')
    axes('position', [0 0.85 0.08 0.1]);
    axis off;
    coordend1 = sqrt(-1)*3;
    coordend2 = -3;
    coordend1 = coordend1*exp(sqrt(-1)*rotate);
    coordend2 = coordend2*exp(sqrt(-1)*rotate);
    
    line([5 5+round(real(coordend1))]', [5 5+round(imag(coordend1))]', 'color', 'k');
    line([5 5+round(real(coordend2))]', [5 5+round(imag(coordend2))]', 'color', 'k');
    if round(real(coordend2))<0
         text( 5+round(real(coordend2))*1.2, 5+round(imag(coordend2))*1.2-2, '+Y');
    else text( 5+round(real(coordend2))*1.2, 5+round(imag(coordend2))*1.2, '+Y');
    end
    if round(real(coordend1))<0
         text( 5+round(real(coordend1))*1.2, 5+round(imag(coordend1))*1.2+1.5, '+X');
    else text( 5+round(real(coordend1))*1.2, 5+round(imag(coordend1))*1.2, '+X');
    end
    set(gca, 'xlim', [0 10], 'ylim', [0 10]);
end

%
%%%%%%%%%%%%% Set EEGLAB background color to match head border %%%%%%%%%%%%%%%%%%%%%%%%
%
try, 
  set(gcf, 'color', BACKCOLOR); 
  catch, 
end; 

hold off
axis off
return

% 
% X(2:size(X,1)-1,2:size(X,2)-1) = NaN;
% X(isnan(X(:))) = [];
% X(1:2:end) = [];
% X(1:2:end) = [];
% X(1:2:end) = [];

function vq = gdatav4(x,y,v,xq,yq)
%GDATAV4 MATLAB 4 GRIDDATA interpolation

%   Reference:  David T. Sandwell, Biharmonic spline
%   interpolation of GEOS-3 and SEASAT altimeter
%   data, Geophysical Research Letters, 2, 139-142,
%   1987.  Describes interpolation using value or
%   gradient of value in any dimension.

xy = x(:) + 1i*y(:);

% Determine distances between points
d = abs(xy - xy.');

% Determine weights for interpolation
g = (d.^2) .* (log(d)-1);   % Green's function.
% Fixup value of Green's function along diagonal
g(1:size(d,1)+1:end) = 0;
weights = g \ v(:);

[m,n] = size(xq);
vq = zeros(size(xq));
xy = xy.';

% Evaluate at requested points (xq,yq).  Loop to save memory.
for i=1:m
    for j=1:n
        d = abs(xq(i,j) + 1i*yq(i,j) - xy);
        g = (d.^2) .* (log(d)-1);   % Green's function.
        % Value of Green's function at zero
        g(d==0) = 0;
        vq(i,j) = g * weights;        
    end
end

%
%%%%%%%%%%%%% Draw circle %%%%%%%%%%%%%%%%%%%%%%%%
%
function h2 = disk(X, Y, radius, colorfill, oriangle, endangle, segments)
	A = linspace(oriangle/180*pi, endangle/180*pi, segments-1);
    if endangle-oriangle == 360
     	 A  = linspace(oriangle/180*pi, endangle/180*pi, segments);
         h2 = patch( [X   + cos(A)*radius(1)], [Y   + sin(A)*radius(end)], zeros(1,segments)+3, colorfill);
    else A  = linspace(oriangle/180*pi, endangle/180*pi, segments-1);
         h2 = patch( [X X + cos(A)*radius(1)], [Y Y + sin(A)*radius(end)], zeros(1,segments)+3, colorfill);
    end
    set(h2, 'FaceColor', colorfill);
	set(h2, 'EdgeColor', 'none'); 
