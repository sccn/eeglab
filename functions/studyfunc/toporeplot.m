% toporeplot() - re-plot a saved topoplot() output image (a square matrix) 
%              in a 2-D circular scalp map view (as looking down at the top 
%              of the head).  May also be used to re-plot a mean topoplot() 
%              map for a number of subjects and/or components without all 
%              the constitutent maps having the same channel montage.
%              Nose is at top of plot. Left = left. See topoplot().
% Usage:
%          >>  toporeplot(topoimage,'plotrad',val1, 'intrad',val2);   
%                 % Use an existing (or mean) topoplot output image. Give the 
%                 % original 'intrad' and 'plotrad' values used to in topoimage.
%          >>  toporeplot(topoimage,'plotrad',val1,'xsurface', Xi, 'ysurface',Yi );   
%                 % Use an existing (or mean) topoplot output image. Give the same 
%                 % 'plotrad' value used to create it. Since 'intrad' was not input 
%                 % to topoplot(), give the grid axes Xi and Yi as inputs.
%          >>  [hfig val ]= toporeplot(topoimage,'plotrad',val1, 'Param1','Value1', ...); 
%                 % Give one of the two options above plus other optional parameters.
% Required inputs:
%   topoimage    - output image matrix (as) from topoplot(). For maximum flexibility, 
%                  create topoimage using topoplot() options: 'plotrad',1, 'intrad',1 
%   'plotrad'    - [0.15<=float<=1.0] plotting radius = max channel arc_length to plot.
%                  If plotrad > 0.5, chans with arc_length > 0.5 (i.e. below ears,eyes) 
%                  are plotted in a circular 'skirt' outside the cartoon head. 
%                  The topoimage depends on 'plotrad', so 'plotrad' is required to 
%                  reproduce the 'topoplot' image. 
% Optional inputs:
%   'chanlocs'   - name of an EEG electrode position file (see >> topoplot example).
%                  Else, an EEG.chanlocs structure        (see >> help pop_editset).
%   'maplimits'  - 'absmax'   -> scale map colors to +/- the absolute-max (makes green 0); 
%                  'maxmin'   -> scale colors to the data range (makes green mid-range); 
%                  [lo,hi]    -> use user-definined lo/hi limits {default: 'absmax'}
%   'style'      - 'map'      -> plot colored map only
%                  'contour'  -> plot contour lines only
%                  'both'     -> plot both colored map and contour lines {default: 'both'}
%                  'fill'     -> plot constant color between contour lines
%                  'blank'    -> plot electrode locations only, requires electrode info. 
%   'electrodes' - 'on','off','labels','numbers','ptslabels','ptsnumbers' See Plot detail 
%                  options below. {default: 'on' -> mark electrode locations with points
%                  unless more than 64 channels, then 'off'}. Requires electrode info. 
%   'intrad'     - [0.15<=float<=1.0] radius of the interpolation area used in topoplot() 
%                  to get the grid. 
%   'headrad'    - [0.15<=float<=1.0] drawing radius (arc_length) for the cartoon head. 
%                  NB: Only headrad = 0.5 is anatomically correct! 0 -> don't draw head; 
%                  'rim' -> show cartoon head at outer edge of the plot {default: 0.5}. 
%                  Requires electrode information. 
%   'noplot'     - [rad theta] are coordinates of a (possibly missing) channel. 
%                  Do not plot but return interpolated value for channel location.
%                  Do not plot but return interpolated value for this location. 
%   'xsurface'   - [Xi- matrix] the Xi grid points for the surface of the plotting
%                  an output of topoplot().
%   'ysurface'   - [Yi- matrix] the Yi grid points for the surface of the plotting,
%                  an output of topoplot().
% Dipole plotting:
%   'dipole'     - [xi yi xe ye ze] plot dipole on the top of the scalp map
%                  from coordinate (xi,yi) to coordinates (xe,ye,ze) (dipole head 
%                  model has radius 1). If several rows, plot one dipole per row.
%                  Coordinates returned by dipplot() may be used. Can accept
%                  an EEG.dipfit.model structure (See >> help dipplot).
%                  Ex: ,'dipole',EEG.dipfit.model(17) % Plot dipole(s) for comp. 17.
%   'dipnorm'    - ['on'|'off'] normalize dipole length {default: 'on'}.
%   'diporient'  - [-1|1] invert dipole orientation {default: 1}.
%   'diplen'     - [real] scale dipole length {default: 1}.
%   'dipscale'   - [real] scale dipole size {default: 1}.
%   'dipsphere'  - [real] size of the dipole sphere. {default: 85 mm}.
%   'dipcolor'   - [color] dipole color as Matlab code code or [r g b] vector
%                  {default: 'k' = black}.
% Plot detail options:
%   'electcolor' {'k'}|'emarker' {'.'}|'emarkersize' {14} ...
%  |'emarkersize1chan' {40}|'efontsize' {var} - electrode marking details and {defaults}. 
%   'shading'    - 'flat','interp'  {default: 'flat'}
%   'colormap'   - (n,3) any size colormap {default: existing colormap}
%   'numcontour' - number of contour lines {default: 6}
%   'ccolor'     - color of the contours {default: dark grey}
%   'hcolor'|'ecolor' - colors of the cartoon head and electrodes {default: black}
%   'circgrid'   - [int > 100] number of elements (angles) in head and border circles {201}
%   'verbose'    - ['on'|'off'] comment on operations on command line {default: 'on'}.
%
% Outputs:
%         hfig   - plot axes handle
%         val    - single interpolated value at the specified 'noplot' arg channel 
%                   location ([rad theta]).
%
% Notes: - To change the plot map masking ring to a new figure background color,
%            >> set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
%        - Topoplots may be rotated from the commandline >> view([deg 90]) {default:[0 90])
%
% Authors: Hilit Serby, Andy Spydell, Colin Humphries, Arnaud Delorme & Scott Makeig
%          CNL / Salk Institute, 8/1996-/10/2001; SCCN/INC/UCSD, Nov. 2001- Nov. 2004
%
% See also: topoplot(), timtopo(), envtopo()

% Deprecated but still usable;
%   'interplimits'    - ['electrodes'|'head'] 'electrodes'-> interpolate the electrode grid; 
%                       'head'-> interpolate the whole disk {default: 'head'}.

% Copyright (C) UCSD
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

function [handle,chanval] = toporeplot(grid,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10)

%
%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

icadefs                 % read defaults MAXTOPOPLOTCHANS and DEFAULT_ELOC and BACKCOLOR
if ~exist('BACKCOLOR')  % if icadefs.m does not define BACKCOLOR
   BACKCOLOR = [.93 .96 1];  % EEGLAB standard
end
GRID_SCALE = length(grid);       
noplot  = 'off';
handle = [];
chanval = NaN;
rmax = 0.5;             % head radius - don't change this!
INTERPLIMITS = 'head';  % head, electrodes
MAPLIMITS = 'absmax';   % absmax, maxmin, [values]
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
plotrad      = [];      % plotting radius ([] = auto, based on outermost channel location)
intrad       = [];      % default interpolation square is to outermost electrode (<=1.0)
headrad      = [];      % default plotting radius for cartoon head is 0.5
MINPLOTRAD = 0.15;      % can't make a topoplot with smaller plotrad (contours fail)
VERBOSE = 'off';
MASKSURF = 'off';

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
if ~mod(nargs,2)
    error('Optional inputs must come in Key - Val pairs')
end
if ~isnumeric(grid) || size(grid,1) ~= size(grid,2)
    error('topoimage must be a square matrix');
end
for i = 2:2:nargs
	Param = eval(['p',int2str((i-2)/2 +1)]);
	Value = eval(['v',int2str((i-2)/2 +1)]);
	if ~ischar(Param)
      error('Flag arguments must be strings')
	end
	Param = lower(Param);
	switch lower(Param)
        case  'chanlocs' 
            loc_file = Value;
        case 'colormap'
		  if size(Value,2)~=3
              error('Colormap must be a n x 3 matrix')
		  end
		  colormap(Value)
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
                  || strcmpi(ELECTRODES,'labelspts') || strcmpi(ELECTRODES,'ptlabels') ...
                  || strcmpi(ELECTRODES,'labelpts') 
                 ELECTRODES = 'labelpoint'; % backwards compatability
             end
             if strcmpi(ELECTRODES,'pointnumbers') || strcmpi(ELECTRODES,'ptsnumbers') ...
                  || strcmpi(ELECTRODES,'numberspts') || strcmpi(ELECTRODES,'ptnumbers') ...
                  || strcmpi(ELECTRODES,'numberpts')  || strcmpi(ELECTRODES,'ptsnums')  ...
                  || strcmpi(ELECTRODES,'numspts') 
                 ELECTRODES = 'numpoint'; % backwards compatability
             end
             if strcmpi(ELECTRODES,'nums') 
                 ELECTRODES = 'numbers'; % backwards compatability
             end
             if strcmpi(ELECTRODES,'pts') 
                 ELECTRODES = 'on'; % backwards compatability
             end
             if ~strcmpi(ELECTRODES,'labelpoint') ...
                && ~strcmpi(ELECTRODES,'numpoint') ...
                && ~strcmp(ELECTRODES,'on') ...
                && ~strcmp(ELECTRODES,'off') ...
                && ~strcmp(ELECTRODES,'labels') ...
                && ~strcmpi(ELECTRODES,'numbers') 
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
		 case 'plotrad'
		  plotrad = Value;
              if ischar(plotrad) || (plotrad < MINPLOTRAD || plotrad > 1)
		     error('plotrad argument should be a number between 0.15 and 1.0');
		  end
		case 'intrad'
		  intrad = Value;
          if ischar(intrad) || (intrad < MINPLOTRAD || intrad > 1)
              error('intrad argument should be a number between 0.15 and 1.0');
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
        case 'xsurface'
            Xi = Value;
            if ~isnumeric(Xi) || size(Xi,1) ~= size(Xi,2) || size(Xi,1) ~= size(grid,1)
                error('xsurface must be a square matrix the size of grid');
            end
        case  'ysurface'
            Yi = Value;
            if ~isnumeric(Yi) || size(Yi,1) ~= size(Yi,2) || size(Yi,1) ~= size(grid,1)
                error('ysurface must be a square matrix the size of grid');
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
            if ~ischar(noplot)
                if length(noplot) ~= 2
                    error('''noplot'' location should be [radius, angle]')
                else
                    chanrad = noplot(1);
                    chantheta = noplot(2);
                    noplot = 'on';
                end
            end
        otherwise
            error(['Unknown input parameter ''' Param ''' ???'])
	end
end

if isempty(plotrad)
    error(' ''plotrad'' must be given')
end
if isempty(intrad)
    if ~exist('Yi') || ~exist('Xi')
        error('either ''intrad'' or the grid axes (Xi and Yi) must be given');
    end
end

%
%%%%%%%%%%%%%%%%%%%% Read the channel location information %%%%%%%%%%%%%%%%%%%%%%%%
% 
if exist('loc_file')
	if ischar(loc_file)
		[tmpeloc labels Th Rd indices] = readlocs(loc_file,'filetype','loc');
	else % a locs struct
		[tmpeloc labels Th Rd indices] = readlocs(loc_file);
            % Note: Th and Rd correspond to indices channels-with-coordinates only
	end
    labels = strvcat(labels);
    Th = pi/180*Th;                              % convert degrees to radians

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
        
	if ischar(plotrad) || plotrad < MINPLOTRAD || plotrad > 1.0
       error('plotrad must be between 0.15 and 1.0');
	end
end
if isempty(plotrad) && ~ exist('loc_file')
    plotrad = 1;                 % default: plot out to the 0.5 head bounda
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
%%%%%%%%%%%%%%%%% Issue warning if headrad ~= rmax  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

if headrad ~= 0.5 && strcmpi(VERBOSE, 'on')
   fprintf('     NB: Plotting map using ''plotrad'' %-4.3g,',plotrad);
   fprintf(    ' ''headrad'' %-4.3g\n',headrad);
   fprintf('Warning: The plotting radius of the cartoon head is NOT anatomically correct (0.5).\n')
end

squeezefac = rmax/plotrad;

%
%%%%%%%%%%%%%%%%%%%%% Find plotting channels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if exist('tmpeloc')
    pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
    [x,y] = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
	%
	%%%%%%%%%%%%%%%%%%%%% Eliminate channels not plotted  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 
	
	allx = x;
	ally = y;
	
	
	Th    = Th(pltchans);              % eliminate channels outside the plotting area
	Rd    = Rd(pltchans);
	x     = x(pltchans);
	y     = y(pltchans);
	
	labels= labels(pltchans,:);
    
	%
	%%%%%%%%%%%%%%% Squeeze channel locations to <= rmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 
	
	Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
                              % to plot all inside the head cartoon
	x    = x*squeezefac;    
	y    = y*squeezefac;   
	allx    = allx*squeezefac;    
	ally    = ally*squeezefac;   
end	
% Note: Now outermost channel will be plotted just inside rmax

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~strcmpi(STYLE,'blank') % if draw interpolated scalp map
  %
  %%%%%%%%%%%%%%%%%%%%%%% Interpolate scalp map data %%%%%%%%%%%%%%%%%%%%%%%%
  %
  if ~isempty(intrad) % intrad specified
      xi = linspace(-intrad*squeezefac,intrad*squeezefac,GRID_SCALE);   % use the specified intrad value 
      yi = linspace(-intrad*squeezefac,intrad*squeezefac,GRID_SCALE);   
      [Xi,Yi] = meshgrid(yi',xi);
  elseif ~exist('Xi') || ~exist('Yi')
      error('toporeplot require either intrad input or both xsurface and ysurface')
  end
  Zi = grid;
  mask = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
  ii = find(mask == 0);
  Zi(ii) = NaN; 
  
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
  if ischar(MAPLIMITS)
    if strcmp(MAPLIMITS,'absmax')
      amin = -max(max(abs(Zi)));
      amax = max(max(abs(Zi)));
    elseif strcmp(MAPLIMITS,'maxmin') || strcmp(MAPLIMITS,'minmax')
      amin = min(min(Zi));
      amax = max(max(Zi));
    else
      error('unknown ''maplimits'' value.');
    end
  else
    amin = MAPLIMITS(1);
    amax = MAPLIMITS(2);
  end
  delta = Xi(1,2)-Xi(1,1); % length of grid entry

  %
  %%%%%%%%%%%%%%%%%%%%%%%%%% Scale the axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  cla  % clear current axis
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
  
  switch STYLE
      %
      %%%%%%%%%%%%%%%%%%%%%%%% Plot map contours only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      case 'contour' % plot surface contours only
          [cls chs] = contour(Xi,Yi,Zi,CONTOURNUM,'k'); 
      %
      %%%%%%%%%%%%%%%%%%%%%%%% Else plot map and contours %%%%%%%%%%%%%%%%%%%%%%%%%
      %
	  case 'both'   % plot interpolated surface and surface contours
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
          end
          [cls chs] = contour(Xi,Yi,Zi,CONTOURNUM,'k'); 
          for h=chs, set(h,'linecolor',CCOLOR); end
      %
      %%%%%%%%%%%%%%%%%%%%%%%% Else plot map only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      case {'straight', 'map'} % 'straight' was former arg
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
      %
      %%%%%%%%%%%%%%%%%% Else fill contours with uniform colors  %%%%%%%%%%%%%%%%%%
      %
	  case 'fill'
          [cls chs] = contourf(Xi,Yi,Zi,CONTOURNUM,'k');
      otherwise
          error('Invalid style')
  end

  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set color axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  caxis([amin amax]) % set coloraxis
  %
  %%%%%%%%%%%%%%%%%%%%%%% Draw blank head %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  else  % if STYLE 'blank'
      if strcmpi(noplot, 'on') 
          if strcmpi(VERBOSE,'on')
              fprintf('topoplot(): no plot requested.\n')
          end
          return;
      end
      cla
      hold on
	  set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC)
      if ~exist('tmpeloc')
          error('No electrode location information found');
      end
      if strcmp(ELECTRODES,'labelpoint') || strcmp(ELECTRODES,'numpoint')
          text(-0.6,-0.6, [ int2str(length(Rd)) ' of ' int2str(length(tmpeloc)) ' electrode locations shown']);
          text(-0.6,-0.7, [ 'Click on electrodes to toggle name/number']);
          tl = title('Channel locations');
          set(tl, 'fontweight', 'bold');
      end
  end
  
if exist('handle') ~= 1
    handle = gca;
end

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

circ = linspace(0,2*pi,CIRCGRID);
rx = sin(circ); 
ry = cos(circ); 
ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];

if ~strcmpi(STYLE,'blank')
  ringh= patch(ringx,ringy,0.01*ones(size(ringx)),BACKCOLOR,'edgecolor','none'); hold on
end

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

 axis square                                           % make textax square

 pos = get(gca,'position');
 set(plotax,'position',pos);

 xlm = get(gca,'xlim');
 set(plotax,'xlim',xlm);

 ylm = get(gca,'ylim');
 set(plotax,'ylim',ylm);                               % copy position and axis limits again

 %%%%%%%%%%%%%%%%%%%%%%%%%only if electrode info is available  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 if exist('tmpeloc')
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
                   [DIPOLE(index, 2) DIPOLE(index, 2)+DIPOLE(index, 4)]');
        set(hh, 'color', DIPCOLOR, 'linewidth', DIPSCALE*30/7);
      end
    end
end

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
