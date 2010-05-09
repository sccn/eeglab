% brainmovie3d() - generate a sequence of images showing event-related coherence,
%               event-related spectral perturbations, and inter-trial coherence
%               of localized EEG waveforms. Uses outputs of timef() and cross().
% Usage:
%   >> brainmovie3d(ersps,itcs,crossfs_amp,crossfs_phase,times,freqs,selected,...
%                         'keyword1',value1,...); % creates files image0001.eps, etc.
%
% Inputs:
% ersps         - Cell array (components,conditions) of ERSP arrays (freqs,times)
%                 ERSP = event-related spectral perturbation; returned by timef()
% itcs          - Cell array (components,conditions) of ITC arrays (freqs,times)
%                 ITC = inter-trial coherence; returned by timef()  
% crossfs_amp   - Cell array (components,components,conditions) of crossf() 
%                 amplitude output arrays of size (freqs,times).
% crossfs_phase - Cell array (components,components,conditions) of crossf() phase
%                 output arrays of size (freqs,times). (Only the upper diagonal part 
%                 of the matrix is taken into account).
% times         - Array of times returned by timef() or crossf()
% freqs         - Indices into the array of freqs returned by timef() or crossf() 
%                 (e.g., [1:2] means plot the mean of the first two frequencies). 
%                 These indexes determine for which freqs plotting will be performed.
% selected      - Component indices to plot (default all)
%
% Optional 'keyword' parameters:
% 'latency'   - plot only a subset of latencies. The time point closest to the 
%               latency given are plotted. Default = empty, all latencies.
% 'frames'    - vector of frame indices to compute. [1:2] only computes the
%               first two frames.
% 'envelope'  - (2,points,conditions) envelopes of the average data (ERP) in each condition
%               (envelope =  min and max traces of each ERP across all channels and times)
% 'rt'        - cell array of vector containing reaction times of the subject in 
%               each conditions. This will plot a small bar which height will vary
%               based on the probability of response (default {} -> ignored)
% 'flashes'   - vector of time indices at which the background flashes.  Specify the color 
%               of the flash with a cell array of [1,2] cell arrays. 
%               Ex. { { 200 'y' } { 1500 '5' }} will generate two flashes, 
%               yellow at 200 ms and red at 1500 ms 
%
% Movie ITC, Power and Crossf options:
% 'power'     - ['on'|'off'] vary the size of the component disks according to spectral power 
%                                                           {default: on}
% 'itc'       - ['on'|'off'] vary component disk colors according to inter-trial coherence 
%							    {default: on}
% 'crossf'    - ['on'|'off'] plot | do not plot coherence   {default: on}
% 'crossfcoh' - ['on'|'off'] vary the width of the connecting arc 
%                               according to cross-coherence magnitude {def: on}
% 'crossfphasecolor' -['on'|'off'] vary the arc color according to coherence {default: on}
% 'crossfphasespeed' - ['on'|'off'] vary the arc speed according to 
%                                      cross-coherence phase {def: off}
% 'crossfphaseunit'  - ['degree'|'radian']. Coherence phase angle unit {Default is degree}.
% 'colmapcrossf' - colormap array for arcs {default: hsv(64) with green as 0} 
% 'colmapcoh'   - colormap array for disks (according to inter-trial coherence) 
%                      {default: hot(64)}
% 'scalepower'  - [min max] dB range for power (and disk size) variation {default: [-5 5]}  
% 'scalecoher'  - [min max] coherence range {default: [0 1]}
% 'scaleitc'    - [absmax] maximum itc {Default: 1}
% 'polarity'  - ['pos'|'posneg'] polarity for ITC and crossf. 'pos' = only positive values
%               'posneg' = positive and negative values.
%
% Movie coordinates and axis options:
% 'magnify'   - integer magnification factor for graphics. Default is 1.
% 'diskscale'   - numeric value that scales the size of disks {default: [1.0]}   
% 'xlimaxes'    - x-axis limits axis for the component locations {default: [-1 1]}
% 'ylimaxes'    - y-axis limits axis for the component locations {default: [-1 to 1]}
% 'coordinates' - 2-column array of [x y] coordinates of the selected components 
%                 {default: spaced evenly around the head circle boundary}  
% 'square'    - ['on'|'off'] re-square all coordinates (so X and Y width is the same)
%               default is 'on';
% 'project3d' - ['on'|'off'] project disks on each 3-D axis. Default is 'off'.
% 'circfactor'  - (ncomps,ncomps) array of arc curvatures (0=straight; 1=half-round, 
%                 positive or negative values give the sense of rotation) {def: 0s}
% 'envylabel'   - ordinate label for envelope. {Default 'Potential \muV'}
% 'envvert'     - cell array of time indices at which to draw vertical lines.
%                 Can also be a cell array of cell to specify line aspect. For instance
%                 { { 0 'color' 'b' 'linewidth' 2 } {1000 'color' 'r' }} would draw two
%                 lines, one blue thick line at latency 0 and one thin red line at latency 1000.
% 'rthistloc' - location and size of rt histograms in individual axes. 
%               [abscissa ordinate width maxheight].
% 'title'       - (string) main movie title
% 'condtitle'   - (string array) condition titles (one condition title per row)
% 'condtitleformat' - list of title properties. Ex: { 'fontize', 12, 'fontweight', 'bold' }
% 'plotorder'   - [integer vector] component plot order from 1 to the number of selected 
%                 components. 
% 'backcolor' - [float array] background color. Default is [1 1 1] (white).
%
% Picture and movie output options:
% 'moviename'  - ['string'] Movie file name. Default is "output.avi".
% 'movieopts'  - [cell] Movie options for avifile function. See "help avifile".
% 'framesout'  - ['eps'|'ppm'|'fig'|'tiff'|'none'] Default format for saving frames on disk. 
%                Default is 'tiff'.
% 'framefolder' - [string] frames output folder. Default uses current directory.
%               the directory is created if it does not exist.
% 'visible'    - ['on'|'off'] show the images on the screen or keep them hidden {default 'on'}
% 'size'      - [widthcond height] output image size {default [400,400]}
%               widthcond is the width of a single condition plot (in pixels)
% 'view'      - 3D static starting view.  See help view. Default is [1 0 0].
% 'path3d'    - ['on'|'off'|[thetafact phifact]] 'on' activate automatic rotation in 3-D. Use
%               [exttheta extphi] to specify theta and phi multiplicative factor (default is
%               [1 0.75]. Use parameter 'view' to specify starting view point. Default is
% 'stereo'    - [Real] Create a stereo movie. The figure should contain a [left right]
%               display of two identical 3-D plots. The left plot view will follow the 
%               given 'path' (see above). The right plot axis will be 3-D rotated by an 
%               additional horizontal disparity angle specified by the 'stereo' argument:
%               6 (degrees) suggested. Default is [] = mono display.
%               'off'.
%Outputs to disk:
% imageX      - brainmovie3d() saves an output.avi movie (see 'moviename' option above)
%               and a sequence of image files to disk (image0001.eps, as define in the 
%               'framesout' option).
%Example:
%
% % Given ICA activations in array icaact (size ncomps,nframes,ntrials), animate (here) 
% % activity at/between two components at 176 points per epoch (from -100 ms to 600 ms 
% % re stimulus onset) assuming a 250-Hz sampling rate and 100 output frames
%
% >> [ersps{1,1},itcs{1,1},powbase,times,freqs] = ...                          % timef for
%                timef(icaact(1,:),176,[-100 600],'Component
%                1',250,1,32,100); %     1st comp
% >> [ersps{2,1},itcs{2,1},powbase,times,freqs] = ...                          % timef for
%                timef(icaact(2,:),176,[-100 600],'Component 2',250,1,32,100); %     2nd comp
% >> [crossfs_amp{1,2},mcoh,times,freqs,cohboot,crossfs_phase{1,2}] = ...      % crossf for
%      crossf_(icaact(1,:),icaact(2,:),176,[-100 600],'Crossf 1 and 2',250,1,32,100); % both
%
% >> brainmovie3d( ersps, itcs, crossfs_amp, crossfs_phase, times, [1:2] );
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 30 Mai 2003
%
% Note: Better resolution movies can be generated by .eps -> .ppm -> .avi, 
%       (or, under a planned upgrade to brainmovie3d, from Matlab6 to .avi directly).
% >> !/usr/local/bin/convert images*.eps movie.mpg % ImageMagic 'convert' may be
%                                                  % used to generate the movie.
   
% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2003

% This program is free software; you can redistribute it and/or
% modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function [alltimepoints mov] = brainmovie3d(ALLERSP,ALLITC,ALLCROSSF,ALLCROSSFANGLE,times,FREQS,selected,varargin);
    
if nargin < 6
	help brainmovie3d;
	return;
end;	

% create structure for option if necessary
%-----------------------------------------
if ~isempty( varargin ), 
	for index=1:length(varargin)
		if iscell(varargin{index})
			varargin{index} = { varargin{index}};
		end;
	end;
	g=struct(varargin{:}); 
else 
	g= []; 
end;

if nargin < 7
	selected = 1:size(ALLERSP, 1);
end;

nbconditions = size(ALLERSP,2);
nbcomponents = size(ALLERSP,1);

% add defaults
%-------------
defmaxpow = 0;
defmaxitc = 0;
for i=1:length(ALLERSP)
    defmaxpow = max(defmaxpow, max(abs(ALLERSP{  i}(:))));
    defmaxitc = max(defmaxitc, max(abs(ALLITC{   i}(:))));
end;    
defmaxcoh = 0;
for i=1:length(ALLCROSSF(:))
    tmpmax = max(abs(ALLCROSSF{i}(:)));
    if ~isempty(tmpmax), defmaxcoh = max(defmaxcoh, tmpmax); end;
end;
defmaxpow = ceil(defmaxpow*100)/100;
defmaxitc = ceil(defmaxitc*100)/100;
defmaxcoh = ceil(defmaxcoh*100)/100;
try, g.head; 			catch, g.head=''; end;
try, g.visible; 		catch, g.visible='on'; end;
try, g.square; 		    catch, g.square='on'; end;defmaxpow = ceil(defmaxpow*100)/100;
try, g.moviename; 	    catch, g.moviename='output.avi'; end;
try, g.movieopts; 	    catch, g.movieopts={}; end;
try, g.rt; 	        	catch, g.rt={}; end;
try, g.power; 	    	catch, g.power='on'; end;
try, g.latency; 	   	catch, g.latency=[]; end;
try, g.itc; 	    	catch, g.itc='on'; end;
try, g.magnify; 	    catch, g.magnify=1; end;
try, g.crossf; 			catch, g.crossf='on'; end;
try, g.crossfcoh; 		catch, g.crossfcoh='on'; end;
try, g.size; 			catch, g.size=[400 400]; end;
try, g.crossfphasecolor;catch, g.crossfphasecolor='on'; end;
try, g.crossfphasespeed;catch, g.crossfphasespeed='off'; end;
try, g.crossfphaseunit; catch, g.crossfphaseunit='degree'; end;
try, g.scalepower;      catch, g.scalepower = [-defmaxpow defmaxpow]; end;
try, g.scalecoher;      catch, g.scalecoher = [0 defmaxcoh]; end;
try, g.scaleitc;        catch, g.scaleitc = defmaxitc; end;
try, g.diskscale;       catch, g.diskscale = 1; end;
try, g.framefolder;     catch, g.framefolder = ''; end;
try, g.envelope;        catch, g.envelope = []; end; 
try, g.caption;			catch, g.caption = 'on'; end; 
try, g.frames;			catch, g.frames = []; end; 
try, g.envvert;			catch, g.envvert = {}; end; 
try, g.flashes;			catch, g.flashes = []; end; 
try, g.polarity;		catch, g.polarity = 'pos'; end; 
try, g.framesout;	    catch, g.framesout = 'tiff'; end; 
try, g.condtitle;		catch, g.condtitle = []; end; 
try, g.condtitleformat;	catch, g.condtitleformat = {'fontsize', 14', 'fontweight', 'bold' }; end;
try, g.title;			catch, g.title = []; end; 
try, g.envylabel;		catch, g.envylabel = 'Potential \muV'; end; 
try, g.plotorder;       catch, g.plotorder = selected; end;
try, g.coordformat;     catch, g.coordformat = 'spherical'; end;
try, g.stereo;          catch, g.stereo = []; end;
try, g.backcolor;       catch, g.backcolor = [0 0 0]; end;
try, g.path3d;          catch, g.path3d = 'off'; end;
try, g.project3d;       catch, g.project3d = 'off'; end;
try, g.view;            catch, g.view = [43.6650 30.4420]; end;
try, g.colmapcoh;       catch, 
    colormtmp = hot(64);
    colormtmp(end,3) = (colormtmp(end,3)+colormtmp(end-1,3))/2; % white does not come out when the
    g.colmapcoh = colormtmp;                                    % the figure is printed to ppm
    g.colmapcoh(:,1) =  colormtmp(:,2);
    g.colmapcoh(:,2) =  colormtmp(:,3);
    g.colmapcoh(:,3) =  colormtmp(:,1);
    g.colmapcoh = [ g.colmapcoh; colormtmp(end:-1:1,:)];
    g.colmapcoh = jet(64);
end; 
try, g.colmapcrossf; catch,
    g.colmapcrossf = jet(64);
	%g.colmapcrossf = hsv(64); 
	%g.colmapcrossf = [ g.colmapcrossf(55:end,:); 
	%g.colmapcrossf(1:54,:)]; g.colmapcrossf = g.colmapcrossf(linspace(64, 1, 64),:); % reorganize the colormap
	%g.colmapcrossf = hsv(64);
	%g.colmapcrossf = [g.colmapcrossf(16:end,:); g.colmapcrossf(1:5,:)];
end;
try, g.xlimaxes; 		catch, g.xlimaxes = [-1 1]; end;  
try, g.ylimaxes; 		catch, g.ylimaxes = [-1 1]; end;  
try, g.rthistloc; 	    catch, g.rthistloc(1) = (g.xlimaxes(2)-g.xlimaxes(1))*0.74 + g.xlimaxes(1); % abscicia
	                           g.rthistloc(3) = (g.xlimaxes(2)-g.xlimaxes(1))*0.1; % width
	                           g.rthistloc(2) = (g.ylimaxes(2)-g.ylimaxes(1))*0.34 + g.ylimaxes(1); % ordinate
	                           g.rthistloc(4) = (g.ylimaxes(2)-g.ylimaxes(1))*0.1; % max height
end;
try, g.coordinates; catch,    
    % coordinates around a circle
    g.coordinates = zeros( nbcomponents, 2 );
    count = 0;
   	for index = selected
    	if length(selected) > 1
   			g.coordinates( index,:) = [ cos(count/length(selected)*2*pi) sin(count/length(selected)*2*pi) ] * 0.7;
    	else	g.coordinates(index,:) = [ 0.01 0.01];
		end;
		count = count + 1;
    end;
end;
try, g.circfactor; catch, g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end;
if isempty(g.circfactor), g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end;
if isstr(g.path3d)
    switch g.path3d
     case 'on', g.path3d = [ 1 0.75];
     case 'off', g.path3d = [];
    end;
else
    if length(g.path3d) ~= 2, error('path3d length have to be a string or a 2 element vector'); end;
end;  

% messages for defaults
% ---------------------
if g.scalepower(1) == -defmaxpow, fprintf('Power limits set to %1.2f to %1.2f dB\n', -defmaxpow, defmaxpow); end;
if g.scaleitc(1)   ==  defmaxitc, fprintf('ITC limits set to 0 to %1.2f\n', defmaxitc); end;
if g.scalecoher(2) ==  defmaxcoh, fprintf('Coherence limits set to 0 to %1.2f\n', defmaxcoh); end;

% check size of inputs
% --------------------
try
	if ~all(size(ALLERSP) == size(ALLITC))
		disp('Error: ERSP and ITC cells array must be the same size'); return;
	end;	
	if ~isempty(ALLCROSSF)
		if ~all(size(ALLCROSSF) == size(ALLCROSSFANGLE))
			disp('Error: Crossf amplitude and Crossf angle cells array must be the same size'); return;
		end;	
		if ~(size(ALLCROSSF,2) == size(ALLERSP,1))
			disp('Error: number of components different in ERSP and Crossf arrays'); return;
		end;	
		if ~(size(ALLCROSSF,3) == size(ALLERSP,2))
			disp('Error: number of conditions different in ERSP and Crossf arrays'); return;
		end;	
		if ~(size(ALLCROSSF{1,2,1},1) == size(ALLERSP{1,1},1))
			disp('Error: number of frequencies (rows) different in ERSP and Crossf arrays'); return;
		end;	
		if ~(size(ALLCROSSFANGLE{1,2,1},2) == size(ALLITC{1,1},2))
			disp('Error: number of time points (columns) different in ERSP and Crossf arrays'); return;
		end;	
		if ~(size(ALLCROSSF{1,2,1},2) == length(times))
			disp('Error: number of time points (columns) different in times and Crossf arrays'); return;
		end;
	end;
	try, tmp = ALLERSP{1,1}; tmp(FREQS,:); catch, disp('Error: unable to access the defined frequencies in ERSPs (out of bounds) '); return; end;
	try, ALLERSP{selected,1}; catch, disp('Error: unable to access the defined components in ERSPs (out of bounds)'); return; end;
catch
	disp('Error accessing one of the variable. Remember: Except for selected, freqs, times and circfactor, all vars are cell arrays. Check also: dimensions and content.'); return;
end;	 

% check structure content
% -----------------------
if ~isempty(g.rt)
	if length(g.rt) ~= nbconditions
		disp('Error: Rt must be either an array of the size of the number of conditions (might be 0 for some conditions)'); return;
	end;
end;	
switch lower(g.visible)
	case {'on', 'off'} ;  
	otherwise disp('Error: Visibility must be either ''on'' or ''off'''); return;
end;	
switch lower(g.square)
	case {'on', 'off'} ;  
	otherwise disp('Error: Square must be either ''on'' or ''off'''); return;
end;	
switch lower(g.power)
	case {'on', 'off'} ;  
	otherwise disp('Error: Power must be either ''on'' or ''off'''); return;
end;	
switch lower(g.itc)
	case {'on', 'off'} ;  
	otherwise disp('Error: Itc must be either ''on'' or ''off'''); return;
end;	
switch lower(g.crossf)
	case {'on', 'off'} ;  
	otherwise disp('Error: Crossf must be either ''on'' or ''off'''); return;
end;	
switch lower(g.crossfcoh)
	case {'on', 'off'} ;  
	otherwise disp('Error: Crossfcoh must be either ''on'' or ''off'''); return;
end;	
switch lower(g.crossfphasecolor)
	case {'on', 'off'} ;  
	otherwise disp('Error: Crossfphasecolor must be either ''on'' or ''off'''); return;
end;	
switch lower(g.crossfphasespeed)
	case {'on', 'off'} ;  
	otherwise disp('Error: Crossfphasespeed must be either ''on'' or ''off'''); return;
end;
switch lower(g.crossfphaseunit)
	case {'degree', 'radian'} ;  
	otherwise disp('Error: Crossfphaseunit must be either ''degree'' or ''radian'''); return;
end;
switch lower(g.caption)
	case {'on', 'off'} ;  
	otherwise disp('Error: Caption must be either ''on'' or ''off'''); return;
end;
switch lower(g.polarity)
	case {'pos', 'posneg'} ;  
	otherwise disp('Error: Polarity must be either ''pos'' or ''posneg'''); return;
end;
if ~isempty(g.envvert),
    if ~iscell(g.envvert) & ~( isstruct(g.envvert{1}) | isnumeric(g.envvert{1}) )
        disp('Error: Invalid type for Envvert.'); return;
    end
end
if ~isempty(g.latency) & ~isnumeric(g.latency)
	disp('Error: Latency must be a vector'); return;
end;	
if length(g.scalepower) ~= 2
	disp('Error: Scalepower must be a 2-element array'); return;
end;
if length(g.scalecoher) ~= 2
	disp('Error: Scalecoher must be a 2-element array'); return;
end;
if (length(g.diskscale) ~= 1 | g.diskscale < 0)
        disp('Error: Diskscale must be a scalar value >= 0.'); return;
end
if size(g.colmapcoh,2) ~= 3
	disp('Error: Colmapcoh must be a colormap (3 columns)'); return;
end;
if size(g.colmapcrossf,2) ~= 3
	disp('Error: Colmapcrossf must be a colormap (3 columns)'); return;
end;
if size(g.circfactor,1) ~= size(g.circfactor,2)
	disp('Error: Circfactor must be a square matrix'); return;
end;
if ~iscell(g.coordinates) & ~isempty(g.circfactor)
    if size(g.circfactor,1) ~= size(g.coordinates,1)
        disp('Error: Circfactor must have the same number of rows as the number of rows of coordinates'); return;
    end;
    if nbcomponents ~= size(g.coordinates,1)
        disp('Error: The array of selected components must have length nrows of the array coordinates'); return;
    end;
end;
if ~isstr(g.envylabel)
	disp('Error: envelope label must be a string'); return;
end;	
if ~isempty(g.envelope)
	if (size( g.envelope,1 ) ~=2) | (size( g.envelope,2 ) ~= length(times)) | (size( g.envelope,3 ) ~= nbconditions)
		fprintf('Error: Enveloppe array does not have the right size (%s), i.e. 2 x %d (number of time points) x %d (number of conditions)\n', int2str(size( g.envelope)), length(times), nbconditions); return;
	end;
end;
if ~isempty(g.condtitle)
    if iscell(g.condtitle), g.condtitle = strvcat(g.condtitle{:}); end;
	if size( g.condtitle,1 ) ~= nbconditions
		fprintf('Error: The number of rows in the title array(%d) must match the number of conditions (%d)\n', size(g.condtitle,1), nbconditions); return;
	end;
end;
if length(g.plotorder) ~= length(selected)
    error([ 'Error: ''plotorder'' must be the same size as the number of selected components:' int2str(length(selected)) ]);
end;
if max(g.plotorder) > max(selected)
    error([ 'Error: ''plotorder'' must be below the number of selected components:' int2str(max(selected)) ]);
end;
if ~isempty(g.framefolder)
    [tmp1 tmp2] = mkdir('/', g.framefolder(2:end) );
    if g.framefolder(end) == '/', g.framefolder(end) = []; end;
end;  

% create movie
% ------------
disp('A movie is being saved under output.avi (movie parameters shown below):');
mov = avifile(g.moviename, g.movieopts{:});

% other variables
% ---------------
%limits: power -6 to 6
%limits: ITC 0-1
%limits: coherence 0-1
%limits: coherence angle -180 to 180 
g.factproj = [-71 88 -71];
g.projcolor = [0.35 0.35 0.35];
g.rthistcolor  = [1 1 1];
g.resmult = 1;
currentphase   = zeros( length(selected), length(selected), nbconditions);
tmp = ALLERSP{1,1};
nwin = size(tmp,2);
	
%for index = 1:64
%	circle(1+index,1, 0.5, g.colormaphsv(index, :));
%end;

% optional resqure of all coordinates
% -----------------------------------
g.magnify = g.magnify/4;

% compute RT distribution
% -----------------------
if ~isempty(g.rt)
	RTdist = zeros(nbconditions,nwin);
	for index = 1:nbconditions	
		if ~isempty(g.rt{index})
			timestep = (times(2)-times(1))/2;
			for indeximage = 1:nwin
				RTdist(index, indeximage) = length( intersect( find( g.rt{index} > times(indeximage)-timestep ) , ...
                                                               find(  g.rt{index} <= times(indeximage)+timestep ) ) );
			end;
			RTdist(index,:) = RTdist(index,:)/max(RTdist(index,:));
		end;	
	end;	
	RTdist = RTdist/max(RTdist(:));
end;	

figure( 'position', [100, 100, ceil(nbconditions*g.size(1)/4)*4, ceil(g.size(2)/4)*4], ...
        'PaperPositionMode', 'auto', 'papertype', 'A1', 'visible',g.visible); %'paperorientation', 'landscape' );

axis off
if 	strcmpi(g.framesout, 'ppm')
    r = 0.8465;
    pos = get(gcf,'position');
    if floor(pos(3)/r)> 1280
        fact = 1280/(pos(3)/r);
        set(gcf, 'position', [ 0 0 1280  floor(pos(4)/r*fact) ]);
    else
        set(gcf, 'position', [ 0 0 floor(pos(3)/r), floor(pos(4)/r) ]);
    end;
end;
pos = get(gca,'position');
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)];

% compute selected latency point
% ------------------------------
if ~isempty(g.latency)
	alltimepoints = [];
	for index = 1:length(g.latency)
		[tmp tmptimepoint] = min(abs(g.latency(index)-times));
		alltimepoints = [ alltimepoints tmptimepoint];
	end;	
else 
	if isempty(g.frames)
		alltimepoints = 1:nwin;
	else
		alltimepoints = g.frames;
	end;
end;

% make black patch behind figure
% ------------------------------
hback = axes('position' , [0 0 1 1], 'xtick', [], 'ytick', [], 'box', 'off');
hpatch = patch([0 1 1 0], [0 0 1 1], g.backcolor); xlim([0 1]); ylim([0 1]);
set(hpatch, 'facecolor' , g.backcolor, 'edgecolor', 'none');

% compute flashes latency
% -----------------------
if ~isempty(g.flashes)
	if iscell(g.flashes)
		for index = 1:length(g.flashes)
			flasheslat(index) = g.flashes{index}{1};
			flashescol{index} = g.flashes{index}{2};
		end;
	else
		flasheslat = g.flashes;
		for index = 1:length(g.flashes)
			flashescol{index} = [0.5 0.5 0.5];
		end;
	end;
	allflashes = [];
	for index = 1:length(g.flashes)
		[tmp tmptimepoint] = min(abs(flasheslat(index)-times));
		allflashes = [ allflashes tmptimepoint];
	end;
	%hpatch = patch([ 0.02 .11 .11 0.02], [0.05 0.05 0.925 0.925], [0.5 0.5 0.5]); lateral
	%hpatch = patch([ 0 1 1 0], [0 0 1 1], [0.5 0.5 0.5]); full
	%hpatch = patch([ 0.13 0.84 0.84 0.13 ], [0.92 0.92 1 1], [0.5 0.5 0.5]); %up
    hpatch = patch([ 0.13 0.84 0.84 0.13 ], [0.8 0.8 0.93 0.93], [0.5 0.5 0.5]);
	set(hpatch, 'facecolor', 'w', 'edgecolor', 'none');
	xlim([0 1]); ylim([0 1]);
	posf = 0; % used as a counter to preserve color
end;	

% draw axes and display images
% ----------------------------
ordinate = 0.2;
max_ordinate = 1-1.4*ordinate;   % makes space at top for figure title  
maxcoordx    = 1.1-1/nbconditions/4;
coords = g.coordinates;
g.coordinates = {};
for i=1:nbconditions
    
    % plot 3d head (*0.9 added for Nick - Arno).
    % ------------
	hh(i) = axes('position', [0+maxcoordx/nbconditions*(i-1), ordinate, maxcoordx/nbconditions*0.9, max_ordinate].*s+q );
    gr = [ 0.3 0.3 0.3 ];
    g.dipplotopt = { 'coordformat' g.coordformat 'gui', 'off', 'cornermri', 'on', 'color', { gr gr gr gr gr gr gr gr gr } };
    if iscell(coords)
        for index = 1:size(coords{1}, 1);
            dipstruct(index).posxyz = coords{1}(index,:);
            dipstruct(index).momxyz = [0 0 0];
            dipstruct(index).component = index;
            dipstruct(index).rv = 0.1;
        end;
    else
        for index = 1:size(coords, 1);
            dipstruct(index).posxyz = coords(index,:);
            dipstruct(index).momxyz = [0 0 0];
            dipstruct(index).component = index;
            dipstruct(index).rv = 0.1;
        end;
    end;        
    
    dipplot( dipstruct, 'view', g.view, g.dipplotopt{:}); axis off;
    
    %g.maxc = 100;
    %surface([-2 -2; -2 -2]*g.maxc, [-20 20; -20 20]*g.maxc,[-20 -20; 20 20]*g.maxc, repmat(reshape([0 0 0], 1, 1, 3), [2 2 1]), 'facelighting', 'none');
    %surface([-20 20; -20 20]*g.maxc,[2 2; 2 2]*g.maxc, [-20 -20; 20 20]*g.maxc,     repmat(reshape([0 0 0], 1, 1, 3), [2 2 1]), 'facelighting', 'none');

    %camproj('perspective');
    set(gca, 'cameraviewanglemode', 'manual'); % disable change size
    axis vis3d % same as above (for security)
    camlight left
    camlight right
    view(g.view)
    %camzoom(1.2)
   
    for index = 1:length(dipstruct)
        htmp = findobj(gca, 'tag', [ 'dipole' int2str(index) ]);
        for dipindex = 1:length(htmp)
            tmpstruct = get(htmp(dipindex), 'userdata');
            if isstruct(tmpstruct) % look for dipole location % THIS DOES NOT WORK
                if isfield(tmpstruct, 'pos3d')
                    g.coordinates{i}(index, :) = tmpstruct.pos3d;
                elseif isfield(tmpstruct, 'eleccoord')
                    g.coordinates{i}(index, :) = tmpstruct.eleccoord;
                else
                    tmpstruct
                    error('Field not found in tmpstruct');
                end;
            end;
        end;
        delete(htmp);
    end;
    %h = plot3(g.coordinates{i}(:, 1),  g.coordinates{i}(:, 2),  g.coordinates{i}(:, 3), 'r.', 'markersize', 30); 
    %dsaf
    xltmp = xlim;
    yltmp = ylim;
    g.dimratio = (xltmp(2) - xltmp(1)) / (yltmp(2) - yltmp(1));
    
	axis off;
	if ~isempty(g.condtitle)
		h = title(g.condtitle(i,:));
		if ~isempty(g.condtitleformat)
			set(h, g.condtitleformat{:} );
		end;
	end;	

    % this axis is used for the enveloppe but
    % also used to print current time (which is why it is always created
    e(i) = axes('position', [0.1/nbconditions+maxcoordx/nbconditions*(i-1), 0, ...
                maxcoordx/nbconditions-0.1/nbconditions, ordinate].*s+q,'visible', g.visible);
end;

% draw captions if necessary
% --------------------------
countl = 1;
switch lower(g.caption)
 case 'on' , 
  xlimnorm = (1.1-maxcoordx)/(maxcoordx/nbconditions) * g.xlimaxes;
  ylimnorm = 0.45/(1-ordinate) * g.ylimaxes;
  switch g.power, case 'on',
      c(countl) = axes('position', [maxcoordx, -0.1,    (1.1-maxcoordx), 0.45].*s+q, 'xlim', xlimnorm, ...
                  'ylim', ylimnorm,'visible', g.visible, 'color', 'w' );
      % draw 3 spheres
      [xstmp ystmp zs] = sphere(15);
      l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
      normals = reshape([xstmp./l ystmp./l zs./l],[16 16 3]);
      tmpsize = 0.5; xs1 = tmpsize*ystmp; ys1 = tmpsize*xstmp; zs1 = tmpsize*zs;
      tmpsize = 0.9; xs2 = tmpsize*ystmp; ys2 = tmpsize*xstmp; zs2 = tmpsize*zs + 2;
      tmpsize = 0.1; xs3 = tmpsize*ystmp; ys3 = tmpsize*xstmp; zs3 = tmpsize*zs - 1.5;
      colorarray = repmat(reshape([1 1 1],  1,1,3), [size(zs,1) size(zs,2) 1]);
      handles = surf(xs1, ys1, zs1, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                     'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3); hold on;
      handles = surf(xs2, ys2, zs2, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                     'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
      handles = surf(xs3, ys3, zs3, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                     'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
      axis off;
      camlight left
      camlight right
      view([1 0 0])
      lightangle(45,0);
      lighting phong;
      material shiny;
      axis equal;
      set(gca, 'zlim', [-2 4]);
      text(0, 1.3, 2, [ '+' num2str(g.scalepower(2),2) ], 'fontweight', 'bold');
      text(0, 1, 0, '0', 'fontweight', 'bold');
      text(0, 0.5, -1.5, [ num2str(g.scalepower(1),2) ' dB' ], 'fontweight', 'bold');
      %scalepower(mean(xlimnorm), min(ylimnorm)+0.2, g); % see function at the end
      %axis off;
      %countl = countl + 1;
  end;
  switch g.itc, case 'on',
	  c(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.29 , (1.1-maxcoordx)/2, 0.14].*s+q, ...
				  'visible', g.visible, 'color', 'none' );
      countl = countl + 1;
      if strcmpi(g.polarity, 'posneg') % negative ITCs (difference only) ?
          cbar( [-1 1], [-1 1], g.colmapcoh, 'vert', 'circle', g);
          ylabel('ITC', 'fontweight', 'bold');
          set(gca, 'ytick', [-1 0 1], 'yticklabel', [-1 0 1], 'xticklabel', [], 'box', 'off');
	  else 
          cbar( [0 1], [0 1], g.colmapcoh(length(g.colmapcoh)/2:end,:), 'vert', 'circle', g);
          ylabel('ITC', 'fontweight', 'bold');
          set(gca, 'ytick', [0 1], 'yticklabel', [0 1], 'xticklabel', [], 'box', 'off');
      end; 
      axis off;
      text(-0.2, 0, '0', 'fontsize', 10);
      text(-0.2, 1, num2str(g.scaleitc,2), 'fontsize', 10);
      text(-0.8, 0.55, 'ITC', 'fontsize', 11, 'fontweight', 'bold');
  end;
  switch g.crossf, case 'on',
      c(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.47 , (1.1-maxcoordx)/4, 0.14].*s+q, ...
                  'visible', g.visible, 'color', 'none' );
      countl = countl + 1;
      if strcmpi(g.polarity, 'posneg') % negative ITCs (difference only) ?
          cbar( [-1 1], [-1 1], g.colmapcrossf, 'vert', '', g);
          ylabel('Cross-Coh' , 'fontweight', 'bold');
          set(gca, 'ytick', [-1 0 1], 'yticklabel', [g.scalecoher(1) 0 g.scalecoher(2)], 'xticklabel', []);
      else
          cbar( [0 1], [0 1], g.colmapcrossf(length(g.colmapcrossf)/2:end,:), 'vert', '', g);
          ylabel('Cross-Coh' , 'fontweight', 'bold');
          set(gca, 'ytick', [0 1], 'yticklabel', [g.scalecoher(1) g.scalecoher(2)], 'xticklabel', []);      
      end;
      switch g.crossfphasespeed, case 'on',
          c(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.69,(1.1-maxcoordx)/2, 0.25 ].*s+q, ...
                      'visible', g.visible );
          countl = countl + 1;
          scalecoher([0.02 1], [0.04 0.96], 5, g); % see function at the end
      end;
  end;
 case 'off', maxcoordx = 1;
end;	

% draw white axis on envelop if flashes DOES NOT WORK WHEN PRINTING IN EPS
% -------------------------------------
%if ~isempty(g.flashes)
%	if ~isempty(g.envelope) % draw axis for the envelope
%		eflash = axes('position', [0 0 maxcoordx-0.1 ordinate].*s+q, ...
%					  'xtick', [], 'ytick', [], 'box', 'off', 'visible', g.visible, 'color', 'none'); 
%		hpatch2 = patch([ 0 1 1 0], [0 0 1 1], [0.5 0.5 0.5]); set(hpatch2, 'facecolor', 'w', 'edgecolor', 'none');
%	end;
%end;

% scan time windows
% -----------------
set(gcf, 'renderer', 'zbuffer');
for indeximage = alltimepoints
    
	fprintf('Processing image %d\n', indeximage);

	% invert background if necessary
	% ------------------------------
	if ~isempty(g.flashes)
		%axes(hback); set (gcf, 'visible', g.visible);
		if ~isempty(find(indeximage == allflashes))
			posf = find(indeximage == allflashes);
			set(hpatch, 'facecolor', flashescol{posf});
		elseif posf == 0 % allow the color to stay 2 images
			set(hpatch, 'facecolor', 'w');
		else
			posf = 0;
		end;
	end;
	
	for tmpcond=1:nbconditions
          axes(hh(tmpcond)); set (gcf, 'visible', g.visible);
          
          
          % clean images and update view
          % ----------------------------
          if ~isempty(g.path3d)
              angle = (indeximage-1)/length(alltimepoints)*360;
              camorbit( cos(angle/180*pi)*g.path3d(1), sin(angle/180*pi)*g.path3d(2) );
          end;
          delete( findobj( hh(i), 'tag', 'tmpmov') );
          set (gcf, 'visible', g.visible); 
          if ~isempty(g.title) & i == 1
            t = textsc(g.title,'title');
            set(t,'VerticalAlignment','top', 'fontsize', 15);
          end;  
	
          % draw correlations
          % -----------------  
          switch lower(g.crossf), case 'on', 
              for index1 = selected
                  for index2 = selected
                      if index2 > index1
                          
                          tmpcrossfpow = ALLCROSSF     	 { index1, index2, tmpcond };
                          tmpcrossfang = ALLCROSSFANGLE    { index1, index2, tmpcond };
                          tmppower  = mean(tmpcrossfpow( FREQS, indeximage));
                          tmpangle  = mean(tmpcrossfang( FREQS, indeximage));
                          
                          if strcmp(lower(g.crossfphaseunit), 'radian'), tmpangle = tmpangle/pi*180; end;
                          %fprintf('%d-%d -> power %1.1f\n', index1, index2, tmppower);
                          drawconnections( g.coordinates{tmpcond}( index1,: ), g.coordinates{tmpcond}( index2,: ), ...
                                           tmppower, tmpangle, g.circfactor(index1, index2), g);
                      end;	
                  end;	
              end;
          end;
	
          % draw circles
          % ------------     
          for index1 = g.plotorder(:)'
              tmptimef = ALLERSP{ index1, tmpcond};
              tmppow   = mean(tmptimef( FREQS, indeximage)); % size is power
              tmptimef = ALLITC{ index1, tmpcond};
              tmpitc = mean(abs(tmptimef( FREQS, indeximage))); % color is ITC
                                                           %index1, tmpitc, tmppow,
              drawcircle( g.coordinates{tmpcond}( index1,: ), tmppow, tmpitc, g);
          end;


          % draw a bar for time probability
          % -------------------------------
          if ~isempty(g.rt)
              if ~isempty(g.rt{tmpcond}) 
                  ll = line([g.rthistloc(1)-g.rthistloc(3)/2 g.rthistloc(1)+g.rthistloc(3)/2], [g.rthistloc(2) g.rthistloc(2)]);
                  set(ll, 'linewidth', 2*g.resmult, 'color', 'k'); 
                  barheight = RTdist(tmpcond, indeximage)*g.rthistloc(4);
                  x1 = g.rthistloc(1)-0.65*g.rthistloc(3)/2;
                  x2 = g.rthistloc(1)+0.65*g.rthistloc(3)/2;
                  y1 = g.rthistloc(2);
                  y2 = g.rthistloc(2)-barheight;
                  ll = patch([x1 x1 x2 x2], [y1 y2 y2 y1], g.rthistcolor, 'linewidth', 2*g.resmult);
              end;
          end;
    end;

	% draw the enveloppe of the signal if necessary
	% ---------------------------------------------
    axes(e(tmpcond)); cla; axis off; set (gcf, 'visible', g.visible);
	if ~isempty( g.envelope )
          minordinate = min(min(min(g.envelope)));
          maxordinate = max(max(max(g.envelope)));
          for tmpcond = 1:nbconditions
            axes(e(tmpcond)); cla; axis on; set (gcf, 'visible', g.visible);
            plot(times, g.envelope(:,:,tmpcond), 'k', 'linewidth', 2*g.resmult); hold on;
            set(gca, 'ylim', [minordinate maxordinate]);
            set(gca, 'xlim', [times(1) times(end)]);
            plot([times(indeximage) times(indeximage)], [minordinate maxordinate], 'b', 'linewidth', 2*g.resmult);
            xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', 12*g.resmult); set(gca, 'box', 'off');
            set(gca, 'fontsize', 10*g.resmult);
            if tmpcond == 1
              ylabel(g.envylabel, 'fontweight', 'bold', 'fontsize', 12*g.resmult);
            end;

            % draw vertical lines if needed
            % -----------------------------
            if ~isempty(g.envvert)
              drawvert(g.envvert, tmpcond,  [minordinate maxordinate]);
            end;
          end;
          
          % put the time on the ERP axis
          % ----------------------------
          %coordx1 = (g.xlimaxes(2)-g.xlimaxes(1))*0.1 + g.xlimaxes(1);
          %coordy1 = (g.ylimaxes(2)-g.ylimaxes(1))*0.87 + g.ylimaxes(1);
	end;
    
    % put the time in the left bottom corner
    % --------------------------------------
    tt = text(-0.1, -0.25, sprintf('%d ms', round(times(indeximage))), 'unit', 'normalized');
    set(tt, 'fontsize', 12*g.resmult, 'horizontalalignment', 'right', 'tag', 'tmpmov', 'color', 'w');
    
    % last 3-D settings
    % -----------------
    lighting phong;
    material shiny;
    setfont(gcf, 'color', [0.99 0.99 0.99]); % warning, for some reasons white does not print
    for index = 1:length(c)
        axes(c(index)); % bring back legend to front
    end;
        
	% save the file for a movie
	% -------------------------
    movframes = getframe(gcf);
    mov = addframe(mov,movframes);
    if strcmpi(g.framesout, 'tiff')
        command2 = sprintf('print -dtiff %s/image%4.4d.tiff', g.framefolder, indeximage);
        eval(command2);
    elseif strcmpi(g.framesout, 'eps')
        command2 = sprintf('print -depsc -loose %s/image%4.4d.eps', g.framefolder, indeximage);
        eval(command2);
    elseif 	strcmpi(g.framesout, 'ppm')
        command2 = sprintf('print -dppm -loose %s/image%4.4d.ppm', g.framefolder, indeximage);
        eval(command2);
    else % fig format
        hgsave(sprintf('%s/image%4.4d.fig', g.framefolder, indeximage));
        if strcmp(g.visible, 'on')
            drawnow;
        end;
    end;
end;
mov = close(mov);
return;

% function to draw circles
% ------------------------
function [tmpsize, tmpcolor, handles] = drawcircle( tmpcoord, tmpersp, tmpitc, g);
% tmpcoord         coordinate of the circle
% tmpersp          erps power -> radius
% tmpitc           itc -> color
% g                preference

		switch lower(g.power)
			case 'on',  tmpsize = (tmpersp-g.scalepower(1))/(g.scalepower(2)-g.scalepower(1)); % in between 0 and 1 
			case 'off', tmpsize = 0.5;
		end;	
		tmpsize = 0.05 *  tmpsize * (g.xlimaxes(2)-g.xlimaxes(1))+0.1;
        if isnan(tmpitc), tmpitc = 0; end;
        
		switch lower(g.itc)
         case 'on', 
          indexcolor = length(g.colmapcoh)/2+ceil((tmpitc/g.scaleitc)*length(g.colmapcoh)/2);
          if indexcolor < 1 | indexcolor > length(g.colmapcoh)
              error([ 'ITC ' num2str(tmpitc) 'out of bound, use ''scaleitc'' to increase maxitc' ] );
          end;
          tmpcolor = g.colmapcoh( indexcolor,: );
         case 'off', tmpcolor = g.colmapcoh( length(g.colmapcoh)/2,: );
          %case 'on',  tmpcolor = g.colmapcoh( 64-ceil((tmpitc+0.01)*63),: );
          %case 'off', tmpcolor = g.colmapcoh( 64-ceil((0+0.01)*63),: );
		end;
		if tmpersp == 0
			dashed = 1;
		else
			dashed = 0;
		end;		
        
        tmpsize = g.diskscale*tmpsize*100;
		if tmpsize > 0
            if length(tmpcoord) > 2
                [xstmp ystmp zs] = sphere(15);
                l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
                normals = reshape([xstmp./l ystmp./l zs./l],[16 16 3]);
                xs = tmpcoord(1) + tmpsize*ystmp*g.dimratio;
                ys = tmpcoord(2) + tmpsize*xstmp;
                zs = tmpcoord(3) + tmpsize*zs;
                if tmpitc ~= 0
                     colorarray = repmat(reshape(tmpcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
                else colorarray = repmat(reshape([1 1 1],  1,1,3), [size(zs,1) size(zs,2) 1]);
                end;
                %figure;gca; hold on;
                handles = surf(xs, ys, zs, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                               'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
                %axis off; axis equal; lighting phong; camlight left; rotate3d
                if strcmpi(g.project3d, 'on')
                    colorarray = repmat(reshape(g.projcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
                    surf(xs, ys, g.factproj(3)*ones(size(zs)), colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
                    surf(xs, g.factproj(2)*ones(size(ys)), zs, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
                    surf(g.factproj(1)*ones(size(xs)), ys, zs, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
                end;
            else
                circle( tmpcoord(1), tmpcoord(2), tmpsize, tmpcolor, 'k', 0, 360, dashed, fastif(dashed, 2, 1));
            end;
		end;
return;

% function to draw the lines
% --------------------------
function handles = drawconnections( pos1, pos2, crossfpower, crossfangle, circfact, g);
% pos1, pos2		position of the points
% crossfpower       coherence power for width of the line
% crossfangle       coherence angle for color and speed of the line
% cirfact           curvature of the line
% g                 preference

	% normalize values depending on scaling
	% -------------------------------------
	%g.scalecoher = 2 * g.scalecoher / (g.xlimaxes(2)-g.xlimaxes(1));
	%g.scalepower = 2 * g.scalepower / (g.xlimaxes(2)-g.xlimaxes(1));
        
    % if the two circle are too close and do not draw the line
    % --------------------------------------------------------
    distance = sqrt(sum((pos1-pos2).^2));
    if distance < 0.05*(g.ylimaxes(2) - g.ylimaxes(1)), return;
    end;
    
	crossfpowerabs = abs(crossfpower);
    switch lower(g.crossfcoh)
		case 'on', tmpthick   = (crossfpowerabs-g.scalecoher(1))/(g.scalecoher(2)-g.scalecoher(1));	% determine thickness = coherence amplitude
		case 'off', tmpthick  = 0;
	end;

	sizec = size( g.colmapcrossf,1 );
	switch lower(g.crossfphasecolor)
		case 'on',  tmpcolor  = g.colmapcrossf( sizec/2+ ceil(tmpthick*(sizec/2-1)+1)*sign(crossfpower), : );    % determine color = coherence phase
		case 'off', tmpcolor  = g.colmapcrossf( sizec/2, : );
	end;
	%tmpthick = 30 * (tmpthick-0.1); % does not vary with the axis zoom
	tmpthick = 30 * (tmpthick); % adjusted for Nick
	
	% absolute value to 90 degree determine speed
	switch lower(g.crossfphasespeed)
		case 'on',  curphase = (crossfangle+180)/360; % phase from 1 to 0
		case 'off', curphase = 0.5;
	end;
    if crossfpower == 0, tmpthick = 0; end;
    
	if tmpthick > 0        
        [xc yc zc] = cylinder( g.resmult*tmpthick/300*100, 10);
        colorarray = repmat(reshape(tmpcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
        handles = surf(xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
              'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', 'phong', 'ambientstrength', 0.3);
        [xc yc zc] = adjustcylinder2( handles, [pos1(1) pos1(2) pos1(3)], [pos2(1) pos2(2) pos2(3)] );
        
        % compute cylinder normals (have to bias normal closer to sphere
        % to get a specular point
        cx = mean(xc,2); cx = [(3*cx(1)+cx(2))/4; (cx(1)+3*cx(2))/4];
        cy = mean(yc,2); cy = [(3*cy(1)+cy(2))/4; (cy(1)+3*cy(2))/4];
        cz = mean(zc,2); cz = [(3*cz(1)+cz(2))/4; (cz(1)+3*cz(2))/4];
        tmpx = xc - repmat(cx, [1 11]);
        tmpy = yc - repmat(cy, [1 11]);
        tmpz = zc - repmat(cz, [1 11]);
        l=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);
        normals = reshape([tmpx./l tmpy./l tmpz./l],[2 11 3]);
        set( handles, 'vertexnormals', normals);
        
        %figure
        %axis off; axis equal; lighting phong; camlight left; rotate3d
        if strcmpi(g.project3d, 'on')
            colorarray  = repmat(reshape(g.projcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
            surf(xc, yc, g.factproj(3)*ones(size(zc)), colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
            surf(xc, g.factproj(2)*ones(size(yc)), zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
            surf(g.factproj(1)*ones(size(xc)), yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
        end;
        %if round(tmpthick) == 7, asdf; end;
	end;
return;

% ***************************************************************************************
%                              Caption and tests
% ***************************************************************************************
				
% function to draw circles at all power
% -------------------------------------
function scalepower(posx, posy, g);

	NBCIRCLE = 3;
	coordy = posy;
	powerscale = [ ceil( g.scalepower(1) ) 0 floor( g.scalepower(2) ) ];
	xlim = get(gca, 'xlim');
	ylim = get(gca, 'ylim');
	
	for i=1:NBCIRCLE
		[tmpsize] = drawcircle( [posx coordy], powerscale(i), 0, g);
		if i == 1, tmpsizeori = tmpsize; end;
        
        if i == NBCIRCLE
             tt = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1g dB', powerscale(i)));
		else tt = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1g', powerscale(i)));
        end; 
        set(tt, 'fontsize', 10*g.resmult, 'horizontalalignment', 'left', 'fontweight', 'bold');
		coordy = coordy + tmpsize + 0.2*(ylim(2)-ylim(1));
		
		%command2 = sprintf('print -depsc -loose scale%d.eps', i);
		%eval(command2);
		%cla;
	end;
	set(gca, 'xlim', xlim, 'ylim', ylim-tmpsizeori, 'clipping', 'off', 'fontsize', 10*g.resmult);
return;

% function to draw lines at all coherence
% ---------------------------------------
function scalecoher(posx, posy, thickness,g);
	compter = -5;
	for i=linspace( posy(1), posy(2), 11)
		superline( [ posx(1) posx(2) ], [ i i ], 'b', thickness*g.resmult, mod(compter/10, 1));  
		compter = compter + 1;
	end;	
	%ylabel('Phase-Coh', 'fontweight', 'bold', 'fontsize', 12*g.resmult);
	set(gca, 'box', 'on', 'ylim', [0 1], 'ytick', [0 0.5 1], ...
			 'yticklabel', strvcat('-180�','0�','180�'), 'xlim', [0 1], 'xtick', [], 'xticklabel', [], 'fontsize', 10*g.resmult);
	%hold on; ff = fill([0 0.02 0.02 0], [0 0 1 1], 'w'); set(ff, 'edgecolor', 'w');
	%hold on; ff = fill([0 0 1 1], [0 0.02 0.02 0], 'w'); set(ff, 'edgecolor', 'w');
return;

% colorbar special
% ----------------
function cbar( X, Y, colors, orientation, style, g );
% colors = colors to plot
% orientation = 'vert' or 'horiz'
% style = shape of the colorbar, 'circle' = circle, bar otherwise

	NSEGMENTS = size(colors,1)-1;
	compter = 0;
	switch lower(orientation)
		case 'horiz'
			inc = (X(2)-X(1))/NSEGMENTS;
			for i=linspace(X(1),X(2)-inc,NSEGMENTS);
				compter = compter + 1;
				hold on;
				h = fill( [i i i+inc i+inc], [Y(1) Y(2) Y(2) Y(1)], colors(size(colors,1)+1-compter, :)); 
				set(h, 'edgecolor', 'none');
			end;
		case 'vert'
			inc = (X(2)-X(1))/NSEGMENTS;
			for i=linspace(Y(1),Y(2)-(Y(2)-Y(1))/NSEGMENTS,NSEGMENTS);
				hold on;
				switch style
					case 'circle', 
						mid     = (X(2)-X(1))/2;
						angle   = acos( compter/NSEGMENTS*2-1);
						angle1  = acos( (compter+1)/NSEGMENTS*2-1);
						coordx1 = mid - sin( angle )*mid;
						coordx2 = mid + sin( angle )*mid;
						coordx3 = mid + sin( angle1 )*mid;
						coordx4 = mid - sin( angle1 )*mid;
						coordx = real([coordx1 coordx2 coordx3 coordx4]);
					otherwise,	coordx = [X(1) X(2) X(2) X(1)];
				end;	
                compter = compter + 1;
				h = fill( coordx, [i i i+inc i+inc], colors(compter, :));
				set(h, 'edgecolor', 'none');
			end;
		otherwise
			disp('Orientation has to be ''vert'' or ''horiz''');
	end;
	set(gca, 'fontsize', 10*g.resmult);
	if strcmp(style, 'circle'), axis square; end;
return;

% draw vertical lines
% -------------------
function drawvert(tmpev, tmpcond, coords);
    
    if isstruct(tmpev) | isstruct(tmpev{1})
        
        % cooper envert
        %--------------
        if length(tmpev) > 1,
            verts = tmpev{ tmpcond };
        else   verts = tmpev{1};
        end
        
        for v=verts,
            if isstruct(v), ev = v;
            else,           ev.time = v;  ev.color = 'k'; ev.style = '-';
            end
            
            phandle = plot([ev.time ev.time], coords, ev.style, 'linewidth', 1);
            set(phandle,'color',ev.color);
        end;
    else
        % standard envvert
        % ----------------
        for index = 1:length(tmpev)
            if ~iscell(tmpev{index}), 
                plot([tmpev{index} tmpev{index}], coords, 'k');
            else
                phandle = plot([tmpev{index}{1} tmpev{index}{1}], coords, 'k');
                if length(tmpev{index}) > 2
                    set(phandle,tmpev{index}{2:end});
                end;
            end;
        end;
    end;
    return;
    
% check the flux 
% --------------
for indeximage = 1:nwin-7
	index1 = 1;
	index2 = 2;	 
	% determine color = coherence phase
	tmpcrossf = ALLCROSSFANGLE     { index1, index2, 1 };
	tmpvalue  = mean(tmpcrossf( 1:2, indeximage));
	tmpcolor  = colormaphsv( ceil((tmpvalue+180)/360*63) + 1, : );    % index for color

	% absolute value to 90 degree determine speed
	speed = 1 - abs(90 - abs(tmpvalue))/90; % speed from 1 to 0
	currentphase(index1, index2) = currentphase(index1, index2) + sign(tmpvalue)*speed/3; % 1 cycle in 5 images at max speed

	superline( [ 2 1] , [ 1+indeximage 0.8+indeximage], 5, tmpcolor, mod(currentphase(index1, index2),1));
end;
return;
