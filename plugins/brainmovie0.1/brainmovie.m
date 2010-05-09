%brainmovie() - generate a sequence of images showing event-related coherence,
%               event-related spectral perturbations, and inter-trial coherence
%               of localized EEG waveforms. Uses outputs of timef() and cross().
%Usage:
%   >> brainmovie(ersps,itcs,crossfs_amp,crossfs_phase,times,freqs,selected,...
%                         'keyword1',value1,...); % creates files image0001.eps, etc.
%
%Inputs:
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
%Optional 'keyword' parameters:
% 'latency'   - plot only a subset of latencies. The time point closest to the 
%               latency given are plotted. Default = empty, all latencies.
% 'frames'    - vector of frame indices to compute
% 'resolution'- ['low' or 'high'], 'high' -> multiply the size of the image by 3 
%               for subsequent antialiasing and high quality movie generation 
%               {default: 'low'}
% 'framesout' - ['eps'|'ppm'|'fig'] Default format for saving frames on disk. Default is '.eps'.
% 'rt'        - cell array of vector containing reaction times of the subject in 
%               each conditions (default {} -> ignored)
% 'rthistloc' - location and size of rt histograms in individual axes. 
%               [abscissa ordinate width maxheight].
% 'square'    - ['on'|'off'] re-square all coordinates (so X and Y width is the same)
%               default is 'on';
% 'magnify'   - integer magnification factor for graphics. Default is 1.
% 'size'      - [widthcond height] output image size {default [400,400]}
%               widthcond is the width of a single condition plot (in pixels)
% 'head'      - [FILENAME], plot the head background image using .pcx image in FILENAME
% 'polarity'  - ['pos'|'posneg'] polarity for ITC and crossf. 'pos' = only positive values
%               'posneg' = positive and negative values.
% 'visible'   - ['on'|'off'] show the images on the screen or keep them hidden {default 'on'}
% 'mesh'      - ['on'|'off'] show mesh of radius one.
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
%                                      cross-coherence phase {def: on}
% 'crossfphaseunit'  - ['degree'|'radian']. Coherence phase angle unit {Default is degree}.
% 'colmapcrossf' - colormap array for arcs {default: hsv(64) with green as 0} 
% 'colmapcoh'   - colormap array for disks (according to inter-trial coherence) 
%                      {default: hot(64)}
% 'scalepower'  - [min max] dB range for power (and disk size) variation {default: [-5 5]}  
% 'scalecoher'  - [min max] coherence range {default: [0 1]}
% 'diskscale'   - numeric value that scales the size of disks {default: [1.0]}   
%
% Movie coordinates and axis options
% 'xlimaxes'    - x-axis limits axis for the component locations {default: [-1 1]}
% 'ylimaxes'    - y-axis limits axis for the component locations {default: [-1 to 1]}
% 'coordinates' - 2-column array of [x y] coordinates of the selected components 
%                 {default: spaced evenly around the head circle boundary}  
% 'circfactor'  - (ncomps,ncomps) array of arc curvatures (0=straight; 1=half-round, 
%                 positive or negative values give the sense of rotation) {def: 0s}
% 'envelope'    - (2,points,conditions) envelopes of the average data (ERP) in each condition
%                 (envelope =  min and max traces of each ERP across all channels and times)
% 'envylabel'   - ordinate label for envelope. {Default 'Potential \muV'}
% 'envvert'     - cell array of time indices at which to draw vertical lines.
%                 Can also be a cell array of cell to specify line aspect. For instance
%                 { { 0 'color' 'b' 'linewidth' 2 } {1000 'color' 'r' }} would draw two
%                 lines, one blue thick line at latency 0 and one thin red line at latency 1000.
% 'flashes'     - vector of time indices at which the background flashes.  Specify the color 
%                 of the flash with a cell array of [1,2] cell arrays. 
%                 Ex. { { 200 'y' } { 1500 '5' }} will generate two flashes, 
%                 yellow at 200 ms and red at 1500 ms 
% 'title'       - (string) main movie title
% 'condtitle'   - (string array) condition titles (one condition title per row)
% 'condtitleformat' - list of title properties. Ex: { 'fontize', 12, 'fontweight', 'bold' }
% 'plotorder'   - [integer vector] component plot order from 1 to the number of selected 
%                 components. 
%
%Outputs to disk:
% imageX      - brainmovie() saves a sequence of image files to disk (image0001.eps, ...)
%
%Example:
%
% % Given ICA activations in array icaact (size ncomps,nframes,ntrials), animate (here) 
% % activity at/between two components at 176 points per epoch (from -100 ms to 600 ms 
% % re stimulus onset) assuming a 250-Hz sampling rate and 100 output frames
%
% >> [ersps{1,1},itcs{1,1},powbase,times,freqs] = ...                          % timef for
%                timef(icaact(1,:),176,[-100 600],'Component 1',250,1,32,100); %     1st comp
% >> [ersps{2,1},itcs{2,1},powbase,times,freqs] = ...                          % timef for
%                timef(icaact(2,:),176,[-100 600],'Component 2',250,1,32,100); %     2nd comp
% >> [crossfs_amp{1,2},mcoh,times,freqs,cohboot,crossfs_phase{1,2}] = ...      % crossf for
%      crossf_(icaact(1,:),icaact(2,:),176,[-100 600],'Crossf 1 and 2',250,1,32,100); % both
%
% >> brainmovie( ersps, itcs, crossfs_amp, crossfs_phase, times, [1:2] ); % makes files in pwd
%                                                          image0001.eps, ... image0100.eps
%
% >> !/usr/local/bin/convert images*.eps movie.mpg % Now use ImageMagic 'convert' 
%                                                  % to generate the movie.
%
% Note: Better resolution movies can be generated by .eps -> .ppm -> .avi, 
%       (or, under a planned upgrade to brainmovie, from Matlab6 to .avi directly).
   
% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or
% modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function alltimepoints = brainmovie(ALLERSP,ALLITC,ALLCROSSF,ALLCROSSFANGLE,times,FREQS,selected,varargin);
    
if nargin < 6
	help brainmovie;
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
try, g.head; 			catch, g.head=''; end;
try, g.visible; 		catch, g.visible='on'; end;
try, g.square; 		    catch, g.square='on'; end;
try, g.resolution; 		catch, g.resolution='low'; end;
try, g.rt; 	        	catch, g.rt={}; end;
try, g.power; 	    	catch, g.power='on'; end;
try, g.latency; 	   	catch, g.latency=[]; end;
try, g.itc; 	    	catch, g.itc='on'; end;
try, g.magnify; 	    catch, g.magnify=1; end;
try, g.crossf; 			catch, g.crossf='on'; end;
try, g.crossfcoh; 		catch, g.crossfcoh='on'; end;
try, g.size; 			catch, g.size=[400 400]; end;
try, g.crossfphasecolor;catch, g.crossfphasecolor='on'; end;
try, g.crossfphasespeed;catch, g.crossfphasespeed='on'; end;
try, g.crossfphaseunit; catch, g.crossfphaseunit='degree'; end;
try, g.scalepower;      catch, g.scalepower = [-5 5]; end;
try, g.scalecoher;      catch, g.scalecoher = [0 1]; end;
try, g.diskscale;  catch, g.diskscale = 1; end;
try, g.envelope;        catch, g.envelope = []; end; 
try, g.caption;			catch, g.caption = 'on'; end; 
try, g.frames;			catch, g.frames = []; end; 
try, g.envvert;			catch, g.envvert = {}; end; 
try, g.flashes;			catch, g.flashes = []; end; 
try, g.mesh;			catch, g.mesh = 'off'; end; 
try, g.polarity;		catch, g.polarity = 'pos'; end; 
try, g.framesout;	    catch, g.framesout = 'eps'; end; 
try, g.condtitle;		catch, g.condtitle = []; end; 
try, g.condtitleformat;	catch, g.condtitleformat = {'fontsize', 14', 'fontweight', 'bold' }; end;
try, g.title;			catch, g.title = []; end; 
try, g.envylabel;		catch, g.envylabel = 'Potential \muV'; end; 
try, g.plotorder;       catch, g.plotorder = selected; end;
try, g.colmapcoh;       catch, 
    colormtmp = hot(64);
    colormtmp(end,3) = (colormtmp(end,3)+colormtmp(end-1,3))/2; % white does not come out when the
    g.colmapcoh = colormtmp;                                    % the figure is printed to ppm
    g.colmapcoh(:,1) =  colormtmp(:,2);
    g.colmapcoh(:,2) =  colormtmp(:,3);
    g.colmapcoh(:,3) =  colormtmp(:,1);
    g.colmapcoh = [ g.colmapcoh; colormtmp(end:-1:1,:)];
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
    	else	g.coordinates( index,:) = [ 0.01 0.01];
		end;
		count = count + 1;
    end;
end;
try, g.circfactor; catch, g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end;
if isempty(g.circfactor), g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end;
    
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
	try, ALLERSP{selected,1}; 
    catch, disp('Error: unable to access the defined components in ERSPs (out of bounds)'); return; end;
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
switch lower(g.resolution)
	case {'low', 'high'} ;  
	otherwise disp('Error: Resolution must be either ''low'' or ''high'''); return;
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
switch lower(g.mesh)
	case {'on', 'off'} ;  
	otherwise disp('Error: Mesh must be either ''on'' or ''off'''); return;
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
switch lower(g.framesout)
	case {'eps', 'fig', 'ppm'} ;  
	otherwise disp('Error: Framesout must be either ''eps'', ''ppm'' or ''fig'''); return;
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
if size(g.circfactor,1) ~= size(g.coordinates,1)
	disp('Error: Circfactor must have the same number of rows as the number of rows of coordinates'); return;
end;
if nbcomponents ~= size(g.coordinates,1)
	disp('Error: The array of selected components must have length nrows of the array coordinates'); return;
end;
if ~isstr(g.envylabel)
	disp('Error: envelope label must be a string'); return;
end;	
if ~isempty(g.envelope)
	if (size( g.envelope,1 ) ~=2) | (size( g.envelope,2 ) ~= length(times))
		fprintf('Error: Enveloppe array does not have the right size (%s), instead of (%s) i.e. 2 x %d (number of time points) x %d (number of conditions)\n', int2str([2 length(times) nbconditions]), int2str(size( g.envelope)), length(times), nbconditions); return;
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

% other variables
% ---------------
%limits: power -6 to 6
%limits: ITC 0-1
%limits: coherence 0-1
%limits: coherence angle -180 to 180 

g.rthistcolor  = [1 1 1];
switch lower(g.resolution)
 case 'low', g.resmult = 1;
 case 'high', g.resmult = 3;
end;
currentphase   = zeros( length(selected), length(selected), nbconditions);
tmp = ALLERSP{1,1};
nwin = size(tmp,2);
	
%for index = 1:64
%	circle(1+index,1, 0.5, g.colormaphsv(index, :));
%end;

% optional resqure of all coordinates
% -----------------------------------
g.magnify = g.magnify/4;
g.dimratio = (g.xlimaxes(2) - g.xlimaxes(1)) / (g.ylimaxes(2) - g.ylimaxes(1));
if strcmp(lower(g.square), 'on') 
    disp('Square option disabled');
%	for index = selected
%    	if length(selected) > 1
%			g.coordinates( index,1) = (g.coordinates( index,1) - g.xlimaxes(1))/(g.xlimaxes(2)-g.xlimaxes(1))/g.magnify;
%			g.coordinates( index,2) = (g.coordinates( index,2) - g.ylimaxes(1))/(g.ylimaxes(2)-g.ylimaxes(1))/g.magnify;
%		end;
%	end;
%	g.rthistloc(1) = (g.rthistloc(1) - g.xlimaxes(1))/(g.xlimaxes(2)-g.xlimaxes(1))/g.magnify;
%	g.rthistloc(2) = (g.rthistloc(2) - g.ylimaxes(1))/(g.ylimaxes(2)-g.ylimaxes(1))/g.magnify;
%	g.rthistloc(3) = g.rthistloc(3)/(g.xlimaxes(2)-g.xlimaxes(1))/g.magnify;
%	g.rthistloc(4) = g.rthistloc(4)/(g.ylimaxes(2)-g.ylimaxes(1))/g.magnify;
%	g.xlimaxes = [0 1]/g.magnify;
%	g.ylimaxes = [0 1]/g.magnify;
end;

% compute RT distribution
% -----------------------
if ~isempty(g.rt)
	RTdist = zeros(nbconditions,nwin);
	for index = 1:nbconditions	
		if ~isempty(g.rt{index})
			timestep = (times(2)-times(1))/2;
			for indeximage = 1:nwin
				RTdist(index, indeximage) = length( intersect( find( g.rt{index} > times(indeximage)-timestep ) , find(  g.rt{index} <= times(indeximage)+timestep ) ) );
			end;
			RTdist(index,:) = RTdist(index,:)/max(RTdist(index,:));
		end;	
	end;	
	RTdist = RTdist/max(RTdist(:));
end;	

% create image
% ------------
switch lower(g.resolution)
	case 'high', figure( 'position', [100, 100, nbconditions*g.size(1)*3, g.size(2)*3], 'PaperPositionMode', 'auto', 'papertype', 'A1', 'visible',g.visible); %'paperorientation', 'landscape' );
 otherwise    figure( 'position', ...
                      [100, 100, ceil(nbconditions*g.size(1)/4)*4, ceil(g.size(2)/4)*4], ...
                      'PaperPositionMode', 'auto', 'papertype', 'A1', 'visible',g.visible); %'paperorientation', 'landscape' );
end;

axis off
if strcmpi(g.framesout, 'ppm')
    r = 0.8465;
    pos = get(gcf,'position');
    set(gcf, 'position', [ 0 0 floor(pos(3)/r), floor(pos(4)/r) ]);
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
	hback = axes('position', [0 0 1 1], 'xtick', [], 'ytick', [], 'box', 'off'); set (gcf, 'visible', g.visible);
	%hpatch = patch([ 0.02 .11 .11 0.02], [0.05 0.05 0.925 0.925], [0.5 0.5 0.5]); lateral
	%hpatch = patch([ 0 1 1 0], [0 0 1 1], [0.5 0.5 0.5]); full
	%hpatch = patch([ 0.13 0.84 0.84 0.13 ], [0.92 0.92 1 1], [0.5 0.5 0.5]); %up
        hpatch = patch([ 0.13 0.84 0.84 0.13 ], [0.8 0.8 0.93 0.93], [0.5 0.5 0.5]);
	set(hpatch, 'facecolor', 'w', 'edgecolor', 'none');
	set(hback, 'xlim', [0 1], 'ylim', [0 1]);
	posf = 0; % used as a counter to preserve color
end;	

% draw captions if necessary
% --------------------------
ordinate = 0.2;
switch lower(g.caption)
 case 'on' , 
  maxcoordx = 1-1/nbconditions/4;
  xlimnorm = (1-maxcoordx)/(maxcoordx/nbconditions) * g.xlimaxes;
  ylimnorm = 0.45/(1-ordinate) * g.ylimaxes;
  switch g.power, case 'on',
      c(1) = axes('position', [maxcoordx, -0.1,    (1-maxcoordx), 0.45].*s+q, 'xlim', xlimnorm, ...
                  'ylim', ylimnorm,'visible', g.visible );
      scalepower(mean(xlimnorm), min(ylimnorm)+0.2, g); % see function at the end
      axis off;
  end;
  switch g.itc, case 'on',
	  c(2) = axes('position', [maxcoordx+(1-maxcoordx)/2, 0.29 , (1-maxcoordx)/2, 0.14].*s+q, ...
				  'visible', g.visible );
      if strcmpi(g.polarity, 'posneg') % negative ITCs (difference only) ?
          cbar( [-1 1], [-1 1], g.colmapcoh, 'vert', 'circle', g);
          ylabel('ITC', 'fontweight', 'bold');
          set(gca, 'ytick', [-1 0 1], 'yticklabel', [-1 0 1], 'xticklabel', []);
	  else 
          cbar( [0 1], [0 1], g.colmapcoh(length(g.colmapcoh)/2:end,:), 'vert', 'circle', g);
          ylabel('ITC', 'fontweight', 'bold');
          set(gca, 'ytick', [0 1], 'yticklabel', [0 1], 'xticklabel', []);
      end; 
  end;
  switch g.crossf, case 'on',
      c(3) = axes('position', [maxcoordx+(1-maxcoordx)/2, 0.47 , (1-maxcoordx)/4, 0.14].*s+q, ...
                  'visible', g.visible );
      if strcmpi(g.polarity, 'posneg') % negative ITCs (difference only) ?
          cbar( [-1 1], [-1 1], g.colmapcrossf, 'vert', '', g);
          ylabel('Cross-Coh' , 'fontweight', 'bold');
          set(gca, 'ytick', [-1 0 1], 'yticklabel', [-1 0 1], 'xticklabel', []);
      else
          cbar( [0 1], [0 1], g.colmapcrossf(length(g.colmapcrossf)/2:end,:), 'vert', '', g);
          ylabel('Cross-Coh' , 'fontweight', 'bold');
          set(gca, 'ytick', [0 1], 'yticklabel', [0 1], 'xticklabel', []);      
      end;
      switch g.crossfphasespeed, case 'on',
          c(4) = axes('position', [maxcoordx+(1-maxcoordx)/2, 0.69,(1-maxcoordx)/2, 0.25 ].*s+q, ...
                      'visible', g.visible );
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

% draw axes and display images
% ----------------------------
max_ordinate = 1-1.4*ordinate;   % makes space at top for figure title  

for i=1:nbconditions
	h(i) = axes('position', [0+maxcoordx/nbconditions*(i-1), ordinate, maxcoordx/nbconditions, max_ordinate].*s+q );
	if ~isempty(g.head)
		try, img = imread(g.head, 'pcx'); catch, pwd, g.head, disp('Error: unable to load PCX image file'); return; end;
		imagesc(img); colormap(gray);
	end;
	axis off;
	if ~isempty(g.condtitle)
		xlim = get(gca, 'xlim');
		ylim = get(gca, 'ylim');
		%h = text( (xlim(2)-xlim(1))*0.05+xlim(1), (ylim(2)-ylim(1))*0.05+ylim(1), g.condtitle(i,:));
		h = title(g.condtitle(i,:));
		if ~isempty(g.condtitleformat)
			set(h, g.condtitleformat{:} );
		end;
        axis image;
	end;	
	hh(i) = axes('position', [0+maxcoordx/nbconditions*(i-1), ordinate, maxcoordx/nbconditions, max_ordinate].*s+q, ...
				 'xlim', g.xlimaxes, 'ylim', g.ylimaxes, 'color', 'none', 'ydir', 'reverse', 'visible', g.visible);
	axis off;
	if ~isempty(g.envelope) % draw axis for the envelope
		e(i) = axes('position', [0.1/nbconditions+maxcoordx/nbconditions*(i-1), 0, ...
					maxcoordx/nbconditions-0.1/nbconditions, ordinate].*s+q,'visible', g.visible);
	end;
end;

%anim = imread('animal.pcx');
%dist = imread('distractor.pcx');
%upmouse = imread('mouseup.pcx');
%downmouse = imread('mousedown.pcx');
%hhimg1 = axes('position', [0, 0.7, 0.2, 0.3].*s+q, 'visible', g.visible, 'color', 'none');	 
%hhimg2 = axes('position', [0.5, 0.7, 0.2, 0.3].*s+q, 'visible', g.visible, 'color', 'none'); 
%hhmouse = axes('position', [0.3, 0.6, 0.2, 0.4].*s+q, 'visible', g.visible, 'color', 'none'); 

% scan time windows
% -----------------
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
	
	% clean images
	% ------------
	for i=1:nbconditions
		axes(hh(i)); 	cla; set (gcf, 'visible', g.visible); 
        if strcmpi(g.mesh, 'on')
            hold on; h = circle(0,0,1); hold on;
            set(h, 'color', 'w');
		end;
        if ~isempty(g.title) & i == 1
			%x = (g.xlimaxes(2)-g.xlimaxes(1))*0.2 + g.xlimaxes(1);
			%y = (g.ylimaxes(2)-g.ylimaxes(1))*(-0.06) + g.ylimaxes(1);
			%text(x, y, g.title, 'fontsize', 14, 'fontweight', 'bold' );
            t = textsc(g.title,'title');
            set(t,'VerticalAlignment','top', 'fontsize', 15);
		end;  
	end;

	% draw the images if necessary
	% ----------------------------
	%if abs(times(indeximage)) < 1
	%	axes(hhimg1); 	cla; set (gcf, 'visible', g.visible); imagesc(anim); axis off;
	%	axes(hhimg2); 	cla; set (gcf, 'visible', g.visible); imagesc(dist); axis off;
	%else
	%	axes(hhimg1); 	cla; set (gcf, 'visible', g.visible, 'color', 'none'); axis off;
	%	axes(hhimg2); 	cla; set (gcf, 'visible', g.visible, 'color', 'none'); axis off;
	%end;
	%if abs(RTdist(indeximage) > 0)
	%	axes(hhmouse); 	cla; set (gcf, 'visible', g.visible); imagesc(upmouse); axis off;
	%else	
	%	axes(hhmouse); 	cla; set (gcf, 'visible', g.visible); imagesc(downmouse); axis off;
	%end;
	
	% draw correlations
	% -----------------  
	switch lower(g.crossf), case 'on', 
		for index1 = selected
			for index2 = selected
				if index2 > index1
					for tmpcond = 1:nbconditions
						axes(hh(tmpcond)); set (gcf, 'visible', g.visible);
					
						tmpcrossfpow = ALLCROSSF     	 { index1, index2, tmpcond };
						tmpcrossfang = ALLCROSSFANGLE    { index1, index2, tmpcond };
						tmppower  = mean(tmpcrossfpow( FREQS, indeximage));
						tmpangle  = mean(tmpcrossfang( FREQS, indeximage));
						
						if strcmp(lower(g.crossfphaseunit), 'radian'), tmpangle = tmpangle/pi*180; end;
                        %fprintf('%d-%d -> power %1.1f\n', index1, index2, tmppower);
						drawconnections( g.coordinates( index1,: ), g.coordinates( index2,: ), ...
							tmppower, tmpangle, g.circfactor(index1, index2), g);
					end;	
				end;	
			end;
		end;	
	end;
	
	%axes(hh1); 	cla; set (gcf, 'visible', g.visible);
	%axes(hh2); 	cla; set (gcf, 'visible', g.visible);

	% draw circles
	% ------------     
	for index1 = g.plotorder(:)'
		for tmpcond = 1:nbconditions
			axes(hh(tmpcond)); set (gcf, 'visible', g.visible);

			tmptimef = ALLERSP{ index1, tmpcond};
			tmppow   = mean(tmptimef( FREQS, indeximage)); % size is power
			tmptimef = ALLITC{ index1, tmpcond};
			tmpitc = mean(tmptimef( FREQS, indeximage)); % color is ITC
			%index1, tmpitc, tmppow, 
            drawcircle( g.coordinates( index1,: ), tmppow, tmpitc, g);
		end;
	end;
	
	% put the time
	% ------------ 
	coordx1 = (g.xlimaxes(2)-g.xlimaxes(1))*0.1 + g.xlimaxes(1);
	coordy1 = (g.ylimaxes(2)-g.ylimaxes(1))*0.87 + g.ylimaxes(1);
	tt = text(coordx1 ,coordy1, sprintf('%d ms', round(times(indeximage))) );
	set(tt, 'fontsize', 12*g.resmult, 'horizontalalignment', 'right');
		
	% draw a bar for time probability
	% -------------------------------
	for tmpcond = 1:nbconditions
		if ~isempty(g.rt)
			if ~isempty(g.rt{tmpcond}) 
				axes(hh(tmpcond)); set (gcf, 'visible', g.visible);      
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
	if ~isempty( g.envelope )
		minordinate = min(min(min(g.envelope)));
		maxordinate = max(max(max(g.envelope)));
		for tmpcond = 1:nbconditions
			axes(e(tmpcond)); cla; set (gcf, 'visible', g.visible);
            plot(times, g.envelope(:,:,tmpcond), 'k', 'linewidth', 2*g.resmult); hold on;
            set(gca, 'ylim', [minordinate maxordinate]);
			set(gca, 'xlim', [times(1) times(end)]);
			plot([times(indeximage) times(indeximage)], [minordinate maxordinate], 'b', 'linewidth', 2*g.resmult);
			xlabel('time (ms)', 'fontweight', 'bold', 'fontsize', 12*g.resmult); set(gca, 'box', 'off');
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
	end;		   
			
	% save the file for a movie
	% -------------------------
    if strcmpi(g.framesout, 'eps')
        command2 = sprintf('print -depsc -loose image%4.4d.eps', indeximage);
        eval(command2);
    elseif 	strcmpi(g.framesout, 'ppm')
        command2 = sprintf('print -dppm -loose image%4.4d.ppm', indeximage);
        eval(command2);
    else % fig format
        hgsave(sprintf('image%4.4d.fig', indeximage));
        if strcmp(g.visible, 'on')
            drawnow;
        end;
    end;
end;		 
return;

% function to draw circles
% ------------------------
function [tmpsize, tmpcolor] = drawcircle( tmpcoord, tmpersp, tmpitc, g);
% tmpcoord         coordinate of the circle
% tmpersp          erps power -> radius
% tmpitc           itc -> color
% g                preference

		switch lower(g.power)
			case 'on',  tmpsize = (tmpersp-g.scalepower(1))/(g.scalepower(2)-g.scalepower(1)); % in between 0 and 1 
			case 'off', tmpsize = 0.5;
		end;	
		tmpsize = 0.05 *  tmpsize * (g.xlimaxes(2)-g.xlimaxes(1))+0.1;

        try
            switch lower(g.itc)
             case 'on',  tmpcolor = g.colmapcoh( length(g.colmapcoh)/2+ceil((tmpitc+0.01)*length(g.colmapcoh)/2),: );
             case 'off', tmpcolor = g.colmapcoh( length(g.colmapcoh)/2,: );
              %case 'on',  tmpcolor = g.colmapcoh( 64-ceil((tmpitc+0.01)*63),: );
              %case 'off', tmpcolor = g.colmapcoh( 64-ceil((0+0.01)*63),: );
            end;
        catch,  tmpcolor = g.colmapcoh( length(g.colmapcoh)/2,: ); end;
		if tmpersp == 0
			dashed = 1;
		else
			dashed = 0;
		end;		
		
                tmpsize = g.diskscale*tmpsize;
		if tmpsize > 0
			circle( tmpcoord(1), tmpcoord(2), [tmpsize tmpsize*g.dimratio], tmpcolor, 'k', 0, 360, dashed, fastif(dashed, 2, 1));
		end;
return;

% function to draw the lines
% --------------------------
function newphase = drawconnections( pos1, pos2, crossfpower, crossfangle, circfact, g);
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
    distance = abs(pos1(1)+j*pos1(2)-pos2(1)-i*pos2(2));
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
	tmpthick = 30 * (tmpthick-0.1); % does not vary with the axis zoom
	
	% absolute value to 90 degree determine speed
	switch lower(g.crossfphasespeed)
		case 'on',  curphase = (crossfangle+180)/360; % phase from 1 to 0
		case 'off', curphase = 0.5;
	end;
	
	if tmpthick > 0	
         %fprintf('(%d,%d)->(%d,%d) %3.2f: %3.2f %3.2f %3.2f\n', pos1(1), pos2(1), pos1(2), ...
         %                        pos2(2), distance, crossfpower, crossfangle, circfact);
        if circfact ~= 0
			circpatch( [ pos1(1) pos2(1) ] , [ pos1(2) pos2(2) ], circfact, tmpcolor, g.resmult*tmpthick, 100, mod(curphase,1), 0);
		else
			superline( [ pos1(1) pos2(1) ] , [ pos1(2) pos2(2) ], tmpcolor, g.resmult*tmpthick, mod(curphase,1), 0);
		end;
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

		tt = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1fdB', powerscale(i)));
		set(tt, 'fontsize', 10*g.resmult, 'horizontalalignment', 'right', 'fontweight', 'bold');
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
			 'yticklabel', strvcat('-180º','0º','180º'), 'xlim', [0 1], 'xtick', [], 'xticklabel', [], 'fontsize', 10*g.resmult);
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
				compter = compter + 1;
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
