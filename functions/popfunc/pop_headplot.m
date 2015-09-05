% pop_headplot() - plot one or more spherically-splined EEG field maps 
%                  using a semi-realistic 3-D head model. Requires a
%                  spline file, which is created first if not found.
%                  This may take some time, but does not need to be 
%                  done again for this channel locations montage. A wait 
%                  bar will pop up to indicate how much time remains.
% Usage:
%   To open input GUI:
%   >> EEGOUT = pop_headplot( EEG, typeplot)
%   To run as a script without further GUI input:
%   >> EEGOUT = pop_headplot( EEG, typeplot, ...
%                    latencies/components, title, rowscols, 'key', 'val' ...);
% Required Inputs:
%   EEG        - EEG dataset structure
%   typeplot   - 1=channel, 0=component {Default: 1}
%
% Required Inputs to bypass input GUI
%   latencies/components  - If channels, array of epoch mean latencies (in ms),
%                Else, for components, array of component indices to plot.
%
% Optional inputs:
%   title      - Plot title
%   rowscols   - Vector of the form [m,n] where m is total vertical tiles and n 
%                horizontal tiles per page. If the number of maps exceeds m*n,
%                multiple figures will be produced {def|0 -> 1 near-square page}
%   
% Optional 'Key' 'Value' Paired Inputs
%   'setup'    - ['name_of_file_to_save.spl'] Make the headplot spline file
%   'load'     - ['name_of_file_to_load.spl'] Load the headplot spline file
%   'colorbar' - ['on' or 'off'] Switch to turn colorbar on or off. {Default: 'on'}
%   others...  - All other key-val calls are passed directly to headplot. 
%                See >> help headplot
%
% Output:
%   EEGOUT - EEG dataset, possibly with a new or modified splinefile. 
%
% Note:
%   A new figure is created only when the pop_up window is called or when
%   several channels/components are plotted. Therefore you may call this 
%   command to draw single 3-D topographic maps in an existing figure.
%
%   Headplot spline file is a matlab .mat file with the extension .spl.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 March 2002
%
% See also: headplot(), eegplot(), traditional()

% Copyright (C) 20 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = pop_headplot( EEG, typeplot, arg2, topotitle, rowcols, varargin);

com = '';
if nargin < 1
   help pop_headplot;
   return;
end;   

if isempty(EEG.chanlocs)
    error('Pop_headplot: this dataset does not contain channel locations. Use menu item: Edit > Dataset info');
end;

if nargin < 3 % Open GUI input window
    % remove old spline file
    % ----------------------
    if isfield(EEG, 'splinefile') 
        if ~isempty(EEG.splinefile) && exist(EEG.splinefile, 'file')
            splfile = dir(EEG.splinefile);
            byteperelec = splfile.bytes/EEG.nbchan;
            if byteperelec/EEG.nbchan < 625, % old head plot file
                EEG.splinefile = [];
                disp('Warning: Wrong montage or old-version spline file version detected and removed; new spline file required');
            end;
        end;
    end;
    
    % show the file be recomputed
    % ---------------------------
    compute_file = 0;
    if typeplot == 1 % ********** data plot
        fieldname    = 'splinefile';        
        if isempty(EEG.splinefile) && exist(EEG.splinefile, 'file')          
            if length(EEG.icachansind) == EEG.nbchan & ~isempty(EEG.icasplinefile)
                EEG.splinefile = EEG.icasplinefile;
            else
                compute_file = 1;
            end;
        else
            compute_file = 1;
        end;
    else % ************* Component plot       
        fieldname    = 'icasplinefile';
        if isempty(EEG.icasplinefile) && exist(EEG.icasplinefile, 'file')
            if length(EEG.icachansind) == EEG.nbchan & ~isempty(EEG.splinefile)
                EEG.icasplinefile = EEG.splinefile;
            else
                compute_file = 1;
            end;
        else
            compute_file = 1;
        end;
    end;
            
    if compute_file
        
		 warndlg2( strvcat('headplot() must generate a spline file the first', ...
		                 'time it is called or after changes in the channel location file.', ...
		                 'You must also co-register your channel locations with the', ...
                         'head template. Using a standard 10-20 system montage, default', ...
                         'parameters should allow creating the correct spline file.'), 'Headplot() warning');
    else
    	pop_options      = {};
    end;
    
 	% graphic interface
	% -----------------
    template(1).keywords  = { 'standard-10-5-cap385' };
    template(1).transform = [ -0.355789     -6.33688      12.3705    0.0533239    0.0187461     -1.55264      1.06367     0.987721     0.932694 ];
    %template(1).transform = [ -0.31937  -5.96928 13.1812 0.0509311 0.0172127 -1.55007  1.08221  1.00037  0.923518 ];
    template(2).keywords  = { 'standard_1005' };
    template(2).transform = [ -1.13598      7.75226      11.4527   -0.0271167    0.0155306     -1.54547     0.912338     0.931611     0.806978 ];
%        -0.732155 7.58141 11.8939 -0.0249659 0.0148571 0.0227427 0.932423 0.918943 0.793166 ];
    template(3).keywords  = { 'gsn' 'sfp' };
    %template(3).transform = [ 0 -9 -9 -0.12 0 -1.6 9.7 10.7 11.5 ];
    template(3).transform = [ 0.664455     -3.39403     -14.2521  -0.00241453     0.015519     -1.55584           11      10.1455           12];
    template(4).keywords  = { 'egi' 'elp' };
    template(4).transform =  [ 0.0773 -5.3235 -14.72 -0.1187 -0.0023 -1.5940 92.4 92.5 110.9 ];
    
    transform = [];
    if isfield(EEG.chaninfo, 'filename')
        [tmp transform] = lookupchantemplate(lower(EEG.chaninfo.filename), template);
    end;
            
	if typeplot
		txt = sprintf('Making headplots for these latencies (from %d to %d ms):', round(EEG.xmin*1000), round(EEG.xmax*1000));
	else
		%txt = ['Component numbers (negate index to invert component polarity):' 10 '(NaN -> empty subplot)(Ex: -1 NaN 3)'];
		txt = ['Component numbers to plot (negative numbers invert comp. polarities):' ];
	end;	
    if compute_file
        enableload = 'off';
        enablecomp = 'on';
    else
        enableload = 'on';
        enablecomp = 'off';
    end;
    cb_load = [ 'set(findobj(gcbf, ''tag'', ''load''), ''enable'', ''on'');' ...
                'set(findobj(gcbf, ''tag'', ''comp''), ''enable'', ''off'');' ...
                'set(findobj(gcbf, ''tag'', ''compcb''), ''value'', 0);' ];
    cb_comp = [ 'set(findobj(gcbf, ''tag'', ''load''), ''enable'', ''off'');' ...
                'set(findobj(gcbf, ''tag'', ''comp''), ''enable'', ''on'');' ...
                'set(findobj(gcbf, ''tag'', ''loadcb''), ''value'', 0);' ];  
    cb_browseload = [ '[filename, filepath] = uigetfile(''*.spl'', ''Select a spline file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''load''), ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_browsecomp = [ '[filename, filepath] = uiputfile(''*.spl'', ''Select a spline file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''coregfile''), ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_browsemesh = [ '[filename, filepath] = uigetfile(''*.spl'', ''Select a spline file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''meshfile''), ''style'', ''edit'', ''callback'', '''', ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_browsemeshchan = [ '[filename, filepath] = uigetfile(''*.spl'', ''Select a spline file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''meshchanfile''), ''style'', ''edit'', ''callback'', '''', ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_selectcoreg = [ 'tmpmodel = get( findobj(gcbf, ''userdata'', ''meshfile'')    , ''string''); tmpmodel = tmpmodel{get( findobj(gcbf, ''userdata'', ''meshfile'')        , ''value'')};' ...
                       'tmploc2  = get( findobj(gcbf, ''userdata'', ''meshchanfile''), ''string''); tmploc2  = tmploc2{ get( findobj(gcbf, ''userdata'', ''meshchanfile'')    , ''value'')};' ...
                       'tmploc1  = get( gcbo, ''userdata'');' ...
                       'tmptransf = get( findobj(gcbf, ''userdata'', ''coregtext''), ''string'');' ...
                       '[tmp tmptransf] = coregister(tmploc1{1}, tmploc2, ''mesh'', tmpmodel, ''helpmsg'', ''on'',' ...
                       '                       ''chaninfo1'', tmploc1{2}, ''transform'', str2num(tmptransf));' ...
                       'if ~isempty(tmptransf), set( findobj(gcbf, ''userdata'', ''coregtext''), ''string'', num2str(tmptransf)); end;' ...
                       'clear tmpmodel tmploc2 tmploc1 tmp tmptransf;' ];    
    cb_helpload = [ 'warndlg2(strvcat(''If you have already generated a spline file for this channel location'',' ...
                                   '''structure, you may enter it here. Click on the "Use existing spline file or'',' ...
                                   '''structure" to activate the edit box first.''), ''Load file for headplot()'');' ];
    cb_helpcoreg = [ 'warndlg2(strvcat(''Your channel locations must be co-registered with a 3-D head mesh to be plotted.'',' ...
                                   '''If you are using one of the template location files, the "Talairach transformation matrix"'',' ...
                                   '''field will be filled automatically (just enter an output file name and press "OK").'',' ...
                                   '''Otherwise press the "Manual coreg." button to perform co-registration.''), ''Load file for headplot()'');' ];
    cb_selectmesh = [ 'set(findobj(gcbf, ''userdata'', ''meshchanfile''), ''value'', get(gcbo, ''value''));' ...
                      'set(findobj(gcbf, ''userdata'', ''meshfile'')    , ''value'', get(gcbo, ''value''));' ...
                      'tmpdat = get(gcbf, ''userdata'');' ...
                      'set(findobj(gcbf, ''userdata'', ''coregtext''), ''string'', num2str(tmpdat{get(gcbo, ''value'')}));' ];
    defaultmat = { 'mheadnew.mat' 'colin27headmesh.mat' };
    defaultloc = { 'mheadnew.xyz' 'colin27headmesh.xyz' };
    defaulttransform = { transform [0 -15 -15 0.05 0 -1.57 100 88 110] };
    if iseeglabdeployed
        defaultmat = fullfile(eeglabexefolder, defaultmat);
        defaultloc = fullfile(eeglabexefolder, defaultloc);
    end;
    userdatatmp = { EEG.chanlocs EEG.chaninfo };
	txt = { { 'style' 'text'        'string' 'Co-register channel locations with head mesh and compute a mesh spline file (each scalp montage needs a headplot() spline file)' 'fontweight' 'bold' } ...
            { 'style' 'checkbox'    'string' 'Use the following spline file or structure' 'userdata' 'loadfile' 'tag' 'loadcb' 'callback' cb_load 'value' ~compute_file } ...
            { 'style' 'edit'        'string' fastif(typeplot, EEG.splinefile, EEG.icasplinefile)  'userdata' 'load' 'tag' 'load' 'enable' enableload } ...
            { 'style' 'pushbutton'  'string' 'Browse'        'callback' cb_browseload                               'tag' 'load' 'enable' enableload } ... 
            { 'style' 'pushbutton'  'string' 'Help'          'callback' cb_helpload } ...
            { 'style' 'checkbox'    'string' 'Or (re)compute a new spline file named:' 'tag' 'compcb' 'callback' cb_comp 'value' compute_file } ...
            { 'style' 'edit'        'string' [fullfile(pwd, EEG.filename(1:length(EEG.filename)-3)),'spl'] 'userdata' 'coregfile'  'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Browse'        'callback' cb_browsecomp                               'tag' 'comp' 'enable' enablecomp } ... 
            { 'style' 'pushbutton'  'string' 'Help'          'callback' cb_helpcoreg } ...
            { 'style' 'text'        'string' '            3-D head mesh file'                                       'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'popupmenu'   'string' defaultmat      'userdata' 'meshfile' 'callback' cb_selectmesh         'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Browse other'        'callback' cb_browsemesh                         'tag' 'comp' 'enable' enablecomp } ... 
            { } ... 
            { 'style' 'text'        'string' '            Mesh associated channel file'                             'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'popupmenu'   'string' defaultloc      'userdata' 'meshchanfile'  'callback' cb_selectmesh    'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Browse other'        'callback' cb_browsemeshchan                     'tag' 'comp' 'enable' enablecomp } ... 
            { } ... 
            { 'style' 'text'        'string' '            Talairach-model transformation matrix'                    'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'edit'        'string' num2str(transform) 'userdata' 'coregtext'                              'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Manual coreg.' 'callback' cb_selectcoreg 'userdata' userdatatmp       'tag' 'comp' 'enable' enablecomp } ... 
            { } ...
            { } ...
            { 'style' 'text'        'string' 'Plot interpolated activity onto 3-D head' 'fontweight' 'bold' } ...
            { 'style' 'text' 'string' txt } ...
	        { 'style' 'edit' 'string' fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))] ) } { } ...
	        { 'style' 'text' 'string' 'Plot title:' } ...
	        { 'style' 'edit' 'string' [ fastif( typeplot, 'ERP scalp maps of dataset:', 'Components of dataset: ') ...
                                        fastif(~isempty(EEG.setname), EEG.setname, '') ] } { }  ...
	        { 'style' 'text' 'string' 'Plot geometry (rows,columns): (Default [] = near square)' } ...
	        { 'style' 'edit' 'string' '' } { }  ...
	        { 'style' 'text' 'string' 'Other headplot options (See >> help headplot):' } ...
	        { 'style' 'edit' 'string' '' }  { } };
        
    % plot GUI and protect parameters
    % -------------------------------
    geom = { [1] [1.3 1.6 0.5 0.5 ] [1.3 1.6 0.5 0.5 ] [1.3 1.6 0.6 0.4 ] [1.3 1.6 0.6 0.4 ] [1.3 1.6 0.6 0.4 ] ...
             [1] [1] [1.5 1 0.5] [1.5 1  0.5] [1.5 1  0.5] [1.5 1 0.5] };
    optiongui = { 'uilist', txt, 'title', fastif( typeplot, 'ERP head plot(s) -- pop_headplot()', ...
                       'Component head plot(s) -- pop_headplot()'), 'geometry', geom 'userdata' defaulttransform };
	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', 'noclose', optiongui{:});
    if isempty(result), return; end;
    if ~isempty(get(0, 'currentfigure')) currentfig = gcf; else return; end;
    
    while test_wrong_parameters(currentfig)
    	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
        if isempty(result), return; end;
    end;
    close(currentfig);
    
    % decode setup parameters
    % -----------------------
    options = {};
    if result{1},               options = { options{:} 'load'    result{2} };
    else
        if ~isstr(result{5})    result{5} = defaultmat{result{5}}; end;
        if isempty(result{7})   setupopt = { result{4} 'meshfile' result{5} };  % no coreg
        else                    setupopt = { result{4} 'meshfile' result{5} 'transform' str2num(result{7}) };
                                fprintf('Transformation matrix: %s\n', result{7});
        end;
        options = { options{:} 'setup' setupopt };
        if ~strcmpi(result{5}, 'mheadnew.mat'), EEG.headplotmeshfile = result{5}; 
        else EEG.headplotmeshfile = ''; end;
    end;
    
    % decode other parameters
    % -----------------------
    arg2 = eval( [ '[' result{8} ']' ] );
	if length(arg2) > EEG.nbchan
		tmpbut = questdlg2(['This will draw ' int2str(length(arg2)) ' plots. Continue ?'], '', 'Cancel', 'Yes', 'Yes');
		if strcmp(tmpbut, 'Cancel'), return; end;
	end;
    if length(arg2) == 0, error('please choose a latency(s) to plot'); end
	topotitle  = result{9};
	rowcols    = eval( [ '[ ' result{10} ' ]' ] );
    tmpopts    = eval( [ '{ ' result{11} ' }' ] );
    if ~isempty(tmpopts)
        options    = { options{:} tmpopts{:} };
    end;
	if size(arg2(:),1) == 1, figure; end;
else % Pass along parameters and bypass GUI input
    options = varargin;
end;

% Check if pop_headplot input 'colorbar' was called, and don't send it to headplot
loc = strmatch('colorbar', options(1:2:end), 'exact');
loc = loc*2-1;
if ~isempty(loc)
    colorbar_switch = strcmp('on',options{ loc+1 });
    options(loc:loc+1) = [];
else
    colorbar_switch = 1;
end 

% read or generate file if necessary
% ----------------------------------
pop_options = options;
loc = strmatch('load', options(1:2:end)); loc = loc*2-1;
if ~isempty(loc)
    if typeplot
        EEG.splinefile    = options{ loc+1 };
    else            
        EEG.icasplinefile = options{ loc+1 };
    end;
    options(loc:loc+1) = [];
end;
loc = strmatch('setup', options(1:2:end)); loc = loc*2-1;
if ~isempty(loc)
    if typeplot
        headplot('setup', EEG.chanlocs, options{loc+1}{1}, 'chaninfo', EEG.chaninfo, options{ loc+1 }{2:end});
        EEG.splinefile    = options{loc+1}{1};
    else
        headplot('setup', EEG.chanlocs, options{loc+1}{1}, 'chaninfo', EEG.chaninfo, 'ica', 'on', options{ loc+1 }{2:end});
        EEG.icasplinefile = options{loc+1}{1};
    end;
    options(loc:loc+1) = [];
    compute_file = 1;
else
    compute_file = 0;
end;

% search for existing file if necessary
% -------------------------------------
if typeplot == 1 % ********** data plot
    fieldname    = 'splinefile';        
    if isempty(EEG.splinefile)            
        if length(EEG.icachansind) == EEG.nbchan & ~isempty(EEG.icasplinefile)
            EEG.splinefile = EEG.icasplinefile;
        end;
    end;
else % ************* Component plot       
    fieldname    = 'icasplinefile';
    if isempty(EEG.icasplinefile)
        if length(EEG.icachansind) == EEG.nbchan & ~isempty(EEG.splinefile)
            EEG.icasplinefile = EEG.splinefile;
        end;
    end;
end;

% headplot mesh file
% ------------------
if isfield(EEG, 'headplotmeshfile')
    if ~isempty(EEG.headplotmeshfile)
        options = { options{:} 'meshfile' EEG.headplotmeshfile };
    end;
end;

% check parameters
% ----------------
if ~exist('topotitle')  
    topotitle = '';
end;    
if typeplot
    if isempty(EEG.splinefile)
        error('Pop_headplot: cannot find spline file, aborting...');
    end;
else
    if isempty(EEG.icasplinefile)
        error('Pop_headplot: cannot find spline file, aborting...');
    end;
end;    
SIZEBOX = 150;
nbgraph = size(arg2(:),1);
if ~exist('rowcols') | isempty(rowcols) | rowcols == 0
    rowcols(2) = ceil(sqrt(nbgraph));
    rowcols(1) = ceil(nbgraph/rowcols(2));
end;    

fprintf('Plotting...\n');

% determine the scale for plot of different times (same scales)
% -------------------------------------------------------------
if typeplot
	SIGTMP = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
	pos = round( (arg2/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
	indexnan = find(isnan(pos));
	nanpos = find(isnan(pos));
	pos(nanpos) = 1;
	SIGTMPAVG = mean(SIGTMP(:,pos,:),3);
	SIGTMPAVG(:, nanpos) = NaN;
	maxlim = max(SIGTMPAVG(:));
	minlim = min(SIGTMPAVG(:));
	maplimits = max(maxlim, -minlim);
    maplimits = maplimits*1.1;
    maplimits = [ -maplimits maplimits ];
else
    maplimits = [-1 1];
end;
	
% plot the graphs
% ---------------
counter = 1;
disp('IMPORTANT NOTICE: electrodes are projected to the head surface so their location');
disp('                  might slightly differ from the one they had during coregistration ');
for index = 1:size(arg2(:),1)
	if nbgraph > 1
        if mod(index, rowcols(1)*rowcols(2)) == 1
            if index> 1, a = textsc(0.5, 0.05, topotitle); set(a, 'fontweight', 'bold'); end;
        	figure;
        	pos = get(gcf,'Position');
			posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
			posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
			set(gcf,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
			try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
        end;    
		subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1);
	end;

	if ~isnan(arg2(index))
		if typeplot
            headplot( SIGTMPAVG(:,index), EEG.splinefile, 'maplimits', maplimits, options{:});
			if nbgraph == 1, title( topotitle );
			else title([int2str(arg2(index)) ' ms']);
			end;
		else
            if arg2(index) < 0
                headplot( -EEG.icawinv(:, -arg2(index)), EEG.icasplinefile, options{:});
            else	
                headplot( EEG.icawinv(:, arg2(index)), EEG.icasplinefile, options{:});
            end;    			
			if nbgraph == 1, title( topotitle );
			else title(['' int2str(arg2(index))]);
			end;
		end;
		drawnow;
		axis equal; 
		rotate3d off;
    else
        axis off
    end
end

% Draw colorbar
if colorbar_switch
    if nbgraph == 1
        ColorbarHandle = cbar(0,0,[maplimits(1) maplimits(2)]); 
        pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
        set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
    else
        cbar('vert',0,[maplimits(1) maplimits(2)]);
    end
    if ~typeplot    % Draw '+' and '-' instead of numbers for colorbar tick labels
        tmp = get(gca, 'ytick');
        set(gca, 'ytickmode', 'manual', 'yticklabelmode', 'manual', 'ytick', [tmp(1) tmp(end)], 'yticklabel', { '-' '+' });
    end
end
        
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;


if nbgraph> 1, 
    a = textsc(0.5, 0.05, topotitle); 
    set(a, 'fontweight', 'bold');
    axcopy(gcf, [ 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position'''');' ...
                  'set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); rotate3d(gcf); clear postmp;' ]);
end;

% generate output command
% -----------------------
com = sprintf('pop_headplot(%s, %d, %s, ''%s'', [%s], %s);', inputname(1), typeplot, vararg2str(arg2), ...
              topotitle, int2str(rowcols), vararg2str( pop_options ) );
if compute_file, com = [ 'EEG = ' com ]; end;
if nbgraph== 1,  com = [ 'figure; ' com ]; rotate3d(gcf); end;

return;

% test for wrong parameters
% -------------------------
function bool = test_wrong_parameters(hdl)

    bool = 0;
    loadfile = get( findobj( hdl, 'userdata', 'loadfile')     , 'value' );
    
    textlines = '';
    if ~loadfile
        coreg1   = get( findobj( hdl, 'userdata', 'coregtext')    , 'string' );
        coreg3   = get( findobj( hdl, 'userdata', 'coregfile')    , 'string' );
        if isempty(coreg1)
            textlines = strvcat('You must co-register your channel locations with the head model.', ...
                                'This is an easy process: Press the "Manual coreg." button in the', ...
                                'right center of the pop_headplot() window and follow instructions.',...
                                'To bypass co-registration (not recommended), enter', ...
                                '"0 0 0 0 0 0 1 1 1" as the "Tailairach transformation matrix.');
            bool = 1;
        end;
        if isempty(coreg3)
            textlines = strvcat(textlines, ' ', 'You need to enter an output file name.');
            bool = 1;
        end;
        
        if bool
            warndlg2( textlines, 'Error');
        end;
    end;
