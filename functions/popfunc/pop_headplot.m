% pop_headplot() - plot one or more spherically-splined EEG field maps 
%                  using a semi-realistic 3-D head model. Requires a
%                  spline file, which is created first if not found.
%                  This may take some time, but does not need to be 
%                  done again for this channel locations montage. A wait 
%                  bar will pop up to indicate how much time remains.
% Usage:
%   >> EEGOUT = pop_headplot( EEG, typeplot, ...
%                    latencies/components, title, rowscols, 'key', 'val' ...);
% Inputs:
%   EEG        - EEG dataset structure
%   typeplot   - 1=channel, 0=component {Default: 1}
%   latencies/components  - If channels, array of epoch mean latencies (in ms),
%                Else, for components, array of component indices to plot.
%   title      - Plot title
%   rowscols   - Vector of the form [m,n] where m is total vertical tiles and n 
%                horizontal tiles per page. If the number of maps exceeds m*n,
%                multiple figures will be produced {def|0 -> 1 near-square page}
%
% Optional inputs:
%   'setup'    - ['setupfile.spl'] Make the headplot spline file
%   'load'     - ['setupfile.spl'] Load the headplot spline file
%   others...  - Other headplot options. See >> help headplot
%
% Output:
%   EEGOUT - EEG dataset, possibly with a new or modified splinefile. 
%
% Note:
%   A new figure is created only when the pop_up window is called or when
%   several channels/components are plotted. Therefore you may call this 
%   command to draw single 3-D topographic maps in an existing figure.     
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 March 2002
%
% See also: headplot(), eegplot()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.36  2006/01/24 19:30:50  arno
% name of dataset
%
% Revision 1.35  2006/01/21 00:44:04  arno
% gui error
%
% Revision 1.34  2006/01/21 00:34:14  arno
% fix coregistration etc...
%
% Revision 1.33  2006/01/20 23:47:01  arno
% fix coregistration
%
% Revision 1.32  2006/01/20 22:11:23  arno
% title, coregistration etc...
%
% Revision 1.31  2006/01/10 00:41:29  arno
% trying to allow coregistration - not finished
%
% Revision 1.30  2005/12/01 20:06:16  arno
% search for existing file
%
% Revision 1.29  2005/11/30 23:38:57  arno
% typo
%
% Revision 1.28  2005/11/30 23:29:43  arno
% taking into account icachansind and channel location file orientation; reprogrammed options cell array
%
% Revision 1.27  2004/08/31 13:50:20  scott
% edited printed messages -sm
%
% Revision 1.26  2004/07/26 19:01:57  arno
% debug detection of old channel location file
%
% Revision 1.25  2004/07/01 22:44:45  arno
% detect if using old spline file
%
% Revision 1.24  2003/12/17 01:07:11  arno
% debug last
%
% Revision 1.23  2003/12/17 01:06:20  arno
% processing empty coordinates
%
% Revision 1.22  2003/05/12 16:11:18  arno
% debuging output command
%
% Revision 1.21  2003/05/12 15:59:46  arno
% debug last
%
% Revision 1.20  2003/05/12 15:54:20  arno
% debuging output command if creating spline file
%
% Revision 1.19  2003/05/10 02:36:02  arno
% compress output command
%
% Revision 1.18  2003/05/10 02:16:41  arno
% debug autoscale
%
% Revision 1.17  2003/03/12 06:34:46  scott
% header edits
%
% Revision 1.16  2003/03/12 03:21:38  arno
% help button
%
% Revision 1.15  2002/10/15 17:03:37  arno
% drawnow
%
% Revision 1.14  2002/08/27 00:38:22  arno
% more optimal auto location
%
% Revision 1.13  2002/08/20 00:05:59  arno
% adding test for plotting a large number of components
%
% Revision 1.12  2002/08/19 23:59:09  arno
% reducing message
%
% Revision 1.11  2002/08/12 16:32:22  arno
% inputdlg2
%
% Revision 1.10  2002/08/12 02:46:45  arno
% inputdlg2
%
% Revision 1.9  2002/08/12 02:44:37  arno
% inputdlg2
%
% Revision 1.8  2002/08/12 01:37:37  arno
% color
%
% Revision 1.7  2002/08/11 22:15:06  arno
% color
%
% Revision 1.6  2002/07/25 18:41:02  arno
% same
%
% Revision 1.5  2002/07/25 18:40:06  arno
% debugging
%
% Revision 1.4  2002/07/25 18:29:07  arno
% change 3d rotate options
%
% Revision 1.3  2002/04/18 15:49:09  scott
% editted msgs -sm
%
% Revision 1.2  2002/04/07 20:47:23  scott
% worked on no-spline-file msg -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, com] = pop_headplot( EEG, typeplot, arg2, topotitle, rowcols, varargin);

com = '';
if nargin < 1
   help pop_headplot;
   return;
end;   

if isempty(EEG.chanlocs)
    error('Pop_headplot: this dataset does not contain channel locations. Use menu item: Edit > Dataset info');
end;

if nargin < 3
    % remove old spline file
    % ----------------------
    if isfield(EEG, 'splinefile') 
        if ~isempty(EEG.splinefile)
            splfile = dir(EEG.splinefile);
            byteperelec = splfile.bytes/EEG.nbchan;
            if byteperelec/EEG.nbchan < 625, % old head plot file
                EEG.splinefile = [];
                disp('Warning: Wrong montage or old-versin spline file version detected and removed; new spline file required');
            end;
        end;
    end;
    
    % show the file be recomputed
    % ---------------------------
    compute_file = 0;
    if typeplot == 1 % ********** data plot
        fieldname    = 'splinefile';        
        if isempty(EEG.splinefile)            
            if length(EEG.icachansind) == EEG.nbchan & ~isempty(EEG.icasplinefile)
                EEG.splinefile = EEG.icasplinefile;
            else
                compute_file = 1;
            end;
        end;
    else % ************* Component plot       
        fieldname    = 'icasplinefile';
        if isempty(EEG.icasplinefile)
            if length(EEG.icachansind) == EEG.nbchan & ~isempty(EEG.splinefile)
                EEG.icasplinefile = EEG.splinefile;
            else
                compute_file = 1;
            end;
        end;
    end;
            
    if compute_file
        
		 warndlg2( strvcat('headplot() must generate a spline file the first', ...
		                 'time it is called and after changes in the channel location file.', ...
		                 'You must also coregister your channel location structure with the', ...
                         'head template.'), 'Headplot() warning');
    else
    	pop_options      = {};
    end;
    
 	% graphic interface
	% -----------------
    autocoregbem = [ -0.732155 7.58141 11.8939 -0.0249659 0.0148571 0.0227427 0.932423 0.918943 0.793166 ];
    autocoregsph = [ -0.31937  -5.96928 13.1812 0.0509311 0.0172127 -1.55007  1.08221  1.00037  0.923518 ];
    transform = [];
    if isfield(EEG.chaninfo, 'filename')
        if ~isempty(findstr(lower(EEG.chaninfo.filename), 'standard-10-5-cap385')), transform = autocoregsph; end;
        if ~isempty(findstr(lower(EEG.chaninfo.filename), 'standard_1005')),        transform = autocoregbem; end;
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
    cb_browseload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''load''), ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_browsecomp = [ '[filename, filepath] = uiputfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''coregfile''), ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_browsemesh = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''meshfile''), ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_browsemeshchan = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj( gcbf, ''userdata'', ''meshchanfile''), ''string'', fullfile(filepath,filename));' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    cb_selectcoreg = [ 'tmpmodel = get( findobj(gcbf, ''userdata'', ''meshfile'')    , ''string'');' ...
                       'tmploc2  = get( findobj(gcbf, ''userdata'', ''meshchanfile''), ''string'');' ...
                       'tmploc1  = get( gcbo, ''userdata'');' ...
                       'tmptransf = get( findobj(gcbf, ''userdata'', ''coregtext''), ''string'');' ...
                       '[tmp tmptransf] = coregister(tmploc1, tmploc2, ''mesh'', tmpmodel, ''helpmsg'', ''on'',' ...
                       '                       ''transform'', str2num(tmptransf));' ...
                       'if ~isempty(tmptransf), set( findobj(gcbf, ''userdata'', ''coregtext''), ''string'', num2str(tmptransf)); end;' ...
                       'clear tmpmodel tmploc2 tmploc1 tmp tmptransf;' ];    
    cb_helpload = [ 'warndlg2(strvcat(''If you have already generated a spline file for this channel location'',' ...
                                   '''structure, you may enter it here. Click on the "Use existing spline file or'',' ...
                                   '''structure" to activate the edit box first.''), ''Load file for headplot()'');' ];
    cb_helpcoreg = [ 'warndlg2(strvcat(''Your channel locations must be co-registered with a 3D head to be plotted.'',' ...
                                   '''If you are using one of the template location file, the "Talairach transformation matrix"'',' ...
                                   '''field will be filled automatically (just enter an output file name and press OK).'',' ...
                                   '''Otherwise press the "Manual coreg." button to perform the coregistration.''), ''Load file for headplot()'');' ];
    
	txt = { { 'style' 'text'        'string' 'Coregister channel locations with head and compute spline file (done only once)' 'fontweight' 'bold' } ...
            { 'style' 'checkbox'    'string' 'Use the following spline file or structure' 'userdata' 'loadfile' 'tag' 'loadcb' 'callback' cb_load 'value' ~compute_file } ...
            { 'style' 'edit'        'string' fastif(typeplot, EEG.splinefile, EEG.icasplinefile)  'userdata' 'load' 'tag' 'load' 'enable' enableload } ...
            { 'style' 'pushbutton'  'string' 'Browse'        'callback' cb_browseload                               'tag' 'load' 'enable' enableload } ... 
            { 'style' 'pushbutton'  'string' 'Help'          'callback' cb_helpload } ...
            { 'style' 'checkbox'    'string' 'Or (re)compute a new spline file named:' 'tag' 'compcb' 'callback' cb_comp 'value' compute_file } ...
            { 'style' 'edit'        'string' ''              'userdata' 'coregfile'                                 'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Browse'        'callback' cb_browsecomp                               'tag' 'comp' 'enable' enablecomp } ... 
            { 'style' 'pushbutton'  'string' 'Help'          'callback' cb_helpcoreg } ...
            { 'style' 'text'        'string' '            3-D head mesh file'                                       'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'edit'        'string' 'mheadnew.mat'  'userdata' 'meshfile'                                  'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Browse'        'callback' cb_browsemesh                               'tag' 'comp' 'enable' enablecomp } ... 
            { } ... 
            { 'style' 'text'        'string' '            Mesh associated channel file'                           'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'edit'        'string' 'mheadnew.xyz'  'userdata' 'meshchanfile'                              'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Browse'        'callback' cb_browsemeshchan                           'tag' 'comp' 'enable' enablecomp } ... 
            { } ... 
            { 'style' 'text'        'string' '            Talairach transformation matrix'                          'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'edit'        'string' num2str(transform) 'userdata' 'coregtext'                              'tag' 'comp' 'enable' enablecomp } ...
            { 'style' 'pushbutton'  'string' 'Manual coreg.' 'callback' cb_selectcoreg 'userdata' EEG.chanlocs      'tag' 'comp' 'enable' enablecomp } ... 
            { } ...
            { } ...
            { 'style' 'text'        'string' 'Plot interpolated activity onto 3-D head' 'fontweight' 'bold' } ...
            { 'style' 'text' 'string' txt } ...
	        { 'style' 'edit' 'string' fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))] ) } ...
	        { 'style' 'text' 'string' 'Plot title:' } ...
	        { 'style' 'edit' 'string' [ fastif( typeplot, 'ERP scalp maps of dataset:', 'Components of dataset: ') ...
                                        fastif(~isempty(EEG.setname), EEG.setname, '') ] } ...
	        { 'style' 'text' 'string' 'Plot geometry (rows,columns): (Default [] = near square)' } ...
	        { 'style' 'edit' 'string' '' } ...
	        { 'style' 'text' 'string' '-> headplot() options  (See >> help headplot):' } ...
	        { 'style' 'edit' 'string' '' } };
        
    % plot GUI and protect parameters
    % -------------------------------
    geom = { [1] [1.3 1 0.5 0.5 ] [1.3 1 0.5 0.5 ] [1.3 1 0.5 0.5 ] [1.3 1 0.5 0.5 ] [1.3 1 0.5 0.5 ] [1] [1] [2 1] [2 1] [2 1] [2 1] };
    optiongui = { 'uilist', txt, 'title', fastif( typeplot, 'ERP head plot(s) -- pop_headplot()', ...
                       'Component head plot(s) -- pop_headplot()'), 'geometry', geom };
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
        if result{7}            options = { options{:} 'setup' { result{4} 'meshfile' result{5} } }; % no coreg
        else                    options = { options{:} 'setup' { result{4} 'meshfile' result{5} 'transform' result{6} } };
        end;
    end;
    
    % decode other parameters
    % -----------------------
    arg2 = eval( [ '[' result{8} ']' ] );
	if length(arg2) > EEG.nbchan
		tmpbut = questdlg2(['This will draw ' int2str(length(arg2)) ' plots. Continue ?'], '', 'Cancel', 'Yes', 'Yes');
		if strcmp(tmpbut, 'Cancel'), return; end;
	end;
	topotitle  = result{9};
	rowcols    = eval( [ '[ ' result{10} ' ]' ] );
    tmpopts    = eval( [ '{ ' result{11} ' }' ] );
    if ~isempty(tmpopts)
        options    = { options{:} tmpopts{:} };
    end;
	if size(arg2(:),1) == 1, figure; end;
else
    options = varargin;
end;

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
		if index == size(arg2(:),1)
	        pos = get(gca,'position');
	        q = [pos(1) pos(2) 0 0];
	        s = [pos(3) pos(4) pos(3) pos(4)];
	        col = colormap;
            if nbgraph > 1
                ax = subplot('position', [1.1 0 .05 1].*s+q);
	        else 
                ax = subplot('position', [1 0 .05 1].*s+q);          
            end;
            col = col(1:end-3,:);
            imagesc([], [maplimits(1) 0 maplimits(2)], reshape(col,size(col,1),1,3));
            set(gca, 'xtick', [], 'yaxislocation', 'right', 'ydir', 'normal');
	    end;	   
    else
        axis off
    end;
end;
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
            textlines = strvcat('You must coregister your channel locations with the head model.', ...
                                'It is an easy process (just press the "Manual coreg. button and', ...
                                'follow instructions). To bypass coregistration (not recommended), enter', ...
                                '"0 0 0 0 0 0 1 1 1" in the coregistration edit box.');
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
