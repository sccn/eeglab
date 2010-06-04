% pop_topoplot() - Plot scalp map(s) in a figure window. If number of input
%                  arguments is less than 3, pop up an interactive query window.
%                  Makes (possibly repeated) calls to topoplot().
% Usage:
%   >> pop_topoplot( EEG); % pops up a parameter query window
%   >> pop_topoplot( EEG, typeplot, items, title, plotdip, options...); % no pop-up
%
% Inputs:
%   EEG      - Input EEG dataset structure (see >> help eeglab)
%   typeplot - 1-> Plot channel ERP maps; 0-> Plot component maps {default:1}.
%
% Commandline inputs also set in pop-up window:
%   items    - [array] If typeplot==1 (ERP maps), within-epoch latencies 
%              (ms) at which to plot the maps. If typeplot==0 (component
%              maps), component indices to plot. In this case,
%              negative map indices -> invert map polarity, or 
%              NaN -> leave a blank subplot. (Ex: [1 -3 NaN 4])
%   title    - Plot title.
%   rowscols - Vector of the form [m,n] giving [rows, cols] per page.
%              If the number of maps exceeds m*n, multiple figures 
%              are produced {default|0 -> one near-square page}.
%   plotdip  - [0|1] plot associated dipole(s) for scalp map if present
%              in dataset.
%
% Optional Key-Value Pair Inputs
%   'colorbar' - ['on' or 'off'] Switch to turn colorbar on or off. {Default: 'on'}
%   options  - optional topoplot() arguments. Separate using commas. 
%              Example 'style', 'straight'. See >> help topoplot
%              for further details. {default: none}
%
% Note:
%   A new figure is created automatically only when the pop_up window is 
%   called or when more than one page of maps are plotted. Thus, this 
%   command may be used to draw topographic maps in a figure sub-axis.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: topoplot(), eeglab()

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

% 01-25-02 reformated help & license -ad 
% 02-15-02 text interface editing -sm & ad 
% 02-16-02 added axcopy -ad & sm
% 03-18-02 added title -ad & sm

function com = pop_topoplot( EEG, typeplot, arg2, topotitle, rowcols, varargin);

com = '';
if nargin < 1
   help pop_topoplot;
   return;
end;   
if nargin < 2   
   typeplot = 1;
end;
if typeplot == 0 & isempty(EEG.icasphere)
   disp('Error: no ICA data for this set, first run ICA'); return;
end;   
if isempty(EEG.chanlocs)
   disp('Error: cannot plot topography without channel location file'); return;
end;   

if nargin < 3
	% which set to save
	% -----------------
	if typeplot
		txtwhat2plot1 = 'Plotting ERP scalp maps at these latencies';
        txtwhat2plot2 = sprintf('(range: %d to %d ms, NaN -> empty):', ...
                                round(EEG.xmin*1000), round(EEG.xmax*1000));
        editwhat2plot = [''];
	else
		txtwhat2plot1 = 'Component numbers';
		txtwhat2plot2 = '(negate index to invert component polarity; NaN -> empty subplot; Ex: -1 NaN 3)';
        editwhat2plot = ['1:' int2str(size(EEG.icaweights,1))];
 	end;	
        if EEG.nbchan > 64, 
            elecdef = ['''electrodes'', ''off''']; 
        else, 
            elecdef = ['''electrodes'', ''on''']; 
        end;
    uilist = { { 'style'   'text'     'string'    txtwhat2plot1 } ...
               { 'style'   'edit'     'string'    editwhat2plot } ...
               { 'style'   'text'     'string'    txtwhat2plot2 } ...
               { } ...
               { 'style'   'text'     'string'    'Plot title' } ...
               { 'style'   'edit'     'string'    fastif(~isempty(EEG.setname), [EEG.setname], '') } ...
               { 'style'   'text'     'string'    'Plot geometry (rows,col.); [] -> near square' } ...
               { 'style'   'edit'     'string'    '[]' } ...
               { 'style'   'text'     'string'    'Plot associated dipole(s) (if present)' } ...
               { 'style'   'checkbox' 'string'    '' } { } ...
               { } ...
               { 'style'   'text'     'string'    [ '-> Additional topoplot()' fastif(typeplot,'',' (and dipole)') ...
                                                  ' options (see Help)' ] } ...
               { 'style'   'edit'     'string'    elecdef } };
    uigeom = { [1.5 1] [1] [1] [1.5 1] [1.5 1] [1.55 0.2 0.8] [1] [1] [1] };
    if typeplot
        uilist(9:11) = [];
        uigeom(6) = [];
    end;
    guititle = fastif( typeplot, 'Plot ERP scalp maps in 2-D -- pop_topoplot()', ...
                       'Plot component scalp maps in 2-D -- pop_topoplot()');
    
    result = inputgui( uigeom, uilist, 'pophelp(''pop_topoplot'')', guititle, [], 'normal');
	if length(result) == 0 return; end;
	
    % reading first param
    % -------------------
    arg2   	     = eval( [ '[' result{1} ']' ] );
	if length(arg2) > EEG.nbchan
		tmpbut = questdlg2(...
                  ['This involves drawing ' int2str(length(arg2)) ' plots. Continue ?'], ...
                         '', 'Cancel', 'Yes', 'Yes');
		if strcmp(tmpbut, 'Cancel'), return; end;
	end;
    if isempty(arg2), error('Nothing to plot; enter parameter in first edit box'); end;
    
    % reading other params
    % --------------------
	topotitle   = result{2};
	rowcols     = eval( [ '[' result{3} ']' ] );
    if typeplot
        plotdip = 0;
        try, options      = eval( [ '{ ' result{4} ' }' ]);
        catch, error('Invalid scalp map options'); end;
    else
        plotdip     = result{4};
        try, options      = eval( [ '{ ' result{5} ' }' ]);
        catch, error('Invalid scalp map options'); end;
    end;        
    if length(arg2) == 1, 
      figure('paperpositionmode', 'auto'); curfig=gcf; 
      try, icadefs; 
         set(curfig, 'color', BACKCOLOR); 
      catch, end; 
    end;
else
    if ~isempty(varargin) & isnumeric(varargin{1})
        plotdip = varargin{1};
        varargin = varargin(2:end);
    else
        plotdip = 0;
    end;
    options = varargin;
end;

% additional options
% ------------------
outoptions = { options{:} }; % for command
options    = { options{:} 'masksurf' 'on' };

% find maplimits
% --------------
maplimits = [];
for i=1:2:length(options)
    if isstr(options{i})
        if strcmpi(options{i}, 'maplimits')
            maplimits = options{i+1};
            options(i:i+1) = [];
            break;
        end;
    end;
end;

nbgraph = size(arg2(:),1);
if ~exist('topotitle')
    topotitle = '';
end;    
if ~exist('rowcols') | isempty(rowcols) | rowcols == 0
    rowcols(2) = ceil(sqrt(nbgraph));
    rowcols(1) = ceil(nbgraph/rowcols(2));
end;    

SIZEBOX = 150;

fprintf('Plotting...\n');
if isempty(EEG.chanlocs)
	fprintf('Error: set has no channel location file\n');
	return;
end;

% Check if pop_topoplot input 'colorbar' was called, and don't send it to topoplot
loc = strmatch('colorbar', options(1:2:end), 'exact');
loc = loc*2-1;
if ~isempty(loc)
    colorbar_switch = strcmp('on',options{ loc+1 });
    options(loc:loc+1) = [];
else
    colorbar_switch = 1;
end 

% determine the scale for plot of different times (same scales)
% -------------------------------------------------------------
if typeplot
    SIGTMP = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
    pos = round( (arg2/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
    nanpos = find(isnan(pos));
    pos(nanpos) = 1;
    SIGTMPAVG = mean(SIGTMP(:,pos,:),3);
    SIGTMPAVG(:, nanpos) = NaN;
    if isempty(maplimits)
        maxlim = max(SIGTMPAVG(:));
        minlim = min(SIGTMPAVG(:));
        maplimits = [ -max(maxlim, -minlim) max(maxlim, -minlim)];
    end;
else
    if isempty(maplimits)
        maplimits = 'absmax';
    end;
end;

if plotdip
    if strcmpi(EEG.dipfit.coordformat, 'CTF')
        disp('Cannot plot dipole on scalp map for CTF MEG data');
    end;
end;

% plot the graphs
% ---------------
counter = 1;
countobj = 1;
allobj = zeros(1,1000);
curfig = get(0, 'currentfigure');
if isfield(EEG, 'chaninfo'), options = { options{:} 'chaninfo' EEG.chaninfo }; end

for index = 1:size(arg2(:),1)
	if nbgraph > 1
        if mod(index, rowcols(1)*rowcols(2)) == 1
            if index> 1, figure(curfig); a = textsc(0.5, 0.05, topotitle); set(a, 'fontweight', 'bold'); end;
        	curfig = figure('paperpositionmode', 'auto');
			pos = get(curfig,'Position');
			posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
			posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
			set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
			try, icadefs; set(curfig, 'color', BACKCOLOR); catch, end;
        end;    
		curax = subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1);
        set(curax, 'visible', 'off')
    end;

	% add dipole location if present
    % ------------------------------
    dipoleplotted = 0;
    if plotdip && typeplot == 0
        if isfield(EEG, 'dipfit') & isfield(EEG.dipfit, 'model')
            if length(EEG.dipfit.model) >= index & ~strcmpi(EEG.dipfit.coordformat, 'CTF')
                %curpos = EEG.dipfit.model(arg2(index)).posxyz/EEG.dipfit.vol.r(end);
                curpos = EEG.dipfit.model(arg2(index)).posxyz;
                curmom = EEG.dipfit.model(arg2(index)).momxyz;
                try,
                    select = EEG.dipfit.model(arg2(index)).select;
                catch select = 0;
                end;
                if ~isempty(curpos)
                    if strcmpi(EEG.dipfit.coordformat, 'MNI') % from MNI to sperical coordinates
                        transform = pinv( sph2spm );
                        tmpres = transform * [ curpos(1,:) 1 ]'; curpos(1,:) = tmpres(1:3);
                        tmpres = transform * [ curmom(1,:) 1 ]'; curmom(1,:) = tmpres(1:3);
                        try, tmpres = transform * [ curpos(2,:) 1 ]'; curpos(2,:) = tmpres(1:3); catch, end;
                        try, tmpres = transform * [ curmom(2,:) 1 ]'; curmom(2,:) = tmpres(1:3); catch, end;
                    end;
                    curpos = curpos / 85;
                    if size(curpos,1) > 1 && length(select) == 2
                        dipole_index = find(strcmpi('dipole',options),1);
                        if  ~isempty(dipole_index) % if 'dipoles' is already defined in options{:}
                            options{dipole_index+1} = [ curpos(:,1:2) curmom(:,1:3) ];
                        else
                            options = { options{:} 'dipole' [ curpos(:,1:2) curmom(:,1:3) ] };
                        end
                        dipoleplotted = 1;
                    else
                        if any(curpos(1,:) ~= 0)
                            dipole_index = find(strcmpi('dipole',options),1);
                            if  ~isempty(dipole_index) % if 'dipoles' is already defined in options{:}
                                options{dipole_index+1} = [ curpos(1,1:2) curmom(1,1:3) ];
                            else
                                options = { options{:} 'dipole' [ curpos(1,1:2) curmom(1,1:3) ] };
                            end
                            dipoleplotted = 1;
                        end
                    end
                end
                if nbgraph ~= 1
                    dipscale_index = find(strcmpi('dipscale',options),1);
                    if ~isempty(dipscale_index) % if 'dipscale' is already defined in options{:}
                        options{dipscale_index+1} = 0.6;
                    else
                        options = {  options{:} 'dipscale' 0.6 };
                    end
                end
                %options = { options{:} 'dipsphere' max(EEG.dipfit.vol.r) };
            end
        end
    end
    
	% plot scalp map
    % --------------
    if index == 1
        addopt = { 'verbose', 'on' };
    else 
        addopt = { 'verbose', 'off' };
    end;
    %fprintf('Printing to figure %d.\n',curfig);
    options = {  'maplimits' maplimits options{:} addopt{:} };
    if ~isnan(arg2(index))
		if typeplot
            if nbgraph > 1, axes(curax); end;
            tmpobj = topoplot( SIGTMPAVG(:,index), EEG.chanlocs, options{:});
			if nbgraph == 1, 
                 figure(curfig); if nbgraph > 1, axes(curax); end;
                 title( [ 'Latency ' int2str(arg2(index)) ' ms from ' topotitle]);
			else 
                 figure(curfig); if nbgraph > 1, axes(curax); end; 
                 title([int2str(arg2(index)) ' ms'] );
			end;
		else
            if arg2(index) < 0
                 figure(curfig);  if nbgraph > 1, axes(curax); end;
                 tmpobj = topoplot( -EEG.icawinv(:, -arg2(index)), EEG.chanlocs, options{:} );
            else
                 figure(curfig);  if nbgraph > 1, axes(curax); end;
                 tmpobj = topoplot( EEG.icawinv(:, arg2(index)), EEG.chanlocs, options{:} );
            end;    			
			if nbgraph == 1, texttitle = [ 'IC ' int2str(arg2(index)) ' from ' topotitle];
			else             texttitle = ['' int2str(arg2(index))];
			end;
            if dipoleplotted, texttitle = [ texttitle ' (' num2str(EEG.dipfit.model(arg2(index)).rv*100,2) '%)']; end;
            figure(curfig);  if nbgraph > 1, axes(curax); end; title(texttitle);
		end;
        allobj(countobj:countobj+length(tmpobj)-1) = tmpobj;
        countobj = countobj+length(tmpobj);
		drawnow;
		axis square;
    else
    axis off
    end;
end

% Draw colorbar
if colorbar_switch
    if nbgraph == 1
        if ~isstr(maplimits)
            ColorbarHandle = cbar(0,0,[maplimits(1) maplimits(2)]);
        else
            ColorbarHandle = cbar(0,0,get(gca, 'clim'));
        end;
        pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
        set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
    elseif ~isstr(maplimits)
         cbar('vert',0,[maplimits(1) maplimits(2)]);
    else cbar('vert',0,get(gca, 'clim'));
    end
    if ~typeplot    % Draw '+' and '-' instead of numbers for colorbar tick labels
        tmp = get(gca, 'ytick');
        set(gca, 'ytickmode', 'manual', 'yticklabelmode', 'manual', 'ytick', [tmp(1) tmp(end)], 'yticklabel', { '-' '+' });
    end
end

if nbgraph> 1, 
   figure(curfig); a = textsc(0.5, 0.05, topotitle); 
   set(a, 'fontweight', 'bold'); 
end;
if nbgraph== 1, 
   com = 'figure;'; 
end;
set(allobj(1:countobj-1), 'visible', 'on');

figure(curfig);
axcopy(curfig, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); clear postmp;');

com = [com sprintf('pop_topoplot(%s,%d, %s);', ...
                   inputname(1), typeplot, vararg2str({arg2 topotitle rowcols plotdip outoptions{:} }))];
return;

		
