% pop_topoplot() - plot head of subjects with a pop-up window if only
%                  two arguments.
%
% Usage:
%   >> pop_topoplot( EEG, typeplot, latencies, title, options...);
%
% Inputs:
%   EEG        - dataset structure
%   typeplot   - 1=channel, 0=component (default:1)
%   latencies/components  - for EEG: array of latencies (in millisecond)
%                at which the head should be plotted.  
%                For components: array of index of components to plot. If
%                negative indices are entered, the function will plot the
%                components corresponding to the absolute value of these
%                indices and inverse the polarity of the plot. To leave blank 
%                subplot, enter nan in this vector. For components [1 nan -3 4]
%                plots component 1, blank, 3 in reverse polarity and 4.
%   title      - plot title.
%   rowscols   - Vector of the form [m,n] where m is total vertical tiles and n 
%                is horizontal tiles per page. If the number of maps exceeds m*n,
%                multiple figures will be produced {def|0 -> one near-square page}.
%   options    - TOPOPLOT options. Default is none. Separate the options
%                using comma. Example 'style', 'straight'. See TOPOPLOT  
%                help for further details. Default is no options. 
%
% Note:
%   A new figure is created only when the pop_up window is called or when
%   several channel/components are plotted, so you may call this command 
%   to draw topographic maps in a tiled window.     
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: topoplot(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.5  2002/08/12 16:31:20  arno
% inputdlg2
%
% Revision 1.4  2002/08/12 01:46:44  arno
% color
%
% Revision 1.3  2002/08/11 22:22:19  arno
% color
%
% Revision 1.2  2002/04/18 18:26:14  arno
% typo can not
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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
		%txt = sprintf('ERP scalp map at these latencies (from %d to %d ms):\n(NaN -> empty subplot)(Ex: -100 NaN 100)', round(EEG.xmin*1000), round(EEG.xmax*1000));
		txt = sprintf('ERP scalp map at these latencies (from %d to %d ms):', round(EEG.xmin*1000), round(EEG.xmax*1000));
	else
		%txt = ['Component numbers (negate index to invert component polarity):' 10 '(NaN -> empty subplot)(Ex: -1 NaN 3)'];
		txt = sprintf('ERP scalp map at these latencies (from %d to %d ms):', round(EEG.xmin*1000), round(EEG.xmax*1000));
	end;	
	txt = { txt ...
	        'Plot title:' ...
	        ['Plot geometry (rows,columns):' ...
	        '(Default [] = near square)'] ...
	        '-> Scalp map plotting options  (See >> help topoplot):' };
	inistr       = { fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))]) ...
	               ['ERP scalp maps' fastif(~isempty(EEG.setname), [' of ' EEG.setname ], '') ] ...
	               '' ['''electrodes'', ''off''' ] };
	result       = inputdlg2( txt, fastif( typeplot, 'ERP scalp map(s) -- pop_topoplot()', 'Component scalp map(s) -- pop_topoplot()'), 1,  inistr, 'topoplot');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	arg2   	     = eval( [ '[' result{1} ']' ] );
	topotitle    = result{2};
	rowcols     = eval( [ '[' result{3} ']' ] );
	options      = [ ',' result{4} ];
	if size(arg2(:),1) == 1, figure; try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end; end;
else
	options = [];
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else
			options = [ options ', [' num2str(varargin{i}) ']' ];
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
	maplimits = [ -max(maxlim, -minlim) max(maxlim, -minlim)];
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
        end;    
		subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1);
	end;
	set(gcf,'Position', [pos(1) pos(2) SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);

	if ~isnan(arg2(index))
		if typeplot
			if length( options ) < 2
				topoplot( SIGTMPAVG(:,index), EEG.chanlocs, 'maplimits', maplimits);
			else	
				eval( [ 'topoplot( SIGTMPAVG(:,index), EEG.chanlocs, ''maplimits'', maplimits' options ');' ] );
			end;
			if nbgraph == 1, title( topotitle );
			else title([int2str(arg2(index)) ' ms']);
			end;
		else
			if length( options ) < 2
			    if arg2(index) < 0
	    			topoplot( -EEG.icawinv(:, -arg2(index)), EEG.chanlocs );
	    		else
	    			topoplot( EEG.icawinv(:, arg2(index)), EEG.chanlocs );
	            end;    			
			else	
			    if arg2(index) < 0
				    eval( [ 'topoplot(  -EEG.icawinv(:, -arg2(index)), EEG.chanlocs' options ');' ]);
	            else
	    			eval( [ 'topoplot(  EEG.icawinv(:, arg2(index)), EEG.chanlocs' options ');' ]);
	            end;    			
			end;
			if nbgraph == 1, title( topotitle );
			else title(['' int2str(arg2(index))]);
			end;
		end;
		drawnow;
		axis square; 
		if index == size(arg2(:),1)
	        %pos = get(gca,'position');
	        %q = [pos(1) pos(2) 0 0];
	        %s = [pos(3) pos(4) pos(3) pos(4)];
	        %col = colormap;
	        %ax = subplot('position', [1 0 .05 1].*s+q);
	        cbar('vert');
	    end;	   
    else
    axis off
    end;
end;
if nbgraph> 1, a = textsc(0.5, 0.05, topotitle); set(a, 'fontweight', 'bold'); end;
if nbgraph== 1, com = 'figure;'; end;
axcopy(gcf, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); clear postmp;');

if length( options ) < 2
	com = [com sprintf('pop_topoplot(%s,%d,[%s], ''%s'', [%s]);', inputname(1), typeplot, sprintf('%d ',arg2), topotitle, int2str(rowcols) )];
else
	com = [com sprintf('pop_topoplot(%s,%d,[%s], ''%s'', [%s] %s);', inputname(1), typeplot, sprintf('%d ',arg2), topotitle, int2str(rowcols), options )];
end;
return;

		
