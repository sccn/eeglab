% pop_topoplot() - Plot scalp map(s) in a figure window. If number of 
%                  arguments is 2, pop-up an interactive query window.
%
% Usage:
%   >> pop_topoplot( EEG, typeplot, latencies, title, options...);
%
% Inputs:
%   EEG        - Input dataset (see eeglab())
%   typeplot   - 1-> plot ERP maps, 0-> plot components {Default:1}
%   items      - [array] If typoplot ERP maps, epoch latencies (in ms) 
%                to plot the scalp maps. If typeplot components, 
%                component indices to plot. In this case,
%                negative indices -> invert map polarity; 
%                NaN -> leave a blank subplot. (Ex: [1 -3 NaN 4])
%   title      - plot title.
%   rowscols   - Vector of the form [m,n] giving [rows, cols] per page.
%                If the number of maps exceeds m*n, multiple figures 
%                are produced {Default|0 -> one near-square page}.
%   options    - topoplot() argument options. Separate using commas. 
%                Example 'style', 'straight'. See topoplot() help 
%                for further details. {Default: none}. 
%
% Note:
%   A new figure is created automatically only when the pop_up window is 
%   called or when more than one page of maps are plotted. Thus, this 
%   command may be used to draw topographic maps in a figure sub-axis.
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
% Revision 1.11  2002/08/17 22:26:05  scott
% typo
%
% Revision 1.10  2002/08/17 22:24:15  scott
% tyop
%
% Revision 1.9  2002/08/17 22:21:02  scott
% help msg and menu text
%
% Revision 1.8  2002/08/17 22:16:19  scott
% editing help msg
%
% Revision 1.7  2002/08/17 22:15:51  scott
% *** empty log message ***
%
% Revision 1.6  2002/08/13 17:49:37  arno
% debug color
%
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
		txt = sprintf(['Plotting ERP scalp maps at these latencies:' 10 ' (range: %d to %d ms, NaN -> empty):'], round(EEG.xmin*1000), round(EEG.xmax*1000));
	else
		%txt = ['Component numbers (negate index to invert component polarity):' 10 '(NaN -> empty subplot)(Ex: -1 NaN 3)'];
		txt = sprintf('Plotting scalp maps at these latencies (from %d to %d ms):', round(EEG.xmin*1000), round(EEG.xmax*1000));
	end;	
	txt = { txt ...
	        'Plot title:' ...
	        ['Plot geometry (rows,columns):' ...
	        '(Default [] -> near square)'] ...
	        '-> Scalp map plotting options (see >> help topoplot):' };
        if typeplot
	  inistr       = { fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))]) ...
	                 [fastif(~isempty(EEG.setname), [EEG.setname ' ERP'], '') ] ...
	                 '' ['''electrodes'', ''off''' ] };
        else
	  inistr       = { fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))]) ...
	                 [fastif(~isempty(EEG.setname), [EEG.setname ' components'], '') ] ...
	                 '' ['''electrodes'', ''off''' ] };
        end
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

		
