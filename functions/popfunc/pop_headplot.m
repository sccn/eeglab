% pop_headplot() - plot a spherically-splined EEG field map on 
%                  a semi-realistic 3-D head model.
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
%   others...  - All headplot options. See >> help headplot
%
% Output:
%   EEGOUT - EEG dataset, possibly with modified splinefile. 
%
% Note:
%   A new figure is created only when the pop_up window is called or when
%   several channels/components are plotted. Thus, you may call this command 
%   to draw 3-D topographic maps in a tiled figure.     
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
    error('Pop_headplot: cannot plot without channel locations. Use Edit/Dataset info');
end;

if nargin < 3
    if ~isfield(EEG, 'splinefile') | isempty(EEG.splinefile)
		 ButtonName=questdlg2( strvcat('3D head plot need to generate a spline file', ...
		                      'the first time it is called (and for every new channel', ...
		                      'location file). Do you want to generate this file now ?'), ...
		                      'Spline File', 'Cancel', 'Load file', 'Yes','Yes');
		 switch lower(ButtonName),
		      case 'cancel', return;
		      case 'load file', 
				    [filename, filepath] = uigetfile('*.spl', 'Load a spline file'); 
                    drawnow;
				    if filename == 0 return; end;
		            EEG.splinefile = [ filepath filename ];
		            options = [ ', ''load'',' EEG.splinefile ',' ]; 
                    com = 'EEG = ';
		      case 'yes',
				    [filename, filepath] = uiputfile('*.spl', 'Save spline file with .spl extension');
                    drawnow;
				    if filename == 0 return; end;
		            EEG.splinefile = [ filepath filename ];
		            headplot('setup', EEG.chanlocs, EEG.splinefile);
		            options = [ ', ''setup'',' EEG.splinefile ',' ]; 
                    com = 'EEG = ';
		 end;
		 return;
    else
    	options      = [ ',' ];
    end;
    
	if isempty(EEG.splinefile) | exist(EEG.splinefile) ~= 2
	    errmsg = sprintf('Pop_headplot: cannot find spline file %s. Check path. Aborting...',EEG.splinefile);
	    error(errmsg);
	end;

 	% which set to save
	% -----------------
	if typeplot
		%txt = sprintf('ERP head plot at these latencies (from %d to %d ms):\n(NaN -> empty subplot)(Ex: -100 NaN 100)', round(EEG.xmin*1000), round(EEG.xmax*1000));
		txt = sprintf('ERP head plot at these latencies (from %d to %d ms):', round(EEG.xmin*1000), round(EEG.xmax*1000));
	else
		%txt = ['Component numbers (negate index to invert component polarity):' 10 '(NaN -> empty subplot)(Ex: -1 NaN 3)'];
		txt = ['Component numbers (negate index to invert component polarity):' ];
	end;	
	txt = { txt ...
	        'Plot title:' ...
	        ['Plot geometry (rows,columns):' ...
	        '(Default [] = near square)'] ...
	        '-> headplot() options  (See >> help headplot):' };
	inistr       = { fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))]) ...
	               ['ERP scalp maps' fastif(~isempty(EEG.setname), [' of ' EEG.setname ], '') ] ...
	               '' '' };
	result       = inputdlg2( txt, fastif( typeplot, 'ERP head plot(s) -- pop_headplot()', 'Component head plot(s) -- pop_headplot()'), 1,  inistr, 'pop_headplot');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	arg2   	     = eval( [ '[' result{1} ']' ] );
	if length(arg2) > EEG.nbchan
		tmpbut = questdlg2(['This involve drawing ' int2str(length(arg2)) ' plots. Continue ?'], '', 'Cancel', 'Yes', 'Yes');
		if strcmp(tmpbut, 'Cancel'), return; end;
	end;
	topotitle    = result{2};
	rowcols     = eval( [ '[' result{3} ']' ] );
	options      = [ ',' result{4} ];
	if size(arg2(:),1) == 1, figure; end;
else
	% read or generate file if necessary
	% ----------------------------------
    loc = strmatch('load', varargin(1:2:end)); loc = loc*2-1;
    if ~isempty(loc)
        EEG.splinefile = varargin{ loc+1 };
        varargin(loc:loc+1) = [];
    end;
    loc = strmatch('setup', varargin(1:2:end)); loc = loc*2-1;
    if ~isempty(loc)
        headplot('setup', EEG.chanlocs, EEG.splinefile);
        varargin(loc:loc+1) = [];
    end;
    
	options = [];
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else
			options = [ options ', [' num2str(varargin{i}) ']' ];
		end;
	end;	
end;

if ~exist('topotitle')  
    topotitle = '';
end;    
if ~isfield(EEG, 'splinefile') | isempty(EEG.splinefile)
    error('Pop_headplot: cannot find spline file, aborting...');
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
			if length( options ) < 2
    			headplot( SIGTMPAVG(:,index), EEG.splinefile, 'maplimits', maplimits);
		    else	
			     eval( [ 'headplot( SIGTMPAVG(:,index), EEG.splinefile, ''maplimits'', maplimits ' options ');' ] );
			end;
			if nbgraph == 1, title( topotitle );
			else title([int2str(arg2(index)) ' ms']);
			end;
		else
			if length( options ) < 2
			    if arg2(index) < 0
			         headplot( -EEG.icawinv(:, -arg2(index)), EEG.splinefile);
		        else	
			         headplot( EEG.icawinv(:, arg2(index)), EEG.splinefile);
	            end;    			
			else	
			    if arg2(index) < 0
			         eval( [ 'headplot(  -EEG.icawinv(:, -arg2(index)), EEG.splinefile ' options ');' ] );
	            else
			         eval( [ 'headplot(  EEG.icawinv(:, arg2(index)), EEG.splinefile ' options ');' ] );
	            end;    			
			end;
			if nbgraph == 1, title( topotitle );
			else title(['' int2str(arg2(index))]);
			end;
		end;
		drawnow;
		axis square; 
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
    axcopy(gcf, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); rotate3d(gcf); clear postmp;');
end;
if nbgraph== 1, com = [ 'figure;' com ]; rotate3d(gcf); end;

if length( options ) < 2
	com = [com sprintf('pop_headplot(%s,%d,%s, ''%s'', [%s]);', inputname(1), typeplot, vararg2str(arg2), topotitle, int2str(rowcols) )];
else
	com = [com sprintf('pop_headplot(%s,%d,%s, ''%s'', [%s] %s);', inputname(1), typeplot, vararg2str(arg2), topotitle, int2str(rowcols), options )];
end;
return;
