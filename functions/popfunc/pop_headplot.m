% pop_headplot() - plot spherically-splined EEG field map on a semi-realistic 
%                  3-D head model.
% Usage:
%   >> EEGOUT = pop_headplot( EEG, typeplot, 
%                    latencies/components, title, rowscols, 'key', 'val' ...);
%
% Inputs:
%   EEG        - dataset structure
%   typeplot   - 1=channel, 0=component (default:1)
%   latencies/components  - array of latencies (in msec)
%                at which the head should be plotted.  
%                For components: array of index of components to plot.
%   title      - plot title.
%   rowscols   - Vector of the form [m,n] where m is total vertical tiles and n 
%                is horizontal tiles per page. If the number of maps exceeds m*n,
%                multiple figures will be produced {def|0 -> one near-square page}.
%
% Optional inputs:
%   'setup'    - ['setupfile.spl'] make the headplot spline file
%   'load'     - ['setupfile.spl'] load the headplot spline file
%   others...  - all headplot options. See >> help headplot
%
% Output:
%   EEGOUT - EEG dataset with potentially modified splinefile name.
%
% Note:
%   A new figure is created only when the pop_up window is called or when
%   several channel/components are plotted, so you may call this command 
%   to draw topographic maps in a tiled window.     
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 March 2002
%
% See also: eeglab(), headplot()

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
		 ButtonName=questdlg( ['3D head plot need to generate a spline file' 10 ...
		                      'the first time it is called (and for every new channel' 10 ...
		                      'location file. Do you want to generate this file now ?'], ...
		              'Spline file', 'Cancel', 'Load an existing file', 'Yes','Yes');
		 switch lower(ButtonName),
		      case 'cancel', return;
		      case 'load an existing file', 
				    [filename, filepath] = uigetfile('*.spl', 'Load a spline file');
				    if filename == 0 return; end;
		            EEG.splinefile = [ filepath filename ];
		            options = [ ', ''load'',' EEG.splinefile ',' ]; 
		      case 'yes',
				    [filename, filepath] = uiputfile('*.spl', 'Save spline file with .spl extension');
				    if filename == 0 return; end;
		            EEG.splinefile = [ filepath filename ];
		            headplot('setup', EEG.chanlocs, EEG.splinefile);
		            options = [ ', ''setup'',' EEG.splinefile ',' ]; 
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
		txt = sprintf('ERP head plot at these latencies (from %d to %d ms):\n(NaN -> empty subplot)(Ex: -100 NaN 100)', round(EEG.xmin*1000), round(EEG.xmax*1000));
	else
		txt = ['Component numbers (negate index to invert component polarity):' 10 '(NaN -> empty subplot)(Ex: -1 NaN 3)'];
	end;	
	txt = { txt ...
	        'Plot title:' ...
	        ['Plot geometry (rows,columns):' ...
	        '(Default [] = near square)'] ...
	        '-> headplot() options  (See >> help headplot):' };
	inistr       = { fastif( typeplot, '', ['1:' int2str(size(EEG.data,1))]) ...
	               ['ERP scalp maps' fastif(~isempty(EEG.setname), [' of ' EEG.setname ], '') ] ...
	               '' '' };
    help headplot;
	result       = inputdlg( txt, fastif( typeplot, 'ERP head plot(s) -- pop_headplot()', 'Component head plot(s) -- pop_headplot()'), 1,  inistr);
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	arg2   	     = eval( [ '[' result{1} ']' ] );
	topotitle    = result{2};
	rowcols     = eval( [ '[' result{3} ']' ] );
	options      = [ ',' result{4} ];
	if size(arg2(:),1) == 1, figure; end;
else
	% read or generate file if necessary
	% ----------------------------------
    loc = strmatch('load', varargin);
    if ~isempty(loc)
        EEG.splinefile = varargin{ loc+1 };
        varargin(loc:loc+1) = [];
    end;
    loc = strmatch('setup', varargin);
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
	        set(gcf,'Position', [pos(1) pos(2) SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
        end;    
		subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1);
	end;

	if ~isnan(arg2(index))
		if typeplot
			if length( options ) < 2
    			headplot( SIGTMPAVG(:,index), EEG.splinefile);
		    else	
			     eval( [ 'headplot( SIGTMPAVG(:,index), EEG.splinefile, '', ''' options ');' ] );
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
if nbgraph> 1, 
    a = textsc(0.5, 0.05, topotitle); 
    set(a, 'fontweight', 'bold');
    rotate3d(gcf);
    axcopy(gcf, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); rotate3d(gcf); clear postmp;');
end;
if nbgraph== 1, com = 'figure;'; end;

if length( options ) < 2
	com = [com sprintf('pop_headplot(%s,%d,[%s], ''%s'', [%s]);', inputname(1), typeplot, sprintf('%d ',arg2), topotitle, int2str(rowcols) )];
else
	com = [com sprintf('pop_headplot(%s,%d,[%s], ''%s'', [%s] %s);', inputname(1), typeplot, sprintf('%d ',arg2), topotitle, int2str(rowcols), options )];
end;
return;
