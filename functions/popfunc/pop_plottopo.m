% pop_plottopo() - plot concatenated multichannel data epochs in a
%                  topographic array format using plottopo()
% Usage:
%   >> pop_plottopo( EEG ); % pop-up
%   >> pop_plottopo( EEG, channels );
%   >> pop_plottopo( EEG, channels, title, singletrials);
%   >> pop_plottopo( EEG, channels, title, singletrials, axsize, ...
%                         'color', ydir, vert);
%
% Inputs:
%   EEG          - input dataset
%   channels     - indices of channels to plot
%   title        - plot title. Default is none.
%   singletrials - [0|1], 0 plot average, 1 plot individual
%                  single trials. Default is 0.
%   others...    - plottopo additional arguments 'axsize', 'color', 'ydir'
%                  'vert' } (see >> help plottopo)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 March 2002
%
% See also: plottopo()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 10 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.2  2002/04/18 15:45:11  scott
% editted msgs -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 02-16-02 text interface editing -sm & ad 
% 03-18-02 added title -ad & sm
% 03-30-02 added single trial capacities -ad

function com = pop_plottopo( EEG, channels, plottitle, singletrials, varargin);

com = '';
if nargin < 1
	help pop_plottopo;
	return;
end;	

if isempty(EEG.chanlocs)
	fprintf('Cannot plot without knowing channel locations. Use Edit/Dataset info\n');
	return;
end;

if nargin < 2
	promptstr    = { 'Channels to plot:' ...
					 'Plot title:' ...
					 'Plot single trials instead of average (yes|no)' 'Vertical lines (ms)' };
	inistr       = { [ '1:' num2str( EEG.nbchan ) ] ...
					 ['ERP in scalp order' fastif(isempty(EEG.setname), '',[' of ' EEG.setname])] ...
					 'no' ''};
	result       = inputdlg( promptstr, 'Topographic ERP plot - pop_plottopo()', 1, inistr);
	if size(result,1) == 0 return; end;
	channels     = eval( [ '[' result{1} ']' ] );
	plottitle    = result{2};
	singletrials = strcmp( lower(result{3}), 'yes');
	vert         = eval( [ '[' result{4} ']' ] );
    figure;
	if ~isempty(vert)
		options ={ EEG.chanlocs, EEG.pnts, [EEG.xmin EEG.xmax 0 0]*1000, plottitle, channels, [.07 .07],0,1,vert };
	else 
		options ={ EEG.chanlocs, EEG.pnts, [EEG.xmin EEG.xmax 0 0]*1000, plottitle, channels };
	end;
else 
	options ={ EEG.chanlocs, EEG.pnts, [EEG.xmin EEG.xmax 0 0]*1000, plottitle, channels varargin{:}};
end;
	
if exist('plottitle') ~= 1
    plottitle = '';
end;    
if exist('singletrials') ~= 1
    singletrials = 0;
end;    

if singletrials 
    plottopo( EEG.data, options{:} );
else
    plottopo( mean(EEG.data,3), options{:} );
end;
set(gcf, 'name', 'ERP in scalp order -- pop_plottopo()');

if ~isempty(options(6:end))
	com = sprintf('figure; pop_plottopo(%s, [%s], ''%s'', %d %s);', ...
				  inputname(1), int2str(channels), plottitle, singletrials, vararg2str(options(6:end)));
else
	com = sprintf('figure; pop_plottopo(%s, [%s], ''%s'', %d);', ...
				  inputname(1), int2str(channels), plottitle, singletrials);
end;

return;
	
		
