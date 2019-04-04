% pop_timtopo() - call the timtopo() function for epoched EEG datasets. 
%                 Plots the epoch mean for each channel on a single axis,
%                 plus scalp maps of the data at specified latencies.
% Usage:
%   >> pop_timtopo( EEG, timerange, topotimes, title, 'key', 'val', ...);
%
% Inputs:
%   EEG         - input dataset
%   timerange   - [min max] epoch time range (in ms) to plot 
%   topotimes   - array of times to plot scalp maps {Default: NaN 
%                 = display scalp map at frame of max var()}
%
% Optional inputs:
%   title       - optional plot title
%   'key','val' - optional topoplot() arguments (see >> help topoplot)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: timtopo()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% 01-25-02 reformated help & license -ad 
% 02-16-02 text interface editing -sm & ad 
% 03-15-02 add all topoplot options -ad
% 03-18-02 added title -ad & sm

function com = pop_timtopo( EEG, timerange, topotime, plottitle, varargin);

com = '';
if nargin < 1
	help pop_timtopo;
	return;
end;	

if nargin < 3
	promptstr    = { 'Plotting time range (ms):', ...
			         ['Scalp map latencies (ms, NaN -> max-RMS)'], ...
					 'Plot title:' ...
			         'Scalp map options (see >> help topoplot):' };
	inistr       = { [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)], ...
			         'NaN', ...
	                 ['ERP data and scalp maps' fastif(~isempty(EEG.setname), [' of ' EEG.setname ], '') ], ...
			         ''  };
	result       = inputdlg2( promptstr, 'ERP data and scalp maps -- pop_timtopo()', 1, inistr, 'pop_timtopo');
	if size(result,1) == 0 return; end
	timerange    = eval( [ '[' result{1} ']' ] );
	topotime     = eval( [ '[' result{2} ']' ] );
	plottitle    = result{3};
	options      = [ ',' result{4} ];
	figure;
else
	options = [];
	for i=1:length( varargin )
		if ischar( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else
			options = [ options ', [' num2str(varargin{i}) ']' ];
		end
	end;	
end
try, icadefs; set(gcf, 'color', BACKCOLOR, 'Name', ' timtopo()'); catch, end

if exist('plottitle') ~= 1
    plottitle = ['ERP data and scalp maps' fastif(~isempty(EEG.setname), [' of ' EEG.setname ], '') ];
end
    
if ~isempty(EEG.chanlocs)
    if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end
	SIGTMP = reshape(EEG.data, size(EEG.data,1), EEG.pnts, EEG.trials);
	posi = round( (timerange(1)/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
	posf = round( (timerange(2)/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
	if length( options ) < 2
    	timtopo( mean(SIGTMP(:,posi:posf,:),3), EEG.chanlocs, 'limits', [timerange(1) timerange(2) 0 0], 'plottimes', topotime, 'chaninfo', EEG.chaninfo);
        com = sprintf('figure; pop_timtopo(EEG, [%s], [%s], ''%s'');', num2str(timerange), num2str(topotime), plottitle);
	else
		com = sprintf('timtopo( mean(SIGTMP(:,posi:posf,:),3), EEG.chanlocs, ''limits'', [timerange(1) timerange(2) 0 0], ''plottimes'', topotime, ''chaninfo'', EEG.chaninfo %s);', options);
		eval(com)
	    com = sprintf('figure; pop_timtopo(EEG, [%s], [%s], ''%s'' %s);', num2str(timerange), num2str(topotime), plottitle, options);
	end;		
else
	fprintf('Cannot make plot without channel locations\n');
	return;
end
return;

		
