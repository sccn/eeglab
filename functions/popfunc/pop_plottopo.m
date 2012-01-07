% pop_plottopo() - plot one or more concatenated multichannel data epochs 
%                  in a topographic array format using plottopo()
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
%   others...    - additional plottopo arguments {'axsize', 'color', 'ydir'
%                  'vert'} (see >> help plottopo)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 March 2002
%
% See also: plottopo()

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
	uilist    = { { 'style' 'text' 'string' 'Channels to plot' } ...
                  { 'style' 'edit' 'string' [ '1:' num2str( EEG.nbchan ) ] 'tag' 'chan' } ...
                  { 'style' 'text' 'string' 'Plot title' } ...
                  { 'style' 'edit' 'string' fastif(isempty(EEG.setname), '',EEG.setname)  'tag' 'title' } ...                  
                  { 'style' 'text' 'string' 'Plot single trials' } ...
                  { 'style' 'checkbox' 'string' '(set=yes)' 'tag' 'cbst' } ...
                  { 'style' 'text' 'string' 'Plot in rect. array' } ...
                  { 'style' 'checkbox' 'string' '(set=yes)' 'tag' 'cbra' } ...
                  { 'style' 'text' 'string' 'Other plot options (see help)' } ...
                  { 'style' 'edit' 'string' '''ydir'', 1' 'tag' 'opt' } };
    geometry = { [1 1] [1 1] [1 1] [1 1] [1 1] };
    [result userdata tmphalt restag ] = inputgui( 'uilist', uilist, 'geometry', geometry, 'helpcom', 'pophelp(''pop_plottopo'')', 'title', 'Topographic ERP plot - pop_plottopo()');
	if length(result) == 0 return; end;
    
	channels     = eval( [ '[' restag.chan ']' ] );
	plottitle    = restag.title;
	singletrials = restag.cbst;
    addoptions   = eval( [ '{' restag.opt '}' ] );
    rect         = restag.cbra;
    
    figure('name', ' plottopo()');
    options ={ 'frames' EEG.pnts 'limits' [EEG.xmin EEG.xmax 0 0]*1000 ...
               'title' plottitle 'chans' channels addoptions{:} };
    if ~rect
        options = { options{:} 'chanlocs' EEG.chanlocs };
    end;
else 
	options ={ 'chanlocs' EEG.chanlocs 'frames' EEG.pnts 'limits' [EEG.xmin EEG.xmax 0 0]*1000 ...
               'title' plottitle 'chans' channels varargin{:}};
    addoptions = {};
end;
% adapt frames to time limit.
if any(strcmp(addoptions,'limits'))
    addoptions{end+1} = 'frames';
    ilimits = find(strcmp(addoptions,'limits'))+1;
    timelims = addoptions{ilimits}(1:2);
    addoptions{end+1} = round(diff(timelims/1000)*EEG.srate);
end

try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
	
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

if ~isempty(addoptions)
	com = sprintf('figure; pop_plottopo(%s, %s, ''%s'', %d, %s);', ...
				  inputname(1), vararg2str(channels), plottitle, singletrials, vararg2str(addoptions));
else
	com = sprintf('figure; pop_plottopo(%s, %s, ''%s'', %d);', ...
				  inputname(1), vararg2str(channels), plottitle, singletrials);
end;

return;
	
		
