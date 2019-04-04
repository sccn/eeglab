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
% 03-18-02 added title -ad & sm
% 03-30-02 added single trial capacities -ad

function com = pop_plottopo( EEG, channels, plottitle, singletrials, varargin);

com = '';
if nargin < 1
	help pop_plottopo;
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
	if length(result) == 0 return; end
    
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
    end
else 
	options ={ 'chanlocs' EEG.chanlocs 'frames' EEG.pnts 'limits' [EEG.xmin EEG.xmax 0 0]*1000 ...
               'title' plottitle 'chans' channels varargin{:}};
    addoptions = {};
end
% adapt frames to time limit.
if any(strcmp(addoptions,'limits'))
    addoptions{end+1} = 'frames';
    ilimits = find(strcmp(addoptions,'limits'))+1;
    timelims = addoptions{ilimits}(1:2);
    addoptions{end+1} = round(diff(timelims/1000)*EEG.srate);
end

try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
	
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
end

if ~isempty(addoptions)
	com = sprintf('figure; pop_plottopo(EEG, %s, ''%s'', %d, %s);', ...
				  vararg2str(channels), plottitle, singletrials, vararg2str(addoptions));
else
	com = sprintf('figure; pop_plottopo(EEG, %s, ''%s'', %d);', ...
				  vararg2str(channels), plottitle, singletrials);
end

return;
	
		
