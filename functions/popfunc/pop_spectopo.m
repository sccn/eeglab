% pop_spectopo() - spectrum of all data channels plus scalp maps of power
%                  at specified frequencies. Calls spectopo(). 
%
% Usage:
%   >> pop_spectopo( EEG, dataflag); % pops-up
%   >> pop_spectopo( EEG, dataflag, timerange, 'key', 'val',...);
%
% Inputs:
%   EEG        - input dataset
%   dataflag   - 1 to process EEG data and 0 to process ICA activity.
%                Default is EEG data.
%   timerange  - epoch timerange (in msec) to use in computing the spectra
%
% Optional inputs:
%   'key','val' - optional spectopo() and topoplot() arguments 
%
% Outputs: Those of spectopo(). When nargin<2, a query window pops-up 
%          to ask for additional arguments and no output is returned.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 March 2002
%
% See also: spectopo(), topoplot()

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
% Revision 1.14  2002/07/30 18:13:23  arno
% same
%
% Revision 1.13  2002/07/30 18:11:04  arno
% debugging
%
% Revision 1.12  2002/07/29 22:25:58  arno
% updating message
%
% Revision 1.11  2002/07/29 22:15:50  arno
% reprogramming for spectopo component
%
% Revision 1.10  2002/07/28 21:23:00  arno
% adding message
%
% Revision 1.9  2002/07/28 21:00:36  arno
% debugging
%
% Revision 1.8  2002/07/26 01:59:53  arno
% debugging
%
% Revision 1.7  2002/07/24 18:16:49  arno
% changing default freqs
%
% Revision 1.6  2002/07/20 19:21:39  arno
% new spectopo compatibility
%
% Revision 1.5  2002/07/20 01:31:03  arno
% same
%
% Revision 1.4  2002/07/20 01:23:28  arno
% debuging varargin decoding
%
% Revision 1.3  2002/04/25 17:05:21  scott
% fruther editting -sm
%
% Revision 1.2  2002/04/25 17:00:47  scott
% editted help nd msgs, added % spectopo call to history -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 03-15-02 debugging -ad
% 03-16-02 add all topoplot options -ad
% 04-04-02 added outputs -ad & sm

function varargout = pop_spectopo( EEG, dataflag, timerange, varargin);

varargout{1} = '';
if nargin < 1
	help pop_spectopo;
	return;
end;	

if nargin < 2
	dataflag = 1;
end;

if nargin < 3
	if dataflag
		geometry = { [2 1] [2 1] [2 1] [2 1] };
		promptstr    = { { 'style' 'text' 'string' 'Epoch time range (ms) to include [min max]:' }, ...
						 { 'style' 'edit' 'string' [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)] }, ...
						 { 'style' 'text' 'string' 'Percent windows to sample (0 to 100):'}, ...
						 { 'style' 'edit' 'string' '20' }, ...
						 { 'style' 'text' 'string' 'Scalp map frequencies (Hz):'}, ...
						 { 'style' 'edit' 'string' '10 20 30' }, ...
						 { 'style' 'text' 'string' 'Scalp map plotting options: See >> help topoplot' } ...
						 { 'style' 'edit' 'string' '''electrodes'',''off''' } };
		result       = inputgui( geometry, promptstr, 'pophelp(''spectopo'')', 'Channel spectra and maps -- pop_spectopo()');
		if size(result,1) == 0 return; end;
		timerange    = eval( [ '[' result{1} ']' ] );
		options = [];
		if eval(result{2}) ~= 100, options = [ options ', ''percent'', '  result{2} ]; end;
		if ~isempty(result{3})   , options = [ options ', ''freq'', ['  result{3} ']' ]; end;
		if ~isempty(result{4}),    options = [ options ',' result{4} ]; end;
	else
		if isempty(EEG.chanlocs)
			error('Pop_spectopo: cannot plot component contribution without channel locations');
		end;
		geometry = { [2 1] [2 1] [2 1] [2 1] [2 1] [2 1] [2 1] [2 0.18 0.78] [2 1]  };
		promptstr    = { { 'style' 'text' 'string' 'Epoch time range (ms) to include [min max]:' }, ...
						 { 'style' 'edit' 'string' [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)] }, ...
						 { 'style' 'text' 'string' 'Percent windows to sample (0 to 100):'}, ...
						 { 'style' 'edit' 'string' '20' }, ...
						 { 'style' 'text' 'string' 'Scalp map frequency (Hz):'}, ...
						 { 'style' 'edit' 'string' '10' }, ...
						 { 'style' 'text' 'string' 'Channel number for comp. contributions:', 'tooltipstring', ...
						 ['The component contribution will be plotted at that channel' 10 ...
						 'If no channel is given, the contribution of component at' 10 ...
						 'the channel of max. power is returned. If 0 is entered, then' ...
						 'spectopo finds components that globally contribute the most at one frequency'] }, ...
						 { 'style' 'edit' 'string' '' }, ...
						 { 'style' 'text' 'string' 'Component indices to compute spectra of:' }, ...
						 { 'style' 'edit' 'string' ['1:' int2str(size(EEG.icaweights,1))] }, ...
						 { 'style' 'text' 'string' 'Number of component''s scalp map of  to plot:' }, ...
						 { 'style' 'edit' 'string' '' }, ...
						 { 'style' 'text' 'string' 'Component indices for plotting scalp map:', 'tooltipstring', ...
						 ['By default the map of the components of maximum power at the selected' 10 ...
						 'channel (or power at channel of max. power) are plotted' ] }, ...
						 { 'style' 'edit' 'string' '' }, ...
						 { 'style' 'text' 'string' 'Compute component or (data-component) spectra:', 'tooltipstring' ...
						 ['Either compute the spectra of the components'' activity (set)' 10 ...
						 'or the spectra of the data minus every single slected component (unset)']}, ...
						 { 'style' 'checkbox' 'value' 1 } { }, ...
						 { 'style' 'text' 'string' 'Scalp map plotting options: See >> help topoplot' } ...
						 { 'style' 'edit' 'string' '''electrodes'',''off''' } };
		result       = inputgui( geometry, promptstr, 'pophelp(''spectopo'')', 'Component spectra and maps -- pop_spectopo()');
		if size(result,1) == 0 return; end;
		timerange    = eval( [ '[' result{1} ']' ] );
		if eval(result{2}) ~= 100, options = [ options ', ''percent'', '  result{2} ]; end;
		if ~isempty(result{3})   , options = [ options ', ''freq'', ['  result{3} ']' ]; end;
		if ~isempty(result{4})   , options = [ options ', ''plotchan'', ' result{4} ]; end;
		if ~isempty(result{5})   , options = [ options ', ''icacomps'', [' result{5} ']' ]; end;
		if ~isempty(result{6})   , options = [ options ', ''nicamaps'', ' result{6} ]; end;
		if ~isempty(result{7})   , options = [ options ', ''icamaps'', [' result{7} ']' ]; end;
		if ~result{8}, options = [ options ', ''icamode'', ''sub''' ]; end;
		if ~isempty(result{9}), options      =  [ options ',' result{9} ]; end;
	end;		
	figure;
else
	options = [',' vararg2str(varargin)];
	if isempty(timerange)
		timerange = [ EEG.xmin*1000 EEG.xmax*1000 ];
	end;
	if nargin < 3 
		percent = 100;
	end;
	if nargin < 4 
		topofreqs = [];
	end;
end;

if ~isempty(EEG.chanlocs)
	spectopooptions = [ options ', ''chanlocs'', EEG.chanlocs' ];
	if dataflag == 0 % i.e. components
		spectopooptions = [ spectopooptions ', ''weights'', EEG.icaweights*EEG.icasphere' ];
	end;
else
	spectopooptions = options;
end;

% The programming here is a bit redundant but it tries to optimize 
% memory usage.
% ----------------------------------------------------------------
if timerange(1)/1000~=EEG.xmin | timerange(2)/1000~=EEG.xmax
	posi = round( (timerange(1)/1000-EEG.xmin)*EEG.srate )+1;
	posf = round( (timerange(2)/1000-EEG.xmin)*EEG.srate )+1;
	pointrange = posi:posf;
	if posi == posf, error('pop_spectopo: empty time range'); end;
	fprintf('pop_spectopo: slecting time range %6.2f ms to %6.2f ms (points %d to %d)\n', ...
			timerange(1), timerange(2), posi, posf);
end;
if exist('pointrange') == 1, SIGTMP = EEG.data(:,pointrange,:); totsiz = length( pointrange);
else                         SIGTMP = EEG.data; totsiz = EEG.pnts;
end;
SIGTMP = reshape(SIGTMP, size(SIGTMP,1), size(SIGTMP,2)*size(SIGTMP,3));

% put warning if boundaries in continuous data
% --------------------------------------------

% outputs
% -------
outstr = '';
if nargin >= 2
	for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
	if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the data and generate output and history commands
% ------------------------------------------------------
popcom = sprintf('figure; pop_spectopo(%s, %d, [%s] %s);', inputname(1), dataflag, num2str(timerange), options);
com = sprintf('%s spectopo( SIGTMP, totsiz, EEG.srate %s);', outstr, spectopooptions);
eval(com)
if nargout < 2 & nargin < 3
	varargout{1} = popcom;
end;

return;

		
