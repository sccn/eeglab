% pop_spectopo() - spectrum of selected channels with 
%                  spectrum head plots at specific frequencies.
%
% Usage:
%   >> pop_spectopo( EEG, timerange, percent, topofreqs, 'key', 'val', ...);
%
% Inputs:
%   EEG        - input dataset
%   timerange  - timerange in millisecond to plot the enveloppe
%   percent    - percentage of the data to sample. Allow to speed up
%                the computation. From 1 to 0.01. Default is 1.
%   topofreqs  - array of frequency(ies) for plotting the head
%
% Optional inputs:
%   'key','val' - optional topoplot() arguments (see topoplot())
%
% Outputs: same as spectopo(), no outputs are returned when a
%          window pops-up to ask for additional arguments
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 March 2002
%
% See also: spectopo()

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

% 03-15-02 debuging -ad
% 03-16-02 add all topoplot options -ad
% 04-04-02 added outputs -ad & sm

function varargout = pop_spectopo( EEG, timerange, percent, topofreqs, varargin);

varargout{1} = '';
if nargin < 1
	help pop_spectopo;
	return;
end;	

if nargin < 4
	promptstr    = { 'Time range to plot [min max] (ms):', ...
	                 'Percent of data to sample (0 to 100, Ex: 50):', ...
			         'Frequencies to plot scalp maps (Hz):', ...
			         'Scalp map text options: See >> help topoplot' };
	inistr       = { [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)], ...
			         '100' '' '' };
    help topoplot;
	result       = inputdlg( promptstr, 'Channel spectra and maps -- pop_spectopo()', 1, inistr);
	if size(result,1) == 0 return; end;
	timerange    = eval( [ '[' result{1} ']' ] );
	percent      = eval( [ '[' result{2} ']' ] ) / 100; 
	topofreqs    = eval( [ '[' result{3} ']' ] );
	options      =  [ ',' result{4} ];
	figure;
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

if isempty(topofreqs)
    close(gcf);
    error('Pop_spectopo: you must enter some frequencies for headplots');
end;    

if ~isempty(EEG.chanlocs)
    % the programming here is a bit redundant but it tries to optimize memory reservation
    % -----------------------------------------------------------------------------------
    if timerange(1)/1000~=EEG.xmin | timerange(2)/1000~=EEG.xmax
	   posi = round( (timerange(1)/1000-EEG.xmin)*EEG.srate )+1;
	   posf = round( (timerange(2)/1000-EEG.xmin)*EEG.srate )+1;
	   pointrange = posi:posf;
	end;
    if exist('pointrange') == 1, SIGTMP = EEG.data(:,pointrange,:); totsiz = length( pointrange);
    else                         SIGTMP = EEG.data; totsiz = EEG.pnts;
    end;
	SIGTMP = reshape(SIGTMP, size(SIGTMP,1), size(SIGTMP,2)*size(SIGTMP,3));

	% outputs
	% -------
	outstr = '';
	if nargin >= 4
	    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
	    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
	end;

	% plot the datas and generate output command
	% --------------------------------------------
	if length( options ) < 2
	    options = '';
	end;
	varargout{1} = sprintf('figure; pop_spectopo(%s, [%s], %s, [%s] %s);', inputname(1), num2str(timerange), num2str(percent), num2str(topofreqs), options);
	com = sprintf('%s spectopo( SIGTMP, totsiz, EEG.srate, topofreqs, EEG.chanlocs, [nan nan nan nan nan nan], '''', 4, percent %s);', outstr, options);
	eval(com)
	%title(['Spectrum head plots (time range ' num2str(timerange(1)) '-' num2str(timerange(2)) ')' ]);
else
	fprintf('Can not plot witout channel location\n');
	return;
end;

return;

		
