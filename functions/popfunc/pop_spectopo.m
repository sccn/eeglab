% pop_spectopo() - spectrum of all data channels plus scalp maps of power
%                  at specified frequencies. Calls spectopo(). 
%
% Usage:
%   >> pop_spectopo( EEG, timerange, percent, topofreqs, 'key', 'val',...);
%
% Inputs:
%   EEG        - input dataset
%   timerange  - epoch timerange (in msec) to use in computing the spectra
%   percent    - percentage of possible data windows to sample. Use <<100
%                to speed the computation. Range: 0 to 100. {Default: 100}
%   topofreqs  - array of frequencies tgo plot scalp maps showing power
%
% Optional inputs:
%   'key','val' - optional topoplot() arguments (see topoplot())
%
% Outputs: Those of spectopo(). When nargin<4, a query window pops-up 
%          to ask for additional arguments and no output is returned.
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

function varargout = pop_spectopo( EEG, timerange, percent, topofreqs, varargin);

varargout{1} = '';
if nargin < 1
	help pop_spectopo;
	return;
end;	

if nargin < 2
	promptstr    = { 'Epoch time range (ms) to include [min max]:', ...
	                 'Percent windows to sample (0 to 100):', ...
			         'Scalp map frequencies (Hz):', ...
			         'Scalp map plotting options: See >> help topoplot' };
	inistr       = { [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)], ...
			         '20' '10' '''electrodes'',''off''' };
    help topoplot;
	result       = inputdlg( promptstr, 'Channel spectra and maps -- pop_spectopo()', 1, inistr);
	if size(result,1) == 0 return; end;
	timerange    = eval( [ '[' result{1} ']' ] );
	percent      = eval( [ '[' result{2} ']' ] ); 
	topofreqs    = eval( [ '[' result{3} ']' ] );
	options      =  [ ',' result{4} ];
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

%if isempty(topofreqs)
%    close(gcf);
%    error('Pop_spectopo: you must enter at least one frequency for scalp map plotting.');
%end;    

if ~isempty(EEG.chanlocs)
    % The programming here is a bit redundant but it tries to optimize 
    % memory use.
    % ----------------------------------------------------------------
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

	% plot the data and generate output and history commands
	% ------------------------------------------------------
	if length( options ) < 2
	    options = '';
	end;
	popcom = sprintf('figure; pop_spectopo(%s, [%s], %s, [%s] %s);', inputname(1), num2str(timerange), num2str(percent), num2str(topofreqs), options);
	com = sprintf('%s spectopo( SIGTMP, totsiz, EEG.srate, ''freq'', topofreqs, ''chanlocs'', EEG.chanlocs, ''limits'', [nan nan nan nan nan nan], ''percent'', percent %s);', outstr, options);
	eval(com)
	varargout{1} = [10 popcom 10 '% ' com];

	%title(['Spectrum head plots (time range ' num2str(timerange(1)) '-' num2str(timerange(2)) ')' ]);
else
	fprintf('Cannot plot scalp maps without channel locations...\n');
	return;
end;

return;

		
