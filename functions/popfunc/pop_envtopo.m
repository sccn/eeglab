% pop_envtopo() - Plot envelope of an averaged EEG epoch, plus scalp maps 
%                 of largest or specified components reference to their 
%                 time point maximum amplitude. Calls envtopo().
% Usage:
%   >> pop_envtopo( EEG, 'key', 'val', ...);
%
% Inputs:
%   EEG        - input dataset
%   timerange  - [min max] time range (in msec) to plot 
%
% Optional inputs:
%   'key','val' - optional spectopo() and topoplot() arguments 
%                 (see >> help topoplot())
%
% Outputs: Same as envtopo(). When nargin < 3, a query window pops-up 
%          to ask for additional arguments and no outputs are returned.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: envtopo(), eeglab()

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
% Revision 1.11  2002/10/08 15:55:04  arno
% nothing
%
% Revision 1.10  2002/10/05 01:52:51  arno
% adapt to new syntax (but still using old convension)
%
% Revision 1.9  2002/10/05 01:07:12  arno
% debugging function call
%
% Revision 1.8  2002/08/12 16:27:38  arno
% inputdlg2
%
% Revision 1.7  2002/08/12 01:32:37  arno
% color
%
% Revision 1.6  2002/08/11 22:09:36  arno
% color
%
% Revision 1.5  2002/08/11 20:51:32  arno
% color
%
% Revision 1.4  2002/04/25 17:46:33  scott
% added 'electrodes','off' default -sm
%
% Revision 1.3  2002/04/25 17:40:30  scott
% edited help msg and added com to history -sm
%
% Revision 1.2  2002/04/18 15:53:26  scott
% editted msgs -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-16-02 add all topoplot options -ad
% 03-18-02 added title -ad & sm

function varargout = pop_envtopo( EEG, timerange, varargin);

varargout{1} = '';
if nargin < 1
	help pop_envtopo;
	return;
end;	
if isempty( EEG.icasphere )
	disp('Error: cannot plot without ICA weights. Use Tools/Run ICA.'); return;
end;
if isempty(EEG.chanlocs)
	fprintf('Cannot plot without channel locations. Use Edit/Dataset info.\n');
	return;
end;
if exist('envtitle') ~= 1
	envtitle = 'Largest ERP components';
end;

if nargin < 3
	% which set to save
	% -----------------
	promptstr    = { 'Enter plotting time range (msec):', ...
					 'Enter time range for component contribution (msec):', ...
					 'Plot this many largest components (1-7):', ...
					 'ELSE plot these components only (<=7) (Ex: 2:4,7 9:11):', ...
                     'Indices of component to subtract from data before plotting' ...
					 'Plot title:' ...
			         'Additional spectopo() and topoplot() options:' };
	inistr       = { [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)], ...
	                 [num2str( EEG.xmin*1000) ' ' num2str(EEG.xmax*1000)], ...
                     '7', ...
	                 '', ...
                     '', ...
	                 ['Largest ERP components' fastif(isempty(EEG.setname), '',[' of ' EEG.setname])] ...
	                 '''electrodes'',''off''' };
	result       = inputdlg2( promptstr, 'Components and ERP envelope -- pop_envtopo()', 1, inistr, 'pop_envtopo');
	if length(result) == 0 return; end;

    timerange    = eval( [ '[' result{1} ']' ] );
    options = ',';
	if ~isempty( result{2} ), options = [ options '''limcontrib'',  [' result{2} '],' ]; end;
	if ~isempty( result{4} ), options = [ options '''compnums'',  [' result{4} '],' ]; 
    else                      options = [ options '''compnums'',   -' result{3} ',' ];
    end;
    if ~isempty(result{5}),   options = [ options '''subcomps'',   [' result{5} '],' ]; end;
    if ~isempty(result{6}),   options = [ options '''title'', ''' result{6} ''',' ]; end;
	options      =  [ options result{7} ];
	figure;
    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
else
	options = [',' vararg2str( varargin ) ];
end;
    
sigtmp = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
posi = round( (timerange(1)/1000-EEG.xmin) * EEG.srate) + 1;
posf = round( (timerange(2)/1000-EEG.xmin) * EEG.srate) + 1;

% outputs
% -------
outstr = '';
if nargin >= 4
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the datas and generate output command
% --------------------------------------------
if length( options ) < 2, options = ''; end;
varargout{1} = sprintf('figure; pop_envtopo(%s, [%s] %s);', inputname(1), num2str(timerange), ...
                       options);
com =  sprintf(['%s envtopo(mean(sigtmp(:,posi:posf,:),3), EEG.icaweights*EEG.icasphere, ' ...
                '''chanlocs'', EEG.chanlocs, ''icawinv'', EEG.icawinv, ''limits'',' ...
                '[timerange(1) timerange(2) 0 0] %s);' ] , outstr, options);
eval(com);

return;

		
