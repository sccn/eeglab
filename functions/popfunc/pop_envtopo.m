% pop_envtopo() - Plot envelope of an averaged EEG epoch, plus scalp maps 
%                 of specified or largest contributing components referenced 
%                 to their time point of maximum variance in the epoch or specified
%                 sub-epoch. Calls envtopo(). When nargin < 3, a query window 
%                 pops-up to allow additional arguments.
% Usage:
%   >> pop_envtopo( EEG ); % pop-up window mode
%   >> pop_envtopo( EEG, timerange, 'key', 'val', ...);
%
% Inputs:
%   EEG        - input dataset. Can also be an array of two epoched datasets. 
%                In this case, the epoch mean (ERP) of the second is subtracted 
%                from the epoch mean (ERP) of the first. Note: The ICA weights 
%                must be the same for the two datasets.
%   timerange  - [min max] time range (in ms) in epoch to plot, or if [], from EEG
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
% Revision 1.29  2005/03/07 21:26:02  arno
% chaninfo
%
% Revision 1.28  2005/03/05 00:13:19  arno
% adding chaninfo
%
% Revision 1.27  2004/11/12 07:17:47  scott
% changed entopo
%
% changed envtopo 'times' to 'timerange'
%
% Revision 1.26  2004/11/11 15:06:04  scott
% made default timerange from EEG.xmax, EEG.xmin
%
% Revision 1.25  2004/11/11 14:37:00  scott
% changed 'limits' option to envtopo() to new 'times'
%
% Revision 1.24  2004/05/08 01:47:56  scott
% line 143 EEG.pnts -> EEG(1).pnts
%
% Revision 1.23  2004/04/25 16:28:09  scott
% edit help message
%
% Revision 1.22  2004/03/03 19:29:49  arno
% remove pvaf on since it is the default
%
% Revision 1.21  2004/03/03 19:11:47  arno
% fixed pvaf on
%
% Revision 1.20  2003/08/18 21:22:28  scott
% made option 'pvaf,'on' the default
%
% Revision 1.19  2003/05/12 23:47:44  arno
% removing extra blanks
%
% Revision 1.18  2003/05/12 23:37:28  arno
% verbose off
%
% Revision 1.17  2003/04/15 16:55:43  arno
% allowing to plot 20 components
%
% Revision 1.16  2003/03/14 03:22:03  arno
% documenting dataset subtraction
%
% Revision 1.15  2003/03/14 03:17:55  arno
% allowing to subtract 2 datasets
%
% Revision 1.14  2002/11/12 16:25:13  scott
% warning command edit
%
% Revision 1.13  2002/10/17 02:35:46  arno
% handles nan now
%
% Revision 1.12  2002/10/09 22:00:53  arno
% upodating for new envtopo
%
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

% 01-25-02 reformatted help & license -ad 
% 03-16-02 added all topoplot() options -ad
% 03-18-02 added title -ad & sm

function varargout = pop_envtopo( EEG, timerange, varargin);

varargout{1} = '';
if nargin < 1
	help pop_envtopo;
	return;
end;	
if length(EEG) == 1 & isempty( EEG.icasphere )
	disp('Error: cannot make plot without ICA weights. See "Tools > Run ICA".'); return;
end;
if length(EEG) == 1 & isempty(EEG.chanlocs)
	fprintf('Cannot make plot without channel locations. See "Edit > Dataset info".\n');
	return;
end;
if exist('envtitle') ~= 1
	envtitle = 'Largest ERP components';
end;

options = ',';
if nargin < 3
	% which set to save
	% -----------------
	promptstr    = { 'Enter time range (in ms) to plot:', ...
			 'Enter time range (in ms) to rank component contributions:', ...
			 'Number of largest contributing components to plot (1-20):', ...
			 'Else plot these component numbers only (<21) (Ex: 2:4,7):', ...
                         'Component numbers to remove from data before plotting:' ...
			 'Plot title:' ...
			 'Optional topoplot() and spectopo() arguments:' };
	inistr       = { [num2str( EEG(end).xmin*1000) ' ' num2str(EEG(end).xmax*1000)], ...
	                 [num2str( EEG(end).xmin*1000) ' ' num2str(EEG(end).xmax*1000)], ...
                     '7', ...
	                 '', ...
                     '', ...
	                 ['Largest ERP components' fastif(isempty(EEG(end).setname), '',[' of ' EEG(end).setname])] ...
	                 '''electrodes'',''off''' };
    if length(EEG) > 1
        promptstr = { 'Dataset indices to subtract (Ex: ''1 2''-> 1-2)' promptstr{:} };
        inistr    = { '2 1' inistr{:} };
    end;
    result = inputdlg2( promptstr, 'Plot component and ERP envelopes -- pop_envtopo()', 1, inistr, 'pop_envtopo');
    if length(result) == 0 return; end;

    if length(EEG) > 1
        subindices = eval( [ '[' result{1} ']' ] );
        result(1) = [];
        EEG = EEG(subindices(1:2));
        fprintf('pop_envtopo(): Subtracting the epoch mean of dataset %d from that of dataset %d\n', ...
                   subindices(2), subindices(1));
    end;

    timerange    = eval( [ '[' result{1} ']' ] );
	if ~isempty( result{2} ), options = [ options '''limcontrib'',[' result{2} '],' ]; end;
	if ~isempty( result{4} ), options = [ options '''compnums'',[' result{4} '],' ]; 
    else                      options = [ options '''compnums'',-' result{3} ',' ];
    end;
    if ~isempty(result{5}),   options = [ options '''subcomps'',[' result{5} '],' ]; end;
    if ~isempty(result{6}),   options = [ options '''title'', ''' result{6} ''',' ]; end;
	options      =  [ options result{7} ];
	figure;
    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
else
    if isempty(timerange)
        timerange = [EEG.xmin*1000 EEG.xmax*1000];
    end
    options = [options vararg2str( varargin ) ];
end;
    
if length(EEG) > 2
    error('Cannot process more than two datasets');
end;

sigtmp = reshape(EEG(1).data, EEG(1).nbchan, EEG(1).pnts, EEG(1).trials);
if length(EEG) == 2
    if ~all(EEG(1).icaweights(:) == EEG(2).icaweights(:))
        error('The ICA decomposition must be the same for the two datasets');
    end;
    sigtmp2 = reshape(EEG(2).data, EEG(2).nbchan, EEG(2).pnts, EEG(2).trials);
end;
posi = round( (timerange(1)/1000-EEG(1).xmin) * EEG(1).srate) + 1;

posf = min(round( (timerange(2)/1000-EEG(1).xmin) * EEG(1).srate) + 1, EEG(1).pnts);

% outputs
% -------
outstr = '';
if nargin >= 4
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% generate output command
% ------------------------
if length( options ) < 2, options = ''; end;
if length(EEG) == 1
    varargout{1} = sprintf('figure; pop_envtopo(%s, [%s] %s);', ...
                                   inputname(1), num2str(timerange), options);
else
    if exist('subindices')
        varargout{1} = sprintf('figure; pop_envtopo(%s([%s]), [%s] %s);', ...
                                   inputname(1), int2str(subindices), num2str(timerange), options);
    end;
end;

% plot the data
% --------------
options = [ options ', ''verbose'', ''off''' ];
if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end;
if any(isnan(sigtmp(:)))
    disp('NaN detected: using nan_mean');
    if length(EEG) == 2
        com =  sprintf(['%s envtopo(nan_mean(sigtmp(:,posi:posf,:),3)-nan_mean(sigtmp2(:,posi:posf,:),3),' ...
                        'EEG(1).icaweights*EEG(1).icasphere, ' ...
                        '''chanlocs'', EEG(1).chanlocs, ''chaninfo'', EEG(1).chaninfo, ''icawinv'', EEG(1).icawinv,' ...
                        '''timerange'', [timerange(1) timerange(2)] %s);' ] , outstr, options);
    else % length(EEG) == 1
        com =  sprintf(['%s envtopo(nan_mean(sigtmp(:,posi:posf,:),3), EEG.icaweights*EEG.icasphere, ' ...
                        '''chanlocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo, ''icawinv'', EEG.icawinv,' ...
                        '''timerange'', [timerange(1) timerange(2)] %s);' ] , outstr, options);
    end;
else    
    if length(EEG) == 2
        com =  sprintf(['%s envtopo(mean(sigtmp(:,posi:posf,:),3)-mean(sigtmp2(:,posi:posf,:),3),' ...
                        ' EEG(1).icaweights*EEG(1).icasphere, ' ...
                        '''chanlocs'', EEG(1).chanlocs, ''chaninfo'', EEG(1).chaninfo, ''icawinv'', EEG(1).icawinv,' ...
                        '''timerange'', [timerange(1) timerange(2)] %s);' ] , outstr, options);
    else % length(EEG) == 1
        com =  sprintf(['%s envtopo(mean(sigtmp(:,posi:posf,:),3), EEG.icaweights*EEG.icasphere, ' ...
                        '''chanlocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo, ''icawinv'', EEG.icawinv,' ...
                        '''timerange'', [timerange(1) timerange(2)] %s);' ] , outstr, options);
    end;    
end;
com
eval(com);

return;

		
