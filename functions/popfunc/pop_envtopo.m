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
%   'key','val' - optional envtopo() and topoplot() arguments 
%                 (see >> help topoplot())
%
% Outputs: Same as envtopo(). When nargin < 3, a query window pops-up 
%          to ask for additional arguments and no outputs are returned.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: envtopo(), eeglab()

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
			 'Number of largest contributing components to plot (7):', ...
			 'Else plot these component numbers only (Ex: 2:4,7):', ...
                         'Component numbers to remove from data before plotting:' ...
			 'Plot title:' ...
			 'Optional topoplot() and envtopo() arguments:' };
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
    if ~isempty( result{3} ), options = [ options '''compsplot'',[' result{3} '],' ]; end;
    if ~isempty( result{4} ), options = [ options '''compnums'',[' result{4} '],' ]; end;
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

if timerange(1) < max([EEG.xmin])*1000, timerange(1) =  max([EEG.xmin])*1000; end;
if timerange(2) > min([EEG.xmax])*1000, timerange(2) =  min([EEG.xmax])*1000; end;

EEG1 = eeg_checkset(EEG(1),'loaddata');
sigtmp = reshape(EEG1.data, EEG1.nbchan, EEG1.pnts, EEG1.trials);
if ~isempty(EEG1.icachansind), sigtmp = sigtmp(EEG1.icachansind,:,:); end;
if length(EEG) == 2
    EEG2 = eeg_checkset(EEG(2),'loaddata');
    if ~all(EEG1.icaweights(:) == EEG2.icaweights(:))
        error('The ICA decomposition must be the same for the two datasets');
    end;
    sigtmp2 = reshape(EEG2.data, EEG2.nbchan, EEG2.pnts, EEG2.trials);
    if ~isempty(EEG2.icachansind), sigtmp2 = sigtmp2(EEG2.icachansind,:,:); end;
end;
posi = round( (timerange(1)/1000-EEG1.xmin) * EEG1.srate) + 1;
posf = min(round( (timerange(2)/1000-EEG1.xmin) * EEG1.srate) + 1, EEG1.pnts);

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

% fprintf(['\npop_envtopo(): Issuing command: ' com '\n\n']); % type the evntopo() call

eval(com); % make the plot using envtopo()

return;

		
