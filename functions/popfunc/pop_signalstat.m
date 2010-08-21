% pop_signalstat() - Computes and plots statistical characteristics of a signal,
%                    including the data histogram, a fitted normal distribution,
%                    a normal ditribution fitted on trimmed data, a boxplot, and
%                    the QQ-plot. The estimates value are printed in a panel and
%                    can be read as output. See SIGNALSTAT.
% Usage:
%   >>  OUTEEG = pop_signalstat( EEG, type ); % pops up
%   >>  [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = pop_signalstat( EEG, type, cnum );
%   >>  [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = pop_signalstat( EEG, type, cnum, percent );
%
% Inputs:
%   EEG   - input EEG dataset
%   type  - type of processing
%           1: process the raw  data; 0: the ICA components
%   cnum  - selected channel or component
%    
% Outputs:
%   OUTEEG  - output dataset
%
% Author: Luca Finelli, CNL / Salk Institute - SCCN, 2 August 2002
%
% See also:
%   SIGNALSTAT,  EEGLAB 

% Copyright (C) 2002 Luca Finelli, Salk/SCCN, La Jolla, CA
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

function varargout = pop_signalstat( EEG, typeproc, cnum, percent );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
varargout{1} = '';

% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_signalstat;
	return;
end;	
popup=0;
if nargin < 3
	popup = 1;
end;
if nargin < 4
	percent=5;
end;

% pop up window
% -------------
if (nargin < 3 & typeproc==1)
	promptstr    = { 'Channel number:'; 'Trim percentage (each end):' };
	inistr       = { '1';'5' };
	result       = inputdlg2( promptstr, 'Plot signal statistics -- pop_signalstat()', 1,  inistr, 'signalstat');
	if length( result ) == 0 return; end;
	cnum   	     = eval( [ '[' result{1} ']' ] ); % the brackets allow processing Matlab arrays
	percent      = eval( [ '[' result{2} ']' ] );
elseif (nargin < 3 & typeproc==0)
	promptstr    = { 'Component number:'; 'Trim percentage (each end):' };
	inistr       = { '1'; '5' };
	result       = inputdlg2( promptstr, 'Plot signal statistics -- pop_signalstat()', 1,  inistr, 'signalstat');
	if length( result ) == 0 return; end;
	cnum    	 = eval( [ '[' result{1} ']' ] ); % the brackets allow processing Matlab arrays
    percent      = eval( [ '[' result{2} ']' ] );
end;

if length(cnum) ~= 1 | (cnum-floor(cnum)) ~= 0
	error('pop_signalstat(): Channel/component number must be a single integer');
end

if cnum < 1 | cnum > EEG.nbchan
   error('pop_signalstat(): Channel/component number out of range');
end;   

% call function signalstat() either on raw data or ICA data
% ---------------------------------------------------------
if typeproc == 1
	tmpsig=EEG.data(cnum,:);
%	[M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh]=signalstat( EEG.data(cnum,:),1,[], percent);
	dlabel=[];
	dlabel2=['Channel ' num2str(cnum)];
	map = cnum;
else 
	if ~isempty( EEG.icasphere )
        tmpsig = eeg_getdatact(EEG, 'component', cnum);
        tmpsig = tmpsig(:,:);
	%	[M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh]=signalstat( tmpsig,1,'Component Activity',percent);
		dlabel='Component Activity';
		dlabel2=['Component ' num2str(cnum)];
		map = EEG.icawinv(:,cnum);
	else
		error('You must run ICA first');
	end;	
end;	 

% outputs
% -------
outstr = '';
if ~popup
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% return the string command
% -------------------------
%fprintf('Pop_signalstat: computing statistics...\n');
varargout{1} = sprintf('pop_signalstat( %s, %d, %d );', inputname(1), typeproc, cnum);


plotloc = 0;
if ~isempty(EEG.chanlocs)
    if isfield(EEG.chanlocs, 'theta')
        if ~isempty(EEG.chanlocs(cnum).theta)
            plotloc = 1;
        end;
    end;
end;
if plotloc
    if typeproc == 1
        com = sprintf('%s signalstat( tmpsig, 1, dlabel, percent, dlabel2, map, EEG.chanlocs );', outstr);
    else
        com = sprintf('%s signalstat( tmpsig, 1, dlabel, percent, dlabel2, map, EEG.chanlocs(EEG.icachansind) );', outstr);
    end;
else
    com = sprintf('%s signalstat( tmpsig, 1, dlabel, percent, dlabel2);', outstr);
end;

eval(com)	
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;

return;
