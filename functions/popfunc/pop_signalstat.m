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
end
if nargin < 4
	percent=5;
end

% pop up window
% -------------
if (nargin < 3 && typeproc==1)
	promptstr    = { 'Channel number:'; 'Trim percentage (each end):' };
	inistr       = { '1';'5' };
	result       = inputdlg2( promptstr, 'Plot signal statistics -- pop_signalstat()', 1,  inistr, 'signalstat');
	if length( result ) == 0 return; end
	cnum   	     = eval( [ '[' result{1} ']' ] ); % the brackets allow processing Matlab arrays
	percent      = eval( [ '[' result{2} ']' ] );
elseif (nargin < 3 && typeproc==0)
	promptstr    = { 'Component number:'; 'Trim percentage (each end):' };
	inistr       = { '1'; '5' };
	result       = inputdlg2( promptstr, 'Plot signal statistics -- pop_signalstat()', 1,  inistr, 'signalstat');
	if length( result ) == 0 return; end
	cnum    	 = eval( [ '[' result{1} ']' ] ); % the brackets allow processing Matlab arrays
    percent      = eval( [ '[' result{2} ']' ] );
end

if length(cnum) ~= 1 || (cnum-floor(cnum)) ~= 0
	error('pop_signalstat(): Channel/component number must be a single integer');
end

if cnum < 1 || cnum > EEG.nbchan
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
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end
end

% return the string command
% -------------------------
%fprintf('Pop_signalstat: computing statistics...\n');
varargout{1} = sprintf('pop_signalstat( EEG, %d, %d );', typeproc, cnum);


plotloc = 0;
if ~isempty(EEG.chanlocs)
    if isfield(EEG.chanlocs, 'theta')
        if ~isempty(EEG.chanlocs(cnum).theta)
            plotloc = 1;
        end
    end
end
if plotloc
    if typeproc == 1
        com = sprintf('%s signalstat( tmpsig, 1, dlabel, percent, dlabel2, map, EEG.chanlocs );', outstr);
    else
        com = sprintf('%s signalstat( tmpsig, 1, dlabel, percent, dlabel2, map, EEG.chanlocs(EEG.icachansind) );', outstr);
    end
else
    com = sprintf('%s signalstat( tmpsig, 1, dlabel, percent, dlabel2);', outstr);
end

eval(com)	
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end

return;
