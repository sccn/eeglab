% pop_eventstat() - Computes and plots statistical characteristics of an EEG event,
%                   including the data histogram, a fitted normal distribution,
%                   a normal ditribution fitted on trimmed data, a boxplot, and
%                   the QQ-plot. The estimates value are printed in a panel and
%                   can be read as output. NaNs are omitted. See signalstat().
%
% Usage:
%   >>  OUTEEG = pop_eventstat( EEG ); % pops up
%   >>  [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = pop_eventstat( EEG, eventfield, type );
%   >>  [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = pop_eventstat( EEG, eventfield, type, percent );
%
% Inputs:
%   EEG        - input EEG dataset
%   eventfield - event field to process (i.e. latency)
%   type       - name of the event type(s) to process. Can be a single element or
%                a cell array. Default is all types.
%   latrange   - [min max] event latency range within data epochs in milliseconds.
%                Default is whole epoch.
%   percent    - percentage for trimmed data statistics. Default is 5%. (see signalstat())
%    
% Outputs:
%   OUTEEG  - output dataset
%
% Author: Arnaud Delorme & Luca Finelli, CNL / Salk Institute - SCCN, 15 August 2002
%
% See also: signalstat(), eeg_getepochevent(), eeglab()

% Copyright (C) 2002 Arnaud Delorme & Luca Finelli, Salk/SCCN, La Jolla, CA
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

function varargout = pop_eventstat( EEG, eventfield, type, latrange, percent );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
varargout{1} = '';

% display help if not enough arguments
% ------------------------------------
if nargin < 1
	help pop_eventstat;
	return;
end;	
popup=0;
if nargin < 2
	popup = 1;
end
if nargin < 3
	percent=5;
end

% pop up window
% -------------
if nargin < 2
	promptstr    = { 'Event field to process:' ...
					 strvcat('Event type(s) ([]=all):', ...
							 'Select "Edit > Event values" to see type values') ...
					strvcat('Event latency range (ms)', ...
							'Default is whole epoch or data') ...
					'Percent for trimmed statistics:' };
	inistr       = { 'latency' '' '' '5' };
	result       = inputdlg2( promptstr, 'Plot event statistics -- pop_eventstat()', 1,  inistr, 'signalstat');
	if length( result ) == 0 return; end
	eventfield   = deblank(result{1}); % the brackets allow to process matlab arrays
    if ~isempty(result{2})
        if strcmpi(result{2}(1),'''')
             type = eval( [ '{' result{2} '}' ] );
        else type = parsetxt( result{2});
        end
    else
        disp('WARNING: you should select an event type');
        type = {};
    end
	latrange     = eval( [ '[' result{3} ']' ] );
	percent      = eval( [ '[' result{4} ']' ] );
else
    if nargin < 3
        type = [];
    end
    if nargin < 4
        latrange = [];
    end
    if nargin < 5
        percent = 5;
    end
end

% call function signalstat() either on raw data or ICA data
% ---------------------------------------------------------
[ typevals alltypevals ] = eeg_getepochevent(EEG, type, latrange, eventfield);
% concatenate alltypevals
% -----------------------
typevals = [];
for index = 1:length(alltypevals)
    typevals = [ typevals alltypevals{index} ];
end;   
if isempty(typevals)
    error('No such events found. See Edit > Event values to confirm event type.');
end
dlabel='Event values';
if isempty(type)
    dlabel2=['All event statistics for ''' eventfield ''' info'];
else
    dlabel2=['Event ' vararg2str(type) ' statistics for ''' eventfield ''' info'];
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
fprintf('pop_eventstat: extracting events...\n');
varargout{1} = sprintf('pop_eventstat( EEG, %s );', vararg2str({eventfield type latrange percent}));
com          = sprintf('%s signalstat( typevals, 1, dlabel, percent, dlabel2 ); %s', outstr);

eval(com)	
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end

return;
