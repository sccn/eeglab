% pop_mergeset() - merge two datasets. If only one argument is given
%               a window pop out to ask for other arguments.
%
% Usage:
%   >> OUTEEG = pop_mergeset( INEEG1, INEEG2, keepall);
%
% Inputs:
%  INEEG1  - first input dataset or dataset number
%  INEEG2  - second  input dataset or dataset number
%  keepall - [0|1] 0 remove ICA and 1 preserve ICA of the first dataset
%            and recompute the activity on the merged data.
%
% Outputs:
%  OUTEEG  - merged dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

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
% Revision 1.2  2002/04/10 21:32:26  arno
% debuging event concatenation
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 01-26-02 change format for events and trial conditions -ad

function [EEG, com] = pop_mergeset( EEG, INEEG2, keepall);

com = '';
if nargin < 1
	help pop_mergeset;
	return;
end;
if isempty(EEG.data)
   disp('Pop_merge error: cannot merge empty dataset'); return;
end;    
originput1 = EEG;
if ~isstruct( EEG )
    EEG = eeg_retrieve( EEG );	
end;	
if nargin < 2
   promptstr    = { 'Enter number of dataset to merge with', ...
      				'Preserve ICA of the first dataset ?' };
	inistr       = { '1', 'no' };
	result       = inputdlg( promptstr, 'Merge datasets -- pop_mergeset()', 1,  inistr);
	size_result  = size( result );
	if size_result(1) == 0 return; end;
   
    INEEG2  = eval( result{1} );
    switch lower(result{2})
    	case 'keepall', keepall = 1;
    	otherwise, keepall = 0;
    end;
else
	if nargin < 3
		keepall = 0;
	end;	
end;

originput2 = INEEG2;
if ~isstruct( INEEG2 )
	INEEG2 = eeg_retrieve( INEEG2 );	
end;	

% check consistency
% -----------------
if EEG.nbchan ~= INEEG2.nbchan
	error('The two datasets must have the same number of channels');
end;	
if EEG.srate ~= INEEG2.srate
	error('The two datasets must have the same sampling rate');
end;	
if EEG.trials > 1 | INEEG2.trials > 1
	if EEG.pnts ~= INEEG2.pnts
		error('The two epoched datasets must have the same number of points');
	end;
	if EEG.xmin ~= INEEG2.xmin
		fprintf('Warning: the two epoched datasets do not have the same time onset, adjusted');
	end;
	if EEG.xmax ~= INEEG2.xmax
		fprintf('Warning: the two epoched datasets do not have the same time offset, adjusted');
	end;
end;	

fprintf('Merging the two datasets...\n');		
EEG.data    = [ EEG.data(:,:) INEEG2.data(:,:) ];
EEG.setname	= 'Merge datasets';

% concatenate events
% ------------------
if ~isempty(INEEG2.event)
	if isfield( EEG.event, 'epoch')
		for index = 1:length(INEEG2.event(:))
			INEEG2.event(index).epoch = INEEG2.event(index).epoch + EEG.trials;
		end;    
	end;
	if isfield( EEG.event, 'latency')
		for index = 1:length(INEEG2.event(:))
			INEEG2.event(index).latency = INEEG2.event(index).latency + EEG.trials*EEG.pnts;
		end;    
	end;
			
	EEG.event(end+1:end+length(INEEG2.event)) = INEEG2.event(:);			
end;

% concatenate epoch
% ---------------------
EEG.epoch = [];
% $$$ if ~isempty(INEEG2.epoch) & ~isempty(EEG.epoch)
% $$$     fields1 = fieldnames( EEG.epoch );
% $$$     fields2 = fieldnames( INEEG2.epoch );
% $$$     fieldsdiff = setdiff( fields2, fields1 );
% $$$     if ~isempty(fieldsdiff)
% $$$         for f = 1:length(fieldsdiff)
% $$$             for index = 1:EEG.trials
% $$$                 eval( [ 'EEG.epoch(' int2str(index) ').' fieldsdiff{f} ' = 0;' ] );
% $$$              end;
% $$$         end;
% $$$     end;
% $$$     fieldsdiff = setdiff( fields1, fields2 );
% $$$     fieldscom  = intersect( fields2, fields1 );
% $$$     for f = 1:length(fieldsdiff)
% $$$         for index = 1:INEEG2.trials
% $$$             eval( [ 'EEG.epoch(EEG.trials+' int2str(index) ').' fieldsdiff{f} ' = 0;' ] );
% $$$         end;
% $$$     end;
% $$$     for f = 1:length(fieldscom)
% $$$         for index = 1:INEEG2.trials
% $$$             eval( [ 'EEG.epoch(EEG.trials+' int2str(index) ').' fieldscom{f} ...
% $$$                     '= INEEG2.epoch(' int2str(index) ').' fieldscom{f} ';' ] );
% $$$         end;
% $$$     end;
% $$$ end;             

if EEG.trials > 1 | INEEG2.trials > 1
	EEG.trials  =  EEG.trials + INEEG2.trials;
else
	EEG.pnts = EEG.pnts + INEEG2.pnts;
end;

EEG = rmfield(EEG, 'reject' );
EEG.specicaact = [];
EEG.specdata = [];
if keepall == 0
    EEG.icaact = [];
    EEG.icawinv = [];
    EEG.icasphere = [];
    EEG.icaweights = [];
	EEG = rmfield(EEG, 'stats' );
else
	EEG.icaact = [];
end;
EEG = eeg_checkset(EEG);

% build the command
% -----------------
if isstruct( originput1 )	input1 = inputname(1);
else						input1 = num2str( originput1 );
end;
if nargin < 2
	com = sprintf('EEG = pop_mergeset( %s, %d, %d);', input1, originput2, keepall);		
else
	if isstruct( originput2 )	input2 = inputname(2);
	else						input2 = num2str( originput2 );
	end;
	com = sprintf('EEG = pop_mergeset( %s, %s, %d);', input1, input2, keepall);		
end;

return;
