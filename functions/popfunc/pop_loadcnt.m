% pop_loadcnt() - load a neuroscan CNT file (pop out window if no arguments).
%
% Usage:
%   >> [dat] = pop_loadcnt( filename, varargin);
%
% Inputs:
%   filename       - file name
%   varargin       - see LDCNT input
% 
% Outputs:
%   dat            - EEGLAB data structure
%
% Note: 
%   This function extract all non-null event from the CNT data structure.
% Null events are usually associated with internal signals (recalibrations...).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: loadcnt(), eeglab()

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
% Revision 1.4  2002/10/15 17:01:13  arno
% drawnow
%
% Revision 1.3  2002/08/12 02:40:59  arno
% inputdlg2
%
% Revision 1.2  2002/08/06 21:33:30  arno
% spelling
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, command] = pop_loadcnt(filename, varargin); 
command = '';

if nargin < 1 

	% ask user
	[filename, filepath] = uigetfile('*.CNT', 'Choose a CNT file -- pop_loadcnt()'); 
    drawnow;
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
	promptstr    = { 'Average reference ?' ...
					 'Enter block size in CNT file (1 or 40):' };
	inistr       = { 'YES', '1'  };
	pop_title    = sprintf('Load a CNT dataset');
	result       = inputdlg2( promptstr, pop_title, 1,  inistr, 'pop_loadcnt');
	if length( result ) == 0 return; end;

	% decode parameters
	% -----------------
    avgref = result{1};
    blockread = eval( result{2} );
    options = [ ', ''blockread'', ' int2str(blockread) ', ''avgref'', ''' avgref ''''];
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

% load datas
% ----------
EEG = eeg_emptyset;
if exist('filepath')
	fullFileName = sprintf('%s%s', filepath, filename);
else
	fullFileName = filename;
end;	
if nargin > 0
	if ~isempty(varargin)
		r = loadcnt( fullFileName, varargin{:});
	else
		r = loadcnt( fullFileName);
	end;	
else
	r = loadcnt( fullFileName, 'blockread', blockread, 'avgref', avgref);
end;

EEG.data            = r.dat;
EEG.filename        = filename;
EEG.filepath        = filepath;
EEG.setname 		= 'CNT file';
EEG.nbchan          = r.nchannels; 
I = find( r.event.stimtype > 0);

EEG.event(1:length(I),1) = 0;
EEG.event(1:length(I),2) = r.event.frame(I);
EEG.event = eeg_eventformat (EEG.event, 'struct', { 'type' 'latency' });

EEG.srate           = r.rate;
EEG = eeg_checkset(EEG);

command = sprintf('EEG = pop_loadcnt(''%s'' %s);',fullFileName, options); 

return;
