% pop_loadcnt() - load a neuroscan CNT file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_loadcnt; % pop-up window mode
%   >> EEG = pop_loadcnt( filename, 'key', 'val', ...);
%
% Graphic interface:
%   "Enter block size in CNT file" - [edit box] Neuroscan files can
%                   have different data block sizes. It was not
%                   possible to read this information in the file header
%                   (though it is probably there). Values known to 
%                   work are 1 and 40. If 1 does not work (the data does
%                   not look like EEG), you should try 40. 
%   "loadcnt() 'key', 'val' params" - [edit box] Enter optional loadcnt()
%                   parameters.
%
% Inputs:
%   filename       - file name
%
% Optional inputs:
%   Same as loadcnt() function.
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Note: 
% 1) This function extract all non-null event from the CNT data structure.
% Null events are usually associated with internal signals (recalibrations...).
% 2) The "Average reference" edit box had been remove since the re-referencing
% menu of EEGLAB offers more options to re-reference data.
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
% Revision 1.9  2003/04/10 17:56:52  arno
% debuging function and history
%
% Revision 1.8  2003/04/10 17:35:29  arno
% header and history
%
% Revision 1.7  2003/03/05 19:48:42  arno
% removing matlab warning
%
% Revision 1.6  2002/11/23 21:10:02  arno
% importing type of event
%
% Revision 1.5  2002/10/22 23:57:21  arno
% change default blockread
%
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
EEG = [];

if nargin < 1 

	% ask user
	[filename, filepath] = uigetfile('*.CNT;*.cnt', 'Choose a CNT file -- pop_loadcnt()'); 
    drawnow;
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
   promptstr    = { 'Enter block size for CNT file (1 or 40):' ...
                    'loadcnt() ''key'', ''val'' params' };
	inistr       = { '1'  '' };
	pop_title    = sprintf('Load a CNT dataset');
	result       = inputdlg2( promptstr, pop_title, 1,  inistr, 'pop_loadcnt');
	if length( result ) == 0 return; end;

	% decode parameters
	% -----------------
    blockread = eval( result{1} );
    options = [ ', ''blockread'', ' int2str(blockread) ];
    if ~isempty(result{2}), options = [ options ',' result{2} ]; end;
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
	eval( [ 'r = loadcnt( fullFileName ' options ');' ]);
end;

EEG.data            = r.dat;
EEG.filename        = filename;
EEG.setname 		= 'CNT file';
EEG.nbchan          = r.nchannels; 
I = find( ( r.event.stimtype ~= 0 ) & ( r.event.stimtype ~= 255 ) );

EEG.event(1:length(I),1) = r.event.stimtype(I);
EEG.event(1:length(I),2) = r.event.frame(I);
EEG.event = eeg_eventformat (EEG.event, 'struct', { 'type' 'latency' });

EEG.srate           = r.rate;
EEG = eeg_checkset(EEG);

if length(options) > 2
    command = sprintf('EEG = pop_loadcnt(''%s'' %s);',fullFileName, options); 
else
    command = sprintf('EEG = pop_loadcnt(''%s'');',fullFileName); 
end;
return;
