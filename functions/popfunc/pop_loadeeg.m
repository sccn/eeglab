% pop_loadeeg() - load a neuroscan EEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_loadeeg; % pop-up window mode
%   >> EEG = pop_loadeeg( filename, filepath, range_chan, range_trials, ...
%                  range_typeeeg, range_response);
%
% Graphical interface:
%   "Enter trial range subset" - [edit box] integer array. Command line
%                equivalent: 'range_trials'
%   "Enter type range subset" - [edit box] integer array. Command line
%                equivalent: 'range_typeeeg'
%   "Enter electrode subset" - [edit box] integer array. Command line
%                equivalent: 'range_chan'
%   "Enter response range subset" - [edit box] integer array. Command line
%                equivalent: 'range_response'
%
% Inputs:
%   filename       - file name
%   filepath       - file path
%   range_chan     - [integer array] only import selected electrodes
%                    (ex: 3,4:10; default all)
%   range_trials   - [integer array] only import selected trials (default all)
%   range_typeeeg  - [integer array] only import trials with selected type
%                    (default all)
%   range_response - [integer array] only import trials with selected responses
%                    values (default all)
% 
% Outputs:
%   EEG            - eeglab() data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: loadeeg(), eeglab()

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
% Revision 1.4  2003/01/14 00:30:23  arno
% handling the case where all rts are 0
%
% Revision 1.3  2002/08/12 02:38:13  arno
% [6~[6~inputdlg2
%
% Revision 1.2  2002/05/02 19:39:43  arno
% updating function for new event structure
% ,
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

% uses calls to eeg_emptyset and loadeeg

% popup loadeeg file
% ------------------
function [EEG, command] = pop_loadeeg(filename, filepath, range_chan, range_sweeps, range_typeeeg, range_response); 

command = '';

if nargin < 2

	% ask user
	[filename, filepath] = uigetfile('*.eeg', 'Choose an EEG file -- pop_loadeeg()'); 
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
	promptstr    = { 'Enter trialrange subset:', ...
					 'Enter the type range subset:', ...
					 'Enter electrodes subset:', ...
					 'Enter response range subset:'};
	inistr       = { '' '' '' '' };
	pop_title    = sprintf('Load an EEG dataset');
	result       = inputdlg2( promptstr, pop_title, 1,  inistr, 'pop_loadeeg');
	if size( result,1 ) == 0 return; end;

	% decode parameters
	% -----------------
	range_sweeps    = eval( [ '[' result{1} ']' ] );
	range_typeeeg   = eval( [ '[' result{2}  ']' ] );
	range_chan      = eval( [ '[' result{3}  ']' ] );
	range_response  = eval( [ '[' result{4}  ']' ] );
end;

if ~exist('range_chan')   | isempty(range_chan)      , range_chan     = 'all'; end;
if ~exist('range_sweeps') | isempty(range_sweeps)    , range_sweeps     = 'all'; end;
if ~exist('range_typeeg') | isempty(range_typeeeg)   , range_typeeeg     = 'all'; end;
if ~exist('range_response') | isempty(range_response), range_response     = 'all'; end;

% load datas
% ----------
EEG = eeg_emptyset;
fullFileName = sprintf('%s%s', filepath, filename);
[EEG.data, accept, eegtype, rt, eegresp, EEG.namechan, EEG.pnts, EEG.trials, EEG.srate, EEG.xmin, EEG.xmax] = loadeeg( fullFileName, range_chan, range_sweeps, range_typeeeg, 'all', 'all', range_response);
EEG.filename        = filename;
EEG.filepath        = filepath;
EEG.setname 		= 'rawdatas';
EEG.nbchan          = size(EEG.data,1);
EEG = eeg_checkset(EEG);
if any(rt)
    EEG = pop_importepoch( EEG, [rt(:)*1000 eegtype(:) accept(:) eegresp(:)], { 'RT' 'type' 'accept' 'response'}, {'RT'}, 1E-3, 0, 1);
else
    EEG = pop_importepoch( EEG, [eegtype(:) accept(:) eegresp(:)], { 'type' 'accept' 'response'}, { }, 1E-3, 0, 1);
end;    
command = sprintf('EEG = pop_loadeeg(''%s'', ''%s'', %s);', filename, filepath, ...
			vararg2str({range_chan range_sweeps range_typeeeg range_response }));

return;



