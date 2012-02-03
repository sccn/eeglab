% pop_loadeeg() - load a Neuroscan .EEG file (via a pop-up window if no
%                  arguments). Calls loadeeg().
%
% Usage:
%   >> EEG = pop_loadeeg; % pop-up data entry window
%   >> EEG = pop_loadeeg( filename, filepath, range_chan, range_trials, ...
%                  range_typeeeg, range_response, format); % no pop-up window
%
% Graphic interface:
%   "Data precision in bits..." - [edit box] data binary format length
%                in bits. Command line equivalent: 'format'
%   "Trial range subset" - [edit box] integer array. 
%                Command line equivalent: 'range_trials'
%   "Type range subset" - [edit box] integer array. 
%                Command line equivalent: 'range_typeeeg'
%   "Electrode subset" - [edit box] integer array. 
%                Command line equivalent: 'range_chan'
%   "Response range subset" - [edit box] integer array. 
%                Command line equivalent: 'range_response'
%
% Inputs:
%   filename       - ['string'] file name
%   filepath       - ['string'] file path
%   range_chan     - [integer array] Import only selected electrodes
%                    Ex: 3,4:10; {Default: [] -> import all}
%   range_trials   - [integer array] Import only selected trials
%                    { Default: [] -> import all}
%   range_typeeeg  - [integer array] Import only trials of selected type
%                    {Default: [] -> import all}
%   range_response - [integer array] Import only trials with selected 
%                    response values {Default: [] -> import all}
%   format         - ['short'|'int32'] data binary format (Neuroscan 4.3
%                    saves data as 'int32'; earlier versions save data as
%                    'short'. Default is 'short'.
% Outputs:
%   EEG            - eeglab() data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: loadeeg(), eeglab()

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

% 01-25-02 reformated help & license -ad 

% uses calls to eeg_emptyset and loadeeg

% popup loadeeg file
% ------------------
function [EEG, command] = pop_loadeeg(filename, filepath, range_chan, range_sweeps, range_typeeeg, range_response, datformat); 

EEG = [];
command = '';

if nargin < 1

	% ask user
	[filename, filepath] = uigetfile('*.eeg;*.EEG', 'Choose an EEG file -- pop_loadeeg()'); 
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
	promptstr    = { 'Data precision in bits (16 / 32 bit or Auto for NS v4.3):', ...
                     'Trial range subset:', ...
					 'Type range subset:', ...
					 'Electrodes subset:', ...
					 'Response range subset:'};
	inistr       = { 'Auto' '' '' '' '' };
	pop_title    = sprintf('Load an EEG dataset');
	result       = inputdlg2( promptstr, pop_title, 1,  inistr, 'pop_loadeeg');
	if size( result,1 ) == 0 return; end;

	% decode parameters
	% -----------------
    precision = lower(strtrim(result{1}));
    if strcmpi(precision, '16')
        datformat = 'short';
    elseif strcmpi(precision, '32')
        datformat = 'int32';
    elseif (strcmpi(precision, '0') || strcmpi(precision, 'auto'))
        datformat = 'auto'
    end;
	range_sweeps    = eval( [ '[' result{2} ']' ] );
	range_typeeeg   = eval( [ '[' result{3}  ']' ] );
	range_chan      = eval( [ '[' result{4}  ']' ] );
	range_response  = eval( [ '[' result{5}  ']' ] );
else
    if exist('filepath') ~= 1
        filepath = '';
    end;
end;

if exist('datformat') ~= 1, datformat = 'auto'; end;
if exist('range_chan') ~= 1   | isempty(range_chan)      , range_chan     = 'all'; end;
if exist('range_sweeps') ~= 1 | isempty(range_sweeps)    , range_sweeps     = 'all'; end;
if exist('range_typeeeg') ~= 1 | isempty(range_typeeeg)   , range_typeeeg     = 'all'; end;
if exist('range_response') ~= 1 | isempty(range_response), range_response     = 'all'; end;

% load datas
% ----------
EEG = eeg_emptyset;
if ~isempty(filepath)
    if filepath(end) ~= '/' & filepath(end) ~= '\' & filepath(end) ~= ':'
        error('The file path last character must be a delimiter');
    end;
    fullFileName = sprintf('%s%s', filepath, filename);
else
    fullFileName = filename;
end;    
[EEG.data, accept, eegtype, rt, eegresp, namechan, EEG.pnts, EEG.trials, EEG.srate, EEG.xmin, EEG.xmax] = ...
    loadeeg( fullFileName, range_chan, range_sweeps, range_typeeeg, 'all', 'all', range_response, datformat);

EEG.comments        = [ 'Original file: ' fullFileName ];
EEG.setname 		= 'Neuroscan EEG data';
EEG.nbchan          = size(EEG.data,1);
for index = 1:size(namechan,1)
    EEG.chanlocs(index).labels = deblank(char(namechan(index,:)));
end;
EEG = eeg_checkset(EEG);
if any(rt)
    EEG = pop_importepoch( EEG, [rt(:)*1000 eegtype(:) accept(:) eegresp(:)], { 'RT' 'type' 'accept' 'response'}, {'RT'}, 1E-3, 0, 1);
else
    EEG = pop_importepoch( EEG, [eegtype(:) accept(:) eegresp(:)], { 'type' 'accept' 'response'}, { }, 1E-3, 0, 1);
end;    
command = sprintf('EEG = pop_loadeeg(''%s'', ''%s'', %s);', filename, filepath, ...
			vararg2str({range_chan range_sweeps range_typeeeg range_response datformat }));

return;



