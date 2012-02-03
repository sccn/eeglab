% pop_loadcnt() - load a neuroscan CNT file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_loadcnt; % pop-up window mode
%   >> EEG = pop_loadcnt( filename, 'key', 'val', ...);
%
% Graphic interface:
%   "Data fomat" - [checkbox] 16-bits or 32-bits. We couldn't find in the
%                   data file where this information was stored. Command
%                   line equivalent in loadcnt() 'dataformat'.
%   "Time interval in seconds" - [edit box] specify time interval [min max]
%                   to import portion of data. Command line equivalent
%                   in loadcnt: 't1' and 'lddur'
%   "Import keystrokes" - [checkbox] set this option to import keystroke
%                   event types in dataset. Command line equivalent
%                   'keystroke'.
%   "loadcnt() 'key', 'val' params" - [edit box] Enter optional loadcnt()
%                   parameters.
%
% Inputs:
%   filename       - file name
%
% Optional inputs:
%   'keystroke'    - ['on'|'off'] set the option to 'on' to import 
%                    keystroke event types. Default is off.
%   'memmapfile'   - ['memmapfile_name'] use this option if the .cnt file
%                    is too large to read in conventially.  The suffix of 
%                    the memmapfile_name must be .fdt.  the memmapfile
%                    functions process files based on their suffix and an
%                    error will occur if you use a different suffix.
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
% 3) The 'blockread' has been disabled since we found where this information
% was stored in the file.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: loadcnt(), eeglab()

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
    callback16 = 'set(findobj(gcbf, ''tag'', ''32''), ''value'', ~get(gcbo, ''value'')); set(findobj(gcbf, ''tag'', ''AD''), ''value'', ~get(gcbo, ''value''));';
    callback32 = 'set(findobj(gcbf, ''tag'', ''16''), ''value'', ~get(gcbo, ''value'')); set(findobj(gcbf, ''tag'', ''AD''), ''value'', ~get(gcbo, ''value''));';
    callbackAD = 'set(findobj(gcbf, ''tag'', ''16''), ''value'', ~get(gcbo, ''value'')); set(findobj(gcbf, ''tag'', ''32''), ''value'', ~get(gcbo, ''value''));';
    uigeom       = { [1.3 0.5 0.5 0.5] [1 0.5] [1.09 0.13 0.4] [1 0.5] [1 0.5] 1 } ;
    uilist       = { { 'style' 'text' 'string' 'Data format 16 or 32 bit (Default = Autodetect)' } ...
                     { 'style' 'checkbox' 'tag' '16' 'string' '16-bits' 'value' 0 'callback' callback16 } ...
                     { 'style' 'checkbox' 'tag' '32' 'string' '32-bits' 'value' 0 'callback' callback32 } ...
                     { 'style' 'checkbox' 'tag' 'AD' 'string' 'Autodetect' 'value' 1 'callback' callbackAD } ...
                     { 'style' 'text' 'string' 'Time interval in s (i.e. [0 100]):' } ...
                     { 'style' 'edit' 'string' '' 'callback' 'warndlg2([ ''Events latency might be innacurate when'' 10 ''importing time intervals (this is an open issue)'']);' } ...                  
                     { 'style' 'text' 'string' 'Check to Import keystrokes:' } ...
                     { 'style' 'checkbox' 'string' '' } { } ...                     
                     { 'style' 'text' 'string' 'loadcnt() ''key'', ''val'' params' } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' [ 'Large files, enter a file name for memory mapping (xxx.fdt)' ]  } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' '    Note: requires to enable memory mapping in EEGLAB memory options and only works for 32-bit files' } };                  
	result = inputgui( uigeom, uilist, 'pophelp(''pop_loadcnt'')', 'Load a CNT dataset');    
	if length( result ) == 0 return; end;

	% decode parameters
	% -----------------
    options = [];
    if result{1}, options = [ options ', ''dataformat'', ''int16''' ];
    elseif result{2}, options = [ options ', ''dataformat'', ''int32''' ];
    elseif result{3}, options = [ options ', ''dataformat'', ''auto''' ];
    end;
    if ~isempty(result{4}), 
        timer =  eval( [ '[' result{4} ']' ]);
        options = [ options ', ''t1'', ' num2str(timer(1)) ', ''lddur'', '  num2str(timer(2)-timer(1)) ]; 
    end;   
    if result{5}, options = [ options ', ''keystroke'', ''on''' ]; end;
    if ~isempty(result{6}), options = [ options ',' result{6} ]; end;
    % Conditional pass if ~isempty(result{7}), options = ... 
    % [options ', ''memmapfile''', result{7} ] ; end ;
    % Always pass the memmapfile paramter? 
    options = [ options ', ''memmapfile'', ', 'result{7}' ] ;
else
	options = vararg2str(varargin);
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

if isfield(r, 'dat')
    error('pop_loadcnt is not compatible with current loadcnt version, please use latest loadcnt() version');
end;
% Check to see if data is in memory or in a file.
EEG.data            = r.data;
EEG.comments        = [ 'Original file: ' fullFileName ];
EEG.setname 		= 'CNT file';
EEG.nbchan          = r.header.nchannels;

% inport events
% -------------
I = 1:length(r.event);
if ~isempty(I)
    EEG.event(1:length(I),1) = [ r.event(I).stimtype ];
    EEG.event(1:length(I),2) = [ r.event(I).offset ]+1;
    EEG.event = eeg_eventformat (EEG.event, 'struct', { 'type' 'latency' });
end;

% modified by Andreas Widmann  2005/05/12  14:15:00
try, % this piece of code makes the function crash sometimes - Arnaud Delorme 2006/04/27
    temp = find([r.event.accept_ev1] == 14 | [r.event.accept_ev1] == 11); % 14: Discontinuity, 11: DC reset
    if ~isempty(temp)
        disp('pop_loadcnt note: event field ''type'' set to ''boundary'' for data discontinuities');
        for index = 1:length(temp)
            EEG.event(temp(index)).type = 'boundary';
        end;
    end
catch, end;
% end modification

% process keyboard entries
% ------------------------
if ~isempty(findstr('keystroke', lower(options)))
    tmpkbd  = [ r.event(I).keyboard ];
    tmpkbd2 = [ r.event(I).keypad_accept ];
    for index = 1:length(EEG.event)
        if EEG.event(index).type == 0
            if r.event(index).keypad_accept,
                EEG.event(index).type = [ 'keypad' num2str(r.event(index).keypad_accept) ];
            else
                EEG.event(index).type = [ 'keyboard' num2str(r.event(index).keyboard) ];
            end;
        end;
    end;
else
    % removeing keystroke events
    % --------------------------
    rmind = [];
    for index = 1:length(EEG.event)
        if EEG.event(index).type == 0
            rmind = [rmind index];
        end;
    end;
    if ~isempty(rmind)
        fprintf('Ignoring %d keystroke events\n', length(rmind));
        EEG.event(rmind) = [];
    end;
end;

% import channel locations (Neuroscan coordinates are not wrong)
% ------------------------
%x            = celltomat( { r.electloc.x_coord } );
%y            = celltomat( { r.electloc.y_coord } );
for index = 1:length(r.electloc)
    names{index} = deblank(char(r.electloc(index).lab'));
    if size(names{index},1) > size(names{index},2), names{index} = names{index}'; end;
end;
EEG.chanlocs  = struct('labels', names);
%EEG.chanlocs = readneurolocs( { names x y } );
%disp('WARNING: Electrode locations imported from CNT files may not reflect true locations');

% Check to see if data is in a file or in memory
% If in memory, leave alone
% If in a file, use values set in loadcnt.m for nbchan and pnts.
EEG.srate    = r.header.rate;
EEG.nbchan   = size(EEG.data,1) ;
EEG.nbchan   = r.header.nchannels ;
% EEG.nbchan       = size(EEG.data,1);
EEG.trials   = 1;
EEG.pnts     = r.ldnsamples ;
%size(EEG.data,2)

%EEG.pnts     = r.header.pnts 
%size(EEG.data,2);
EEG          = eeg_checkset(EEG, 'eventconsistency');
EEG          = eeg_checkset(EEG, 'makeur');

if ((size(EEG.data,1) ~= EEG.nbchan) && (size(EEG.data,2) ~= EEG.pnts))
   % Assume a data file
   EEG      = eeg_checkset(EEG, 'loaddata');
end
if length(options) > 2
    command = sprintf('EEG = pop_loadcnt(''%s'' %s);',fullFileName, options); 
else
    command = sprintf('EEG = pop_loadcnt(''%s'');',fullFileName); 
end;
return;
