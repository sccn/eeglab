% pop_importevent() - Import events into an EEG dataset. If the EEG dataset
%              is the only input, a window pops up to ask for the relevant 
%              parameter values.
%
% Usage: >> EEG = pop_importevent( EEG ); % pop-up window mode
%        >> EEG = pop_importevent( EEG, 'key1', 'value1', ...);
%
% Graphic interface:
%   "Event indices" - [edit box] Enter indices of events to modify. 
%               Leave this field blank to import new events. 
%               Command line equivalent: 'indices'.
%   "Append events?" - [checkbox] Check this checkbox to clear prior
%               event information. In addition, see the "Align event latencies ..."
%               edit box. Command line equivalent: 'append'.
%   "Event file or array" - [edit box] Enter event file name. Use "Browse" 
%               button to browse for a file. If a file with the given name
%               cannot be found, the function search for a variable with
%               this name in the global workspace. 
%               Command line equivalent: 'filename'.
%   "Input field (column) name" - [edit box] Enter a name for each of the
%               columns in the event text file. If column names are defined 
%               in the text file, they cannnot be used and you must copy 
%               the names into this edit box (and skip the name row). Must
%               provide a name for each column. The keywords "type",
%               "latency", and "duration" are recognized EEGLAB keywords and 
%               should be used to define the event log file columns containing 
%               event types, latencies, and durations. Column names can be
%               separated by commas, quoted or not. 
%               Command line equivalent: fields.
%  "Latency time unit (sec)" - [edit box] Specify the time unit for the 
%               latency column defined above relative to seconds. 
%               Command line equivalent: 'timeunit'.
%  "Number of header lines to ignore" - [edit box] For some text files, the 
%               first rows do not contain epoch information and need to be
%               skipped. Command line equivalent: 'skiplines'.
%  "Align event latencies to data events" - [edit box] For most EEG datasets,
%               basic event information is defined along with the EEG, and
%               a more detailed file is recorded separately. This option 
%               helps fuse the two sources of information by realigning the 
%               imported data text file information into the existing event. 
%               A value of 0 indicates that the first events of the pre-defined 
%               events and imported events will be aligned. A positive value (num)
%               aligns the first event to the num-th pre-existing event. 
%               A negative value can also be used; then event number (-num) 
%               is aligned to the first pre-existing event.  Default is 0. 
%               (NaN-> no alignment). Command line equivalent is 'align'.
%  "Auto adjust event sampling rate" - [checkbox] When checked, the function
%               automatically adjusts the sampling rate of the new events so
%               they best align with the closest old events. This may account
%               for small differences in sampling rate that could lead to 
%               big differences at the end of the experiement (e.g., A 0.01%
%               clock difference over an hour would lead to a 360-ms difference 
%               if not corrected). Command line line equivalent is 'optimalim'.
% Input:
%   EEG      - input dataset
%
% Optional file or array input:
%  'event'    - [ 'filename'|array ] Filename of a text file, or name of s
%               Matlab array in the global workspace containing an
%               array of events in the folowing format: The first column
%               is the type of the event, the second the latency. 
%               The others are user-defined. The function can read 
%               either numeric or text entries in ascii files.
%  'fields'   - [Cell array] List of the name of each user-defined column, 
%               optionally followed by a description. Ex: { 'type', 'latency' }
%  'skipline' - [Interger] Number of header rows to skip in the text file 
%  'timeunit' - [ latency unit rel. to seconds ]. Default unit is 1 = seconds. 
%  'delim'    - [string] String of delimiting characters in the input file. 
%               Default is tab|space.
%
% Optional oldevent input:
%  'append'   - ['yes'|'no'] 'yes' = Append events to the current events in 
%               the EEG dataset {default}: 'no' = Erase the previous events.
%  'indices'  - {integer vector] Vector indicating the indices of the events to
%               modify. 
%  'align'    - [num] Align the first event latency to the latency of existing 
%               event number (num), and check latency consistency.
%  'optimalign' - ['on'|'off'] Optimize the sampling rate of the new events so 
%               they best align with old events. Default is 'on'.
%
% Outputs:
%   EEG          - EEG dataset with updated event fields
%
% Example: >> [EEG, eventnumbers] = pop_importevent(EEG, 'event', ...
%         'event_values.txt', 'fields', {'type', 'latency','condition' }, ...
%         'append', 'no', 'align', 0, 'timeunit', 1E-3 );
%
%         This loads the ascii file 'event_values.txt' containing 3 columns 
%         (event_type, latency, and condition). Latencies in the file are
%         in ms (1E-3). The first event latency is re-aligned with the 
%         beginning of the dataset ('align', 0). Any previous events 
%         are erased ('append', 'no').
%
% Author: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 9 Feb 2002
%
% See also: importevent(), pop_editeventfield(), pop_selectevent()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 9 Feb 2002, arno@salk.edu
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

function [EEG, com] = pop_importevent(EEG, varargin);

com ='';
if nargin < 1
   help pop_importevent;
   return;
end;	

if isempty(EEG.data)
   disp('pop_importevent(): error: cannot process empty dataset'); return;
end;    

I = [];

% warning if data epochs
% ----------------------
if nargin<2 && EEG.trials > 1
		questdlg2(strvcat('Though epoch information is defined in terms of the event structure,', ...
				  'this function is usually used to import events into continuous data.', ...
				  'For epoched data, use menu item ''File > Import epoch info'''), ...
				'pop_importevent warning', 'OK', 'OK');
end
	
% remove the event field
% ----------------------
if ~isempty(EEG.event), allfields = fieldnames(EEG.event);
else                    allfields = {}; end
    
if nargin<2
        commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
		helpfields = ['latency field (lowercase) must be present; field ''type'' and' 10 ...
				   '''duration'' are also recognized keywords and it is recommended to define them'];
        uilist = { ...
         { 'Style', 'text', 'string', 'Event indices', 'fontweight', 'bold' }, ...
         { 'Style', 'text', 'string', 'Append events?', 'fontweight', 'bold' } };
        geometry    = { [ 1 1.1 1.1 1] [ 1 1 2] [1 1 2] [ 1.2 1 1] };
        uilist = { uilist{:}, ...     
         { 'Style', 'text', 'string', 'Event file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''globfile'';' commandload ] }, ...
         { 'Style', 'edit' } ...
         { 'Style', 'checkbox', 'string', 'Yes/No', 'value', 0 }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { }, { 'Style', 'text', 'string', 'NB: No = overwrite', 'value', 0 }, { }, ...
         { 'Style', 'text', 'string', 'Input field (column) names       ', 'fontweight', 'bold', 'tooltipstring', helpfields } ...
         { 'Style', 'edit', 'string', '' } { 'Style', 'text', 'string', 'Ex: type latency duration', 'tooltipstring', helpfields } };
         geometry = { geometry{:} [1.2 1 1] [1.2 1 1] [1.2 1 1] [1.2 0.2 1.8] };
         uilist   = { uilist{:}, ...
                { 'Style', 'text', 'string', 'Number of file header lines', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '0' }, ...
					  { 'Style', 'text', 'string', '(latency field required above)', 'tooltipstring', helpfields },...
                { 'Style', 'text', 'string', 'Time unit (sec)', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '1' } ...
					  { 'Style', 'text', 'string', 'Ex: If ms, 1E-3; if points, NaN' },...
                { 'Style', 'text', 'string', 'Align event latencies to data events', 'horizontalalignment', 'left' }, ...
					  { 'Style', 'edit', 'string', fastif(isempty(EEG.event),'NaN','0') } { 'Style', 'text', 'string', 'See Help' },...
                { 'Style', 'text', 'string', 'Auto adjust new events sampling rate', 'horizontalalignment', 'left' }, ...
					  { 'Style', 'checkbox', 'value' 1 } { },...
               };
        results = inputgui( geometry, uilist, 'pophelp(''pop_importevent'');', 'Import event info -- pop_importevent()' );
        if length(results) == 0, return; end

	    % decode top inputs
	    % -----------------
	    args = {};
	    if ~isempty( results{1} ), args = { args{:}, 'indices', eval( [ '[' results{1} ']' ]) }; end
	    if results{2} == 0 && ~isempty(EEG.event), args = { args{:}, 'append', 'no' }; end
	    if ~isempty( results{3} ), 
            if ischar( results{3} ) && ~exist(results{3})
                args = { args{:}, 'event', evalin('base', results{3}) }; 
            else
                args = { args{:}, 'event', results{3} }; 
            end
        end
	    if ~isempty( results{4} ), args = { args{:}, 'fields', parsetxt(results{4}) }; end
        
	    % handle skipline 
	    % ---------------     
	    if ~isempty(eval(results{end-3})), if eval(results{end-3}) ~= 0,  args = { args{:}, 'skipline', eval(results{end-3}) }; end; end

	    % handle timeunit 
	    % -------------     
	    if ~isempty(eval(results{end-2})), if eval(results{end-2}) ~= 0,  args = { args{:}, 'timeunit', eval(results{end-2}) }; end; end

	    % handle alignment 
	    % ----------------     
	    if ~isempty(eval(results{end-1})), if ~isnan(eval(results{end-1})),  args = { args{:}, 'align', eval(results{end-1}) }; end; end
	    if ~results{end} ~= 0,  args = { args{:}, 'optimalign', 'off' }; end
        
else % no interactive inputs
    args = varargin;
    % scan args to modify array/file format
    % array are transformed into string 
    % files are transformed into string of string
    % (this is usefull to build the string command for the function)
    % --------------------------------------------------------------
    for index=1:2:length(args)
        if iscell(args{index+1}), if iscell(args{index+1}{1}) args{index+1} = args{index+1}{1}; end; end; % double nested 
        if ischar(args{index+1}) && length(args{index+1}) > 2 && args{index+1}(1) == '''' && args{index+1}(end) == ''''             
            args{index+1} = args{index+1}(2:end-1); end
        %else if ~isempty( inputname(index+2) ), args{index+1} = inputname(index+2); end
        %end
    end;                
end

EEG.event = importevent( [], EEG.event, EEG.srate, args{:});

% generate ur variables
% ---------------------
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% generate the output command
% ---------------------------
com = sprintf('EEG = pop_importevent( EEG, %s);', vararg2str(args));
