% pop_importevent() - Import events into an EEG dataset. If the EEG dataset
%              is the only input, a window pops up to ask for the relevant 
%              parameter values.
%
% Usage: >> EEG = pop_importevent( EEG ); % pop-up window mode
%        >> EEG = pop_importevent( EEG, 'key1', 'value1', ...);
%
% Graphic interface:
%   "Event indices" - [edit box] Enter indices of event to modify. 
%               Leave this field blank to import new events. 
%               Command line equivalent: 'indices'.
%   "Append events?" - [checkbox] Check this checkbox to clear prior
%               event information. See also the "Align event latencies ..."
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
%               "latency", and "duration" are recognized EEGLAB keyword and 
%               should be used to define the columns containing event types, 
%               event latencies, and event durations. Columns names can be
%               separated by comas, quoted or not. 
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
%               event number (num), and check latency consistency. See also the 
%               GUI help above.
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

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 9 Feb 2002, arno@salk.edu
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

% $Log: pop_importevent.m,v 
% Revision 1.30  2004/01/31 01:23:24  arn
% nothin

% Revision 1.29  2004/01/31 01:16:53  arn
% message chang
%
% Revision 1.28  2003/12/12 00:26:55  arno
% debuging oldevents and align
%
% Revision 1.27  2003/12/11 19:42:45  arno
% debuging auto alignment
%
% Revision 1.26  2003/12/11 02:54:10  arno
% automatic alignment
%
% Revision 1.25  2003/11/19 19:28:46  arno
% nothing
%
% Revision 1.24  2003/11/07 02:14:26  arno
% more detailed message
%
% Revision 1.23  2003/11/07 02:10:56  arno
% remove events with no latencies
%
% Revision 1.22  2003/11/07 01:29:25  arno
% [Anothing
%
% Revision 1.21  2003/11/04 01:11:14  arno
% fixing reading presentation file problem
%
% Revision 1.20  2003/11/04 00:10:49  arno
% g.delim integer OK
%
% Revision 1.19  2003/11/03 23:23:37  arno
% adapting to new history ...
%
% Revision 1.18  2003/11/01 03:05:34  arno
% removing some double brackets
%
% Revision 1.17  2003/11/01 02:58:46  arno
% implementing finputcheck
%
% Revision 1.16  2003/06/18 22:24:40  arno
% ur events
%
% Revision 1.15  2003/06/09 16:53:28  arno
% nothing
%
% Revision 1.14  2003/04/10 17:34:53  arno
% header edit
%
% Revision 1.13  2003/01/24 02:14:12  arno
% debuging 'delim' option
%
% Revision 1.12  2003/01/03 02:20:47  scott
% header edits -sm
%
% Revision 1.11  2003/01/03 01:40:58  arno
% changing exemple
%
% Revision 1.10  2002/11/18 18:18:01  arno
% removing debugging message
%
% Revision 1.9  2002/10/10 21:18:09  arno
% remove text message
%
% Revision 1.8  2002/09/04 17:53:45  luca
% added ;
%
% Revision 1.7  2002/08/12 18:33:09  arno
% quesdlg2
%
% Revision 1.6  2002/08/06 21:49:46  arno
% spelling
%
% Revision 1.5  2002/06/28 02:27:32  arno
% adding ori fields
%
% Revision 1.4  2002/04/18 18:25:27  arno
% typo can not
%
% Revision 1.3  2002/04/18 00:56:36  arno
% inserting warning for data epochs
%
% Revision 1.2  2002/04/10 03:25:28  arno
% added eeg check set consistency
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

%02/13/2001 fix bug if EEG.event is empty -ad
%03/12/2001 add timeunit option -ad
%03/18/2001 debug rename option -ad & sm
%03/18/2001 correct allignment problem -ad & ja

function [EEG, com] = pop_importevent(EEG, varargin);

com ='';
if nargin < 1
   help pop_importevent;
   return;
end;	

if isempty(EEG.data)
    disp('Setevent error: cannot process empty dataset'); return;
end;    

I = [];

% warning if data epochs
% ----------------------
if nargin<2 & EEG.trials > 1
		questdlg2(strvcat('Though epoch information is defined in terms of event,', ...
				  'this function is usually used to import events into continuous data.', ...
				  'For data epochs you may better use menu /File/Import epoch info/'), ...
				'pop_importevent warning', 'OK', 'OK');
end;
	
% remove the event field
% ----------------------
if ~isempty(EEG.event), allfields = fieldnames(EEG.event);
else                    allfields = {}; end;
    
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
					  { 'Style', 'text', 'string', 'Ex: If ms, 1E-3' },...
                { 'Style', 'text', 'string', 'Align event latencies to data events', 'horizontalalignment', 'left' }, ...
					  { 'Style', 'edit', 'string', fastif(isempty(EEG.event),'NaN','0') } { 'Style', 'text', 'string', 'See Help' },...
                { 'Style', 'text', 'string', 'Auto adjust new events sampling rate', 'horizontalalignment', 'left' }, ...
					  { 'Style', 'checkbox', 'value' 1 } { },...
               };
        results = inputgui( geometry, uilist, 'pophelp(''pop_importevent'');', 'Import event info -- pop_importevent()' );
        if length(results) == 0, return; end;

	    % decode top inputs
	    % -----------------
	    args = {};
	    if ~isempty( results{1} ), args = { args{:}, 'indices', eval( [ '[' results{1} ']' ]) }; end;
	    if results{2} == 0       , args = { args{:}, 'append', 'no' }; end;
	    if ~isempty( results{3} ), args = { args{:}, 'event', results{3} }; end; 
	    if ~isempty( results{4} ), args = { args{:}, 'fields', parsetxt(results{4}) }; end;
        
	    % handle skipline 
	    % ---------------     
	    if ~isempty(eval(results{end-3})), if eval(results{end-3}) ~= 0,  args = { args{:}, 'skipline', eval(results{end-3}) }; end; end;

	    % handle timeunit 
	    % -------------     
	    if ~isempty(eval(results{end-2})), if eval(results{end-2}) ~= 0,  args = { args{:}, 'timeunit', eval(results{end-2}) }; end; end;

	    % handle alignment 
	    % ----------------     
	    if ~isempty(eval(results{end-1})), if ~isnan(eval(results{end-1})),  args = { args{:}, 'align', eval(results{end-1}) }; end; end;
	    if ~results{end} ~= 0,  args = { args{:}, 'optimalign', 'off' }; end;
        
else % no interactive inputs
    args = varargin;
    % scan args to modify array/file format
    % array are transformed into string 
    % files are transformed into string of string
    % (this is usefull to build the string command for the function)
    % --------------------------------------------------------------
    for index=1:2:length(args)
        if iscell(args{index+1}), if iscell(args{index+1}{1}) args{index+1} = args{index+1}{1}; end; end; % double nested 
        if isstr(args{index+1}) & length(args{index+1}) > 2 & args{index+1}(1) == '''' & args{index+1}(end) == ''''             
            args{index+1} = args{index+1}(2:end-1); end;
        %else if ~isempty( inputname(index+2) ), args{index+1} = inputname(index+2); end;
        %end;
    end;                
end;

EEG.event = importevent( [], EEG.event, EEG.srate, args{:});

% generate ur variables
% ---------------------
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% generate the output command
% ---------------------------
com = sprintf('%s = pop_importevent( %s, %s);', inputname(1), inputname(1), vararg2str(args));
