% pop_importevent() - Import events into an EEG dataset. If the EEG dataset
%              is the only input, a window pops up to ask for the relevant 
%              parameter values.
%
% Usage: >> EEG = pop_importevent( EEG ); % pop-up window mode
%        >> [EEG,eventindices] = pop_importevent( EEG, 'key1', 'value1', ...);
%
% Graphic interface:
%   "Event indices" - [edit box] enter indices of event to modify. Let
%               this field blank to import new events. Command line
%               equivalent: 'indices'.
%   "Append events?" - [checkbox] check this checkbox to clear prior
%               event information. See also the "Align event latencies ..."
%               edit box. Command line equivalent: 'append'.
%   "Event file or array" - [edit box] enter event file name. Use "Browse" 
%               button to browse for a file. If a file with the given name
%               can not be found, the function search for a variable with
%               this name in the global workspace. Command line 
%               equivalent: 'filename'.
%   "Input field (column) name" - [edit box] enter a name for each of the
%               column in the event text file. If columns names are defined 
%               in the text file, they cannnot be used and you must copy 
%               their name in this edit box (and skip the rows). One column 
%               name for each column must be provided. The keywords "type"
%               and "latency" should be used to define the columns containing
%               event types and event latencies. Columns names can be
%               separated by comas, quoted or not. Command line 
%               equivalent: fields.
%  "Latency time unit (sec)" - [edit box] specify the time unit for the 
%               latency column defined above relative to seconds. Command
%               line equivalent: 'timeunit'.
%  "Number of header lines to ignore" - [edit box] for some text files, the 
%               first rows do not contain epoch information and have to be
%               skipped. Command line equivalent: 'skiplines'.
%  "Align event latencies to data events" - [edit box] For most EEG datasets,
%               basic event information are defined along with the EEG and
%               a more detailed file is recorded separatelly. This option 
%               help fuse the two information by realigning the imported
%               data text file to the existing event. A value of 0 indicates
%               that the first events of the pre-defined events and imported
%               events will be aligned. A positive value align the first 
%               event to the number num pre-existing event. A negative value 
%               can also be used; then event number -num is aligned with the 
%               first pre-existing event. Default is 0. (NaN-> no alignment).
% Input:
%   EEG      - input dataset
%
% Optional input
%  'append'   - ['yes'|'no'] 'yes' = Append events to the current events in 
%               the EEG dataset {default}: 'no' = Erase the previous events.
%  'event'    - [ 'filename'|array ] Filename of a text file, or name of s
%               Matlab array in the global workspace containing an
%               array of events in the folowing format: The first column
%               is the type of the event, the second the latency. 
%               The others are user-defined. The function can read 
%               either numeric or text entries in ascii files.
%  'fields'   - [Cell array] list of the name of each user-defined column, 
%               optionally followed by a description. Ex: { 'type', 'latency' }
%  'skipline' - [Interger] Number of header rows to skip in the text file 
%  'indices'  - {integer vector] Vector indicating the indices of the events to
%               modify. 
%  'timeunit' - [ latency unit rel. to seconds ]. Default unit is 1 = seconds. 
%  'delim'    - [string] String of delimiting characters in the input file. 
%               Default is tab|space.
%  'align'    - [num] Align the first event latency to the latency of existing 
%               event number num, and check latency consistency. See also GUI
%               help above.
% Outputs:
%   EEG          - EEG dataset with updated event fields
%   eventindices - Indexes of the appended events
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
% See also: pop_editeventfield(), eeg_eventformat(), pop_selectevent()

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

% $Log: not supported by cvs2svn $
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
		helpfields = ['latency field (lowercase) must be present' 10 ...
				   'field type is also a recognized keyword and it is recommended to define it'];
        uilist = { ...
         { 'Style', 'text', 'string', 'Event indices', 'fontweight', 'bold' }, ...
         { 'Style', 'text', 'string', 'Append events?', 'fontweight', 'bold' } };
        geometry    = { [ 1 1.1 1.1 1] [ 1 1 2] [1 1 2] [ 2 1 1] };
        uilist = { uilist{:}, ...     
         { 'Style', 'text', 'string', 'Event file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''globfile'';' commandload ] }, ...
         { 'Style', 'edit' } ...
         { 'Style', 'checkbox', 'string', 'Yes/No', 'value', 0 }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { }, { 'Style', 'text', 'string', 'NB: No = overwrite', 'value', 0 }, { }, ...
         { 'Style', 'text', 'string', 'Input field (column) names       ', 'fontweight', 'bold', 'tooltipstring', helpfields } ...
         { 'Style', 'edit', 'string', '' } { 'Style', 'text', 'string', 'Ex: latency type position', 'tooltipstring', helpfields } };
         geometry = { geometry{:} [2 1 1] [2 1 1] [2 1 1] };
         uilist   = { uilist{:}, ...
                { 'Style', 'text', 'string', 'Number of file header lines', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '0' }, ...
					  { 'Style', 'text', 'string', 'Note: latency required', 'tooltipstring', helpfields },...
                { 'Style', 'text', 'string', 'Latency time unit (sec)', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '1' } ...
					  { 'Style', 'text', 'string', 'Ex: If ms, 1E-3' },...
                { 'Style', 'text', 'string', 'Align event latencies to data events', 'horizontalalignment', 'left' }, ...
					  { 'Style', 'edit', 'string', fastif(isempty(EEG.event),'NaN','0') } { 'Style', 'text', 'string', 'See Help' },...
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
	    if ~isempty(eval(results{end-2})), if eval(results{end-2}) ~= 0,  args = { args{:}, 'skipline', eval(results{end-2}) }; end; end;

	    % handle timeunit 
	    % -------------     
	    if ~isempty(eval(results{end-1})), if eval(results{end-1}) ~= 0,  args = { args{:}, 'timeunit', eval(results{end-1}) }; end; end;

	    % handle alignment 
	    % ----------------     
	    if ~isempty(eval(results{end})), if eval(results{end}) ~= 0,  args = { args{:}, 'align', eval(results{end}) }; end; end;
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

g = finputcheck( args, { 'fields'    'cell'     []         {};
                         'skipline'  'integer'  [0 Inf]    0;
                         'indices'   'integer'  [1 Inf]    [];
                         'append'    'string'   {'yes' 'no' '''yes''' '''no''' }         'yes';
                         'timeunit'  'real'     [0 Inf]    1;
                         'event'     { 'real' 'string' }     []    [];
                         'align'     'integer'  []         NaN;
                         'delim'     {'integer' 'string'}   []         char([9 32])}, 'pop_importevent');
if isstr(g), error(g); end;
if ~isempty(g.indices), g.append = 'yes'; end;
g.delim = char(g.delim);    

% test the presence of variables
% ------------------------------
%try, g.fields; 	 	  catch, g.fields = {}; end;
%try, g.skipline;      catch, g.skipline = 0; end;
%try, g.indices;  g.append = '''yes'''; catch, g.indices = []; end;
%try, g.append; 	      catch, g.append = '''yes'''; end;
%try, g.timeunit; 	  catch, g.timeunit = 1; end;
%try, g.align; 	      catch, g.align = NaN; end;
%try, g.delim; 	      catch, g.delim = char([9 32]); end;

% determine latency for old event alignment
% -----------------------------------------
g.align.val = g.align;
if ~isnan(g.align.val)
    if isempty(EEG.event)
        error('Setevent: no pre-existing event, cannot perform alignment');
    end;    
    if ~isfield(EEG.event, 'latency')
        error('Setevent: pre-existing events do not have a latency field for re-alignment');
    end;    
    switch g.append
        case {'yes' '''yes'''}, disp('Setevent warning: using align, events should not be appended but erased');
    end;
    if g.align.val < 0
        g.align.event = EEG.event(1).latency;
    else
        g.align.event = EEG.event(g.align.val+1).latency;
    end
    g.align.nbevent = length(EEG.event);
    g.align.txt = sprintf('Check alignment between pre-existing (old) and loaded event latencies:\nOld event latencies (10 first): %s ...\n', int2str(cell2mat({ EEG.event(1:min(10, length(EEG.event))).latency })));
end;

tmpfields = fieldnames(g);
% scan all the fields of g
% ------------------------
for curfield = tmpfields'
    if ~isempty(EEG.event), allfields = fieldnames(EEG.event);
    else                    allfields = {}; end;
    switch lower(curfield{1})
        case {'append', 'fields', 'skipline', 'indices', 'timeunit', 'align', 'delim' }, ; % do nothing now
        case 'event', % load an ascii file
            switch g.append 
                case { '''no''' 'no' } % for backward compatibility
                      EEG.event = load_file_or_array( g.event, g.skipline, g.delim );
                      allfields = g.fields(1:min(length(g.fields), size(EEG.event,2)));
                      EEG.event = eeg_eventformat(EEG.event, 'struct', allfields);
					  % generate ori fields
					  % -------------------
					  for index = 1:length(EEG.event)
						  EEG.event(index).init_index = index;
						  EEG.event(index).init_time  = EEG.event(index).latency*g.timeunit;
					  end;
					  EEG.event = recomputelatency( EEG.event, 1:length(EEG.event), EEG.srate, g.timeunit, g.align);
                case { '''yes''' 'yes' }
                      % match existing fields
                      % ---------------------
                      tmparray = load_file_or_array( g.event, g.skipline, g.delim );
                      if isempty(g.indices) g.indices = [1:size(tmparray,1)] + length(EEG.event); end;
                      if length(g.indices) ~= size(tmparray,1)
                            error('Set error: number of row in file does not match the number of event given as input'); 
                      end;

                      % add field
                      % ---------
                      g.fields = getnewfields( g.fields, size(tmparray,2)-length(g.fields));
                      
                      % add new values
                      % ---------------------
                      for eventfield = 1:size(tmparray,2)
                          EEG.event = setstruct( EEG.event, g.fields{eventfield}, g.indices, cell2mat(tmparray(:,eventfield)));
                      end;      
					  % generate ori fields
					  % -------------------
					  offset = length(EEG.event)-size(tmparray,2);
					  for index = 1:size(tmparray,2)
						  EEG.event(index+offset).init_index = index;
						  EEG.event(index+offset).init_time  = EEG.event(index+offset).latency*g.timeunit;
					  end;
                      EEG.event = recomputelatency( EEG.event, g.indices, EEG.srate, g.timeunit, g.align);
            end;
      end;
end;

if isempty(EEG.event) % usefull 0xNB empty structure
    EEG.event = [];
end;

% remove the events which latency are out of boundary
% ---------------------------------------------------
if isfield(EEG.event, 'latency')
	alllatencies = cell2mat( { EEG.event.latency } );
	I1 = find(alllatencies < 0);
	I2 = find(alllatencies > EEG.pnts*EEG.trials);
	if (length(I1) + length(I2)) > 0 
	    fprintf('Setevent warning: %d/%d events had out-of-bounds latencies and were removed\n', length(I1) + length(I2), length(EEG.event));
	    EEG.event(union(I1, I2)) = [];
	end;
end;

% generate ur variables
% ---------------------
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% generate the output command
% ---------------------------
com = sprintf('%s = pop_importevent( %s, %s);', inputname(1), inputname(1), vararg2str(args));

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline, delim );
    if isstr(varname) & exist(varname) == 2  % mean that it is a filename
                                             % --------------------------
        array = loadtxt( varname, 'skipline', skipline, 'delim', delim );
        
    else % variable in the global workspace
         % --------------------------
         if isstr(varname)
             array = evalin('base', varname);
             if ~iscell(array)
                 array = mat2cell(array, ones(1, size(array,1)), ones(1, size(array,2)));
             end;    
         else
             array = varname;
         end;
    end;     
return;

% update latency values
% ---------------------
function event = recomputelatency( event, indices, srate, timeunit, align);
    if ~isfield(event, 'latency'), return; end;
    for index = indices
        event(index).latency = round(event(index).latency*srate*timeunit);
    end;
    if ~isnan( align.val )
        if align.val >= 0, alignlatency = event(1).latency;
        else               alignlatency = event(-align.val+1).latency;
        end;
        for index = indices
             event(index).latency = event(index).latency-alignlatency+align.event;
        end;
        if length(event) ~= align.nbevent
            disp([ 'Setevent warning: the number of pre-existing events do not correspond to the ' ...
                   'number of event that were read, so their latencies may have been wrongly re-aligned' ]);
        end;           
        fprintf(align.txt);
        fprintf('New event latencies (10 first): %s ...\n', int2str(cell2mat({ event(1:min(10, length(event))).latency })));
    end;
         
% create new field names
% ----------------------
function epochfield = getnewfields( epochfield, nbfields )
   count = 1;
   while nbfields > 0
       if isempty( strmatch([ 'var' int2str(count) ], epochfield ) )
               epochfield =  { epochfield{:} [ 'var' int2str(count) ] };
               nbfields = nbfields-1;
       else    count = count+1;
       end;                    
   end;     
return;
    
function var = setstruct( var, fieldname, indices, values )
    if exist('indices') ~= 1, indices = 1:length(var); end;
    for index = 1:length(indices)
        var = setfield(var, {indices(index)}, fieldname, values(index));
    end;
return;
