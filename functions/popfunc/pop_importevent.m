% pop_importevent() - Import events in a dataset. If the EEG dataset
%              is the only inputs, a window pops up to ask for the relevant 
%              parameter values.
%
% Usage: >> [EEG,eventnumbers] = pop_importevent( EEG, 'key1', 'value1', ...);
%
% Input:
%   EEG      - input dataset
%
% Optional input
%  'append'   - ['yes'|'no'] 'yes' = Append to the current events of the
%               EEG dataset {default}: 'no' = Erase the previous events.
%  'event'    - [ 'filename'|array ] filename of a text file or
%               Matlab array of the global workspace containing an
%               array of events in the folowing format: The first column
%               is the type of the event, the second the latency and
%               the other ones are user defined fields. This function
%               can handle text entries in ascii files.
%  'fields'   - cell array of name for each user-defined column, 
%               optionally followed by a description. Ex: { 'type', 'latency' }
%  'skipline' - number of rows to skip for text files 
%  'indices'  - vector indicating the indices of the events to modify 
%  'timeunit' - [ latency unit in second ]. Default unit is 1 second. 
%  'delim'    - delimiting characters in file. Default is tab and space.
%  'align'    - [num] align the first event latency to existing 
%               event latency (number num) and check latency
%               consistency. Negative values can also be used and then
%               event number -num is aligned with the first pre-existing
%               event. Default is 0. (NaN-> no alignment)
%
% Outputs:
%   EEG          - dataset with updated event field
%   eventnumbers - indexes of the appended events
%
% Example: [EEG, eventnumbers] = pop_importevent(EEG, 'event', ...
%         'event_values.txt', 'fields', {'type', 'latency','condition' }, ...
%         'append', 'no', 'align', 0, 'timeunit', 1E-3 );
%
%         This loads the ascii file 'event_values.txt' containing 3 columns 
%         (event type, latency and condition). Latencies in the file are
%         in ms. The first event latency is re-aligned with the first 
%         pre-existing latencies in the dataset ('align', 0) and old
%         events are erased ('append' no).
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
		questdlg(['Though epoch information is defined in terms of event,' 10 ...
				  'this function is usually used to import events into continuous data.' 10 ...
				  'For data epochs you may better use menu /File/Import epoch info/'], ...
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
         { 'Style', 'edit', 'string', '' } { 'Style', 'text', 'string', 'Ex: type latency position', 'tooltipstring', helpfields } };
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
	    if results{2} == 0       , args = { args{:}, 'append', '''no''' }; end;
	    if ~isempty( results{3} ) 
	        if exist(results{3}) == 2,  args = { args{:}, 'event', [ '''' results{3} '''' ] }; 
	        else                        args = { args{:}, 'event', results{3} }; end; 
	    end;
	    if ~isempty( results{4} )  % FIELDS
	        args = { args{:}, 'fields', { parsetxt(results{4}) } };
	    end;
	    % handle skipline 
	    % ---------------     
	    if ~isempty(eval(results{end-2})), if eval(results{end-2}) ~= 0,  args = { args{:}, 'skipline', eval(results{end-2}) }; end; end;

	    % handle timeunit 
	    % -------------     
	    if ~isempty(eval(results{end-1})), if eval(results{end-1}) ~= 0,  args = { args{:}, 'timeunit', eval(results{end-1}) }; end; end;

	    % handle alignment 
	    % ----------------     
	    if ~isempty(eval(results{end})), if eval(results{end}) ~= 0,  args = { args{:}, 'align', eval(results{end}) }; end; end;
        args
else % no interactive inputs
    args = varargin;
    % scan args to modify array/file format
    % array are transformed into string 
    % files are transformed into string of string
    % (this is usefull to build the string command for the function)
    % --------------------------------------------------------------
    for index=1:2:length(args)
        if iscell(args{index+1}), args{index+1} = { args{index+1} }; end; % double nested 
        if isstr(args{index+1})                 args{index+1} = [ '''' args{index+1} '''' ]; % string 
        else if ~isempty( inputname(index+2) ), args{index+1} = inputname(index+2); end;
        end;
    end;                
end;
     
% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('Setevent: wrong syntax in function arguments'); return; end;
else
    g = [];
end;
g

% test the presence of variables
% ------------------------------
try, g.fields; 	 	  catch, g.fields = {}; end;
try, g.skipline;      catch, g.skipline = 0; end;
try, g.indices;  g.append = '''yes'''; catch, g.indices = []; end;
try, g.append; 	      catch, g.append = '''yes'''; end;
try, g.timeunit; 	  catch, g.timeunit = 1; end;
try, g.align; 	      catch, g.align = NaN; end;
try, g.delim; 	      catch, g.delim = char([9 32]); end;

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
        case 'yes', disp('Setevent warning: using align, events should not be appended but erased');
    end;
    if g.align.val < 0
        g.align.event = EEG.event(1).latency;
    else
        g.align.event = EEG.event(g.align.val+1).latency;
    end
    g.align.nbevent = length(EEG.event);
    g.align.txt = sprintf('Check alignment between pre-existing (old) and loaded event latencies:\nOld event latencies (10 first): %s ...\n', int2str(cell2mat({ EEG.event(1:min(10, length(EEG.event))).latency })));
    g.align
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
                case '''no'''
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
                case '''yes'''
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

% generate the output command
% ---------------------------
EEG = eeg_checkset(EEG, 'eventconsistency');
com = sprintf('%s = pop_importevent( %s', inputname(1), inputname(1));
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if isstr( args{i+1} ) com = sprintf('%s, ''%s'', %s', com, args{i}, args{i+1} );
        else    if ~iscell( args{i+1} ) com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
                else 
                    com = sprintf('%s, ''%s'', { ', com, args{i});
                    for index = 1:length( args{i+1}{1} ), com = sprintf('%s ''%s'' ', com, args{i+1}{1}{index}); end;
                    com = sprintf('%s }', com);
                end;                    
        end;
    else
        com = sprintf('%s, ''%s'', []', com, args{i} );
    end;
end;
com = [com ');'];

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline, delim );
    if varname(1) == '''', % mean that it is a filename
                          % --------------------------
        varname = eval(varname);
        array = loadtxt( varname, 'skipline', skipline );
        
    else % variable in the global workspace
         % --------------------------
         array = evalin('base', varname);
         if ~iscell(array)
             array = mat2cell(array, ones(1, size(array,1)), ones(1, size(array,2)));
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
