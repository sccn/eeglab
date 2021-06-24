% pop_editeventfield() - Add/remove/rename/modify a field in the event structure 
%              of an EEG dataset. Can also be used to append new events to the end of the 
%              event structure or to delete all current events. If the dataset is 
%              the only input, a window pops up to ask for relevant parameter values.
%
% Usage: >> [EEG] = pop_editeventfield( EEG, 'key1', 'value1', ...);
%
% Input:
%   EEG      - input dataset
%
% Optional inputs:
%  'FIELDNAME_X' - [ 'filename'|vector ]. Name of a current or new 
%               user-defined event field. The ascii file, vector variable, 
%               or explicit numeric vector should contain values for this field 
%               for all events specified in 'indices' (below) or in new events
%               appended to the dataset (if 'indices' not specified). If one arg 
%               value is given, it will be used for all the specified events.
%               If the arg is [], the named field is *removed* from all the 
%               specified (or new) events. Use this option to add a new field to 
%               all events in the dataset, or to modify the field values of specified 
%               events, to specify field information for new events appended to
%               the event structure, or to remove an event field. For example,
%               'FIELDNAME_X' may be 'latency', 'type', 'duration'.
%  'latency'  - [ 'filename'|vector ] example of field name (see description above).
%  'type'     - [ 'filename'|vector ] example of field name (see description above).
%  'duration' - [ 'filename'|vector ] example of field name (see description above). 
%               Units must be in seconds (s).
%  'FIELDNAME_X_info' - new comment string for field FIELDNAME_X.
%  'latency_info'     - description string for the latency field.
%  'type_info'        - description string for the type field.
%  'duration_info'    - description string for the duration field.
%  'indices'  - [vector of event indices] The indices of the events to modify. 
%               If adding a new event field, events not listed here 
%               will have an empty field value IF they are not in an epoch
%               containing events whose field value is specified. However,
%               if adding a new (FIELDNAME) field to an epoched dataset,
%               and the field value for only one event in some data epoch is 
%               specified, then the other events in the same epoch will be given 
%               the specified value. If field values of more than one, but not all 
%               the events in an epoch are specified, then unspecified events at 
%               the beginning of the epoch will be given the value of the first 
%               specified event, unspecified epoch events after this event will 
%               be given the field value of the second specified epoch event, etc.
%               {default|[]: modify all events in the dataset}
%  'rename'   - ['FIELDNAME1->FIELDNAME2'] rename field 'FIELDNAME1' to
%               field 'FIELDNAME2'. Ex: { 'rename', 'blocktype->condition' }.    
%  'delold'   - ['yes'|'no'] 'yes' = delete ALL previous events.  {default: 'no'}
%  'timeunit' - [latency field time unit in fraction of seconds]. Ex: 1e-3 -> 
%               read specified latencies as msec {default: 1 (-->seconds)}
%  'skipline' - number of leading text file lines to skip in named text files.
%               {default: 0}.
%  'delim'    - delimiting characters in named text files {default: tabs and spaces}
%
% Outputs:
%   EEG          - dataset with updated event field
%
% Example: EEG = pop_editeventfield(EEG, 'type', 1, ...
%                         'sleepstage', 'sleepstage_values.txt', ...
%                         'latency', [100 130 123 400],'timeunit',1e-3);
%          % Append 4 events to the EEG struct, all of type 1, at the latencies
%          % given (in msec) with user-defined field ('sleepstage') values read 
%          % from a one-column ascii file ('sleepstage_values.txt').
%
% Example: EEG = pop_editeventfield(EEG, 'indices', 1:2:3277, ...
%                          'sleepstage', 'sleepstage_values.txt');
%          % Add a new 'sleepstage' field to all events in the EEG struct.
%          % Read 'sleepstage' field values for odd-numbered events from a one-column 
%          % ascii file ('sleepstage_values.txt') -- even-numbered % events will have 
%          % an empty 'sleepstage' value (unless the data are epoched, see 'indices' above).
%
% Note: To save events into a readable table, first convert them to 'array' format,
%        then save the array.
%            >> events = eeg_eventformat(EEG.event, 'array');
%            >> save -ascii myevents.txt events 
%
% Author: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 9 Feb 2002-
%
% See also: pop_importevent(), eeg_eventformat(), pop_selectevent()

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

%02/13/2001 fix bug if EEG.event is empty -ad
%03/12/2001 add timeunit option -ad
%03/18/2001 debug rename option -ad & sm
%03/18/2001 correct allignment problem -ad & ja

function [EEG, com] = pop_editeventfield(EEG, varargin);

com ='';
if nargin < 1
   help pop_editeventfield;
   return;
end;	

if isempty(EEG.data)
    error('Setevent error: cannot process empty dataset');
end;    

I = [];

% remove the event field
% ----------------------
if ~isempty(EEG.event), 
    allfields = fieldnames(EEG.event); 
    ind1 = strmatch('urevent', allfields, 'exact');
    ind2 = strmatch('epoch', allfields, 'exact');
    allfields([ind1 ind2]) = [];
    try, EEG.eventdescription([ind1 ind2]) = []; catch, end
else
    allfields = { 'type' 'latency' }; 
end
    
if nargin<2
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    uilist = { ...
        { 'Style', 'text', 'string', 'Edit fields:', 'fontweight', 'bold'  }, ...
        { 'Style', 'text', 'string', 'Edit description', 'fontweight', 'bold'  }, ...
        { 'Style', 'text', 'string', sprintf('New values (file or array containing %d values)', length(EEG.event)), 'fontweight', 'bold'  }, ...
        { 'Style', 'text', 'string', 'Delete field', 'fontweight', 'bold'  } ...
             };
    geometry = { [1.05 1.05 2 0.8] };

    listboxtext = { 'No field selected' };
    txt_warn = 'warndlg2(strvcat(''Warning: deleting/renaming this field might cause EEGLAB'', ''to be unstable. Some functionalities will also be lost.''));';
    cb_warn2 = [ 'strtmp = get(gcbo, ''string''); if ~isempty(strmatch(strtmp(get(gcbo, ''value'')), { ''latency'' ''type''}, ''exact'')),' txt_warn 'end;' ];
    for index = 1:length(allfields) 
        geometry = { geometry{:} [1 1 1 0.7 0.2 0.32 0.2] };
        description = '';
        try, 
            description = fastif(isempty(EEG.eventdescription{index}), '', EEG.eventdescription{index});
            description = description(1,:);
            tmplines = find(description == 10);
            if ~isempty(tmplines), description = description(1:tmplines(1)-1); end
        catch, end
        if strcmp(allfields{index}, 'latency') || strcmp(allfields{index}, 'type')
             cb_warn  = { 'callback' [ 'if get(gcbo, ''value''),' txt_warn 'end;' ] };
        else cb_warn = { };
        end
        if strcmp(allfields{index}, 'latency')
            tmpfield = [ allfields{index} '(s)' ];
        elseif strcmp(allfields{index}, 'duration')
            tmpfield = [ allfields{index} '(s)' ];
        else
            tmpfield = allfields{index};
        end
        uilist   = { uilist{:}, ...
                     { 'Style', 'text', 'string', tmpfield }, ...
                     { 'Style', 'pushbutton', 'string', description, 'callback', ...
                       [ 'tmpuserdata = get(gcf, ''userdata'');' ...
                         'tmpuserdata{' int2str(index) '} = pop_comments(tmpuserdata{' int2str(index) ...
                         '}, ''Comments on event field: ' allfields{index} ''');' ...
                         'set(gcbo, ''string'', tmpuserdata{' int2str(index) '});' ...
                         'set(gcf, ''userdata'', tmpuserdata); clear tmpuserdata;' ] }, ...
                     { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  allfields{index} }, ...
                     { 'Style', 'pushbutton', 'string', 'Browse', 'callback', ['tagtest = ''' allfields{index} ''';' commandload ] }, ...
                     { }, { 'Style', 'checkbox', 'string', '    ', cb_warn{:} }, { } };
        listboxtext = { listboxtext{:}  allfields{index} }; 
    end
    index = length(allfields) + 1;
    uilist   = { uilist{:}, ...
                 { 'Style', 'edit', 'string', ''}, ...
                 { 'Style', 'pushbutton', 'string', '', 'callback', ...
                   [ 'tmpuserdata = get(gcf, ''userdata'');' ...
                     'tmpuserdata{' int2str(index) '} = pop_comments(tmpuserdata{' int2str(index) ...
                     '}, ''Comments on new event field:'');' ...
                     'set(gcbo, ''string'', tmpuserdata{' int2str(index) '});' ...
                     'set(gcf, ''userdata'', tmpuserdata); clear tmpuserdata;' ] }, ...
                 { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'newfield' }, ...
                 { 'Style', 'pushbutton', 'string', 'Browse', 'callback', ['tagtest = ''newfield'';' commandload ] }, ...
                 { 'Style', 'text', 'string', '-> add field'} ...
                 { } ...
                 { 'Style', 'text', 'string', 'Rename field', 'fontweight', 'bold' }, ...
                 { 'Style', 'popupmenu', 'string', listboxtext 'callback' cb_warn2 }, ...
                 { 'Style', 'text', 'string', 'as', 'fontweight', 'bold' }, ...
                 { 'Style', 'edit', 'string', '' } ...
                 fastif(isunix,{ 'Style', 'text', 'string', '(Click on field name to select it!)' },{ })};
    geometry = { geometry{:} [1 1 1 0.7 0.72] [1] [1 1.2 0.6 1 2] };

    descriptions = EEG.eventdescription;
    if isempty(descriptions), descriptions = { '' '' }; end
    [results userdat ]= inputgui( geometry, uilist, 'pophelp(''pop_editeventfield'');', ...
                                  'Edit event field(s) -- pop_editeventfield()', { descriptions{:} '' } );
    if length(results) == 0, return; end

    % decode top inputs
    % -----------------
    args = { };
    
    % dealing with existing fields
    %-----------------------------
    for index = 1:length(allfields) 
        if results{index*2} == 1, args = { args{:}, allfields{index}, [] };
        else 
            if ~isempty( results{index*2-1} )
                if exist(results{index*2-1}) == 2,  args = { args{:}, allfields{index}, [ results{index*2-1} ] }; % file
                else                                args = { args{:}, allfields{index}, results{index*2-1} }; end
            end
            try, 
                if ~strcmp( userdat{index}, EEG.eventdescription{index})
                    args = { args{:}, [ allfields{index} 'info' ], userdat{index} }; 
                end
            catch, end
        end;     
    end
    
    % dealing with the new field
    %---------------------------
    sub = 3;
    if ~isempty( results{end-sub} )
        args = { args{:}, results{end-sub}, results{end-sub+1} }; % file
    end
    
    % handle rename 
    % -------------
    if results{end-1} ~= 1, args = { args{:}, 'rename', [ allfields{results{end-1}-1} '->' results{end} ] }; end;  
else % no interactive inputs
    args = varargin;
    % scan args to modify array/file format
    % array are transformed into string 
    % files are transformed into string of string
    % (this is usefull to build the string command for the function)
    % --------------------------------------------------------------
    for index=1:2:length(args)
        if iscell(args{index+1}), args{index+1} = { args{index+1} }; end; % double nested 
        if ischar(args{index+1})   args{index+1} = args{index+1}; % string 
        end
    end;                
end

% create structure
% ----------------
try, g = struct(args{:});
catch, disp('pop_editeventfield(): wrong syntax in function arguments'); return; 
end

% test the presence of variables
% ------------------------------
try, g.skipline;      catch, g.skipline = 0; end
try, g.indices;       catch, g.indices = [1:length(EEG.event)]; end
try, g.delold; 	      catch, g.delold = 'no'; end
try, g.timeunit; 	  catch, g.timeunit = 1; end
try, g.delim; 	      catch, g.delim = char([9 32]); end
if ischar(g.indices), g.indices = eval([ '[' g.indices ']' ]); end

tmpfields = fieldnames(g);
% scan all the fields of g
% ------------------------
for curfield = tmpfields'
    if ~isempty(EEG.event), allfields = fieldnames(EEG.event);
    else                    allfields = {}; end
    switch lower(curfield{1})
       case { 'append' 'delold', 'fields', 'skipline', 'indices', 'timeunit', 'delim' }, ; % do nothing now
       case 'rename',
            if isempty( findstr('->',g.rename) ), 
                disp('warning pop_editeventfield() bad syntax for ''rename'', ignoring input'); 
            else
                oldname = g.rename(1:findstr('->',g.rename)-1);
                newname = g.rename(findstr('->',g.rename)+2:end);
                indexmatch = strmatch(oldname, allfields);
                if isempty(indexmatch), disp('pop_editeventfield() warning: name not found for rename'); 
                else
                    for index  = 1:length(EEG.event)
                        eval([ 'EEG.event(index).' newname '=EEG.event(index).' oldname ';']);  
                    end;    
                    EEG.event = rmfield(EEG.event, oldname);
                end
                if isfield(EEG, 'urevent')
                    disp('pop_editeventfield() warning: field name not renamed in urevent structure');
                end
            end
       otherwise, % user defined field command
                  % --------------------------
            infofield = findstr(curfield{1}, 'info');
            if ~isempty(infofield) && infofield == length( curfield{1} )-3
                % description of a field
                % ----------------------     
                fieldname = curfield{1}(1:infofield-1);
                indexmatch = strmatch( fieldname, allfields);
                if isempty( indexmatch )
                    disp(['pop_editeventfield() warning: Field ' fieldname ' not found to add description, ignoring']);
                else
                    EEG.eventdescription{indexmatch} = getfield(g, curfield{1});
                end
            else              
                % not an field for description
                % ----------------------------      
	            if isempty( getfield(g, curfield{1}) ) % delete
	                 indexmatch = strmatch( curfield{1}, allfields);
                     if isempty( indexmatch )
                        disp(['pop_editeventfield() warning: Field ''' curfield{1} ''' not found for deletion, ignoring']);
                     else
	                    EEG.event = rmfield(EEG.event, curfield{1}); 
	                    allfields(indexmatch) = [];
                        if isfield(EEG, 'urevent')
                            fprintf('pop_editeventfield() warning: field ''%s'' not deleted from urevent structure\n', curfield{1}  );
                        end
						try,
							EEG.eventdescription(indexmatch) = [];
						catch, end
	                 end;    
	            else % interpret
		            switch g.delold % delete old events
		                case 'yes'
		                      EEG.event = load_file_or_array( getfield(g, curfield{1}), g.skipline, g.delim );
		                      allfields = { curfield{1} };
                              EEG.event = eeg_eventformat(EEG.event, 'struct', allfields);
                              EEG.event = recomputelatency( EEG.event, 1:length(EEG.event), EEG.srate, g.timeunit);
                              EEG = eeg_checkset(EEG, 'makeur');
		                 case 'no' % match existing fields
		                           % ---------------------
		                      tmparray = load_file_or_array( getfield(g, curfield{1}), g.skipline, g.delim );
		                      if isempty(g.indices) g.indices = [1:size(tmparray(:),1)] + length(EEG.event); end
		                      
		                      indexmatch = strmatch(curfield{1}, allfields);
		                      if isempty(indexmatch) % no match
		                          disp(['pop_editeventfield(): creating new field ''' curfield{1} '''' ]);
		                      end
                              try
                                  EEG.event = setstruct(EEG.event, curfield{1}, g.indices, tmparray);
                              catch,
                                  error('Wrong size for input array');
                              end
  							  if strcmp(curfield{1}, 'latency')
								  EEG.event = recomputelatency( EEG.event, g.indices, EEG.srate, g.timeunit);
							  end
 							  if strcmp(curfield{1}, 'duration')
                                  for indtmp = 1:length(EEG.event)
                                      EEG.event(indtmp).duration = EEG.event(indtmp).duration*EEG.srate;
                                  end
							  end
                              if isfield(EEG, 'urevent')
                                  disp('pop-editeventfield(): updating urevent structure');
                                  try
                                      tmpevent = EEG.event;
                                      for indtmp = g.indices(:)'
                                          if ~isempty(EEG.event(indtmp).urevent)
                                              tmpval      = getfield (EEG.event, {indtmp}, curfield{1});
                                              EEG.urevent = setfield (EEG.urevent, { tmpevent(indtmp).urevent }, ...
                                                                                   curfield{1}, tmpval);
                                          end
                                      end
                                  catch,
                                      disp('pop_editeventfield(): problem while updating urevent structure');
                                  end
                              end
		             end
	            end
	        end;    
      end
end

if isempty(EEG.event) % usefull 0xNB empty structure
    EEG.event = [];
end
EEG = eeg_checkset(EEG, 'eventconsistency');

% generate the output command
% ---------------------------
com = sprintf('EEG = pop_editeventfield( EEG, %s);', vararg2str(args));

% interpret the variable name
% ---------------------------
function str = str2str( array )
	str = '';
	for index = 1:size(array,1)
		str = [ str ', ''' array(index,:) '''' ];
	end
	if size(array,1) > 1
		str = [ 'strvcat(' str(2:end) ')'];
	else
		str = str(2:end);
	end;	
return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline, delim );
    if ischar(varname)
        if exist(varname) == 2 % mean that it is a filename
                               % --------------------------
            array = loadtxt( varname, 'skipline', skipline, 'delim', delim);
            
        else % variable in the global workspace
             % --------------------------
             try
                array = evalin('base', varname);
             catch
                 array = { varname };
             end
             if ~iscell(array)
                 array = mattocell(array, ones(1, size(array,1)), ones(1, size(array,2)));
             end
        end
    else
        array = mattocell(varname);
    end
return;

% update latency values
% ---------------------
function event = recomputelatency( event, indices, srate, timeunit);
    if ~isfield(event, 'latency'), return; end
    for index = indices
        event(index).latency = event(index).latency*srate*timeunit+1;
    end
         
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
    if exist('indices') ~= 1, indices = 1:length(var); end
    if iscell(values)
        if length(values) > 1
            for index = 1:length(indices)
                var = setfield(var, {indices(index)}, fieldname, values{index});
            end
        else
            for index = 1:length(indices)
                var = setfield(var, {indices(index)}, fieldname, values{1});
            end
        end; 
    else % Code below may be unused
        if length(values) > 1
            for index = 1:length(indices)
                var = setfield(var, {indices(index)}, fieldname, values(index));
            end
        else
            for index = 1:length(indices)
                var = setfield(var, {indices(index)}, fieldname, values);
            end
        end
    end
return;
