% pop_editeventfield() - Load/remove event fields in a dataset. If the 
%              EEG dataset is the only inputs, a window pops up to ask for 
%              the relevant parameter values.
%
% Usage: >> [EEG] = pop_editeventfield( EEG, 'key1', 'value1', ...);
%
% Input:
%   EEG      - input dataset
%
% Optional input
%  'latency'  - [ 'filename'|array|number ], 1-column file or
%               array containing the latency of the events. If the arg
%               is a single number, it will be used for all events.
%  'type'     - [ 'filename'|array|number ], Type (numbers) of the events.
%               Same format as 'latency' {Default event type is 1}.
%  'USER_VAR' - [ 'filename'|array|number|[] ], here 'USER_VAR'
%               is the name of a user-defined field. If the argument
%               is [], the field is removed from all the events.
%  'skipline' - number of rows to skip for text files 
%  'rename'   - ['USER_VAR1->USER_VAR2'] rename field 'USER_VAR1' into
%               'USER_VAR2'. Ex: { 'rename', 'type->condition' }.    
%  'append'   - ['yes'|'no'] 'yes' = Append to the current events of the
%               EEG dataset {default}: 'no' = Erase the previous events.
%  'indices'  - vector indicating the indices of the events to modify 
%  'timeunit' - [ latency unit in second ]. Default unit is 1 second. 
%  'delim'    - delimiting characters in file. Default is tab and space.
%  'latencyinfo'  - comment string for the latency field
%  'typeyinfo'    - comment string for the type field
%  'USER_VARinfo' - comment string for the variable USERVAR
%
% Outputs:
%   EEG          - dataset with updated event field
%
% Example: [EEG, eventnumbers] = pop_editeventfield(EEG, 'type', 1, ...
%                         'sleepstage', 'sleepstage_values.txt', ...
%                         'latency', [0.100 0.130 0.123 0.400] );
%
%         This appends 4 events to the EEG struct, all of type 1,
%         with user-defined field ('sleepstage') values read from an
%         ascii 1-column file ('sleepstage_values.txt').
%
% Notes: As indicated above, to remove a user=defined field, simply
%           enter an empty array [].
%            Ex: >> EEG = pop_editeventfield(EEG,'sleepstage',[]);
%        To save the events first convert them in 'array' format and
%        then save the array.
%            >> tmp = eeg_eventformat(EEG.event, 'array');
%            >> save -ascii filename.txt tmp 
%
% Author: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 9 Feb 2002
%
% See also: pop_importevent(), eeg_eventformat(), pop_selectevent()

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
% Revision 1.18  2002/08/20 23:38:12  arno
% can not remove epoch field
%
% Revision 1.17  2002/08/14 16:32:34  arno
% constrain description to 1 line
%
% Revision 1.16  2002/08/13 00:03:10  scott
% text
%
% Revision 1.15  2002/08/13 00:01:34  scott
% text
%
% Revision 1.14  2002/08/12 23:59:00  scott
% same
%
% Revision 1.13  2002/08/12 23:58:11  scott
% YES/no button text
%
% Revision 1.12  2002/08/12 22:49:33  arno
% text
%
% Revision 1.11  2002/08/12 21:48:46  arno
% same
%
% Revision 1.10  2002/08/12 21:48:26  arno
% text
%
% Revision 1.9  2002/08/12 21:47:14  arno
% text
%
% Revision 1.8  2002/07/25 01:08:15  arno
% debugging
%
% Revision 1.7  2002/07/23 18:37:41  arno
% removing debug message
%
% Revision 1.6  2002/05/03 22:11:07  arno
% debugging rename
%
% Revision 1.5  2002/04/18 18:23:14  arno
% typo can not
%
% Revision 1.4  2002/04/09 03:51:32  arno
% modifying input comments and history
%
% Revision 1.3  2002/04/09 01:41:49  arno
% debuging array passage, change default
%
% Revision 1.2  2002/04/08 21:55:09  arno
% not edited
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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
    disp('Setevent error: cannot process empty dataset'); return;
end;    

I = [];

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
        uilist = { ...
         { 'Style', 'text', 'string', 'Event indices to replace', 'fontweight', 'bold' }, ...
         { 'Style', 'text', 'string', 'Append events?', 'fontweight', 'bold' } };
        geometry    = { [ 1 1.1 0.2 1 0.2] [ 1.2 1 1.7 0.3] };
        uilist = { uilist{:} {} {} {} ...
         { 'Style', 'edit', 'string', ['1:' int2str(length(EEG.event))] } ...
         { 'Style', 'checkbox', 'string', 'YES / no', 'value', 1 }, ...
         { 'Style', 'text', 'string', 'NB: No (unchecked) -> overwrite events' }, { }, { }, ...   
         { 'Style', 'text', 'string', 'Edit fields:', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Field description', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Values file|array', 'fontweight', 'bold'  }, ...
         {   }, ...
         { 'Style', 'text', 'string', 'Delete field', 'fontweight', 'bold'  } ...
         };
        geometry = { geometry{:} [1] [0.9 0.95 1.2 0.8 1] };

	    listboxtext = 'No field selected';  
	    for index = 1:length(allfields) 
	        geometry = { geometry{:} [1 1 1 1 0.45 0.35 0.45] };
	        description = '';
	        try, 
				description = fastif(isempty(EEG.eventdescription{index}), '', EEG.eventdescription{index});
				description = description(1,:);
				tmplines = find(description == 10);
				if ~isempty(tmplines), description = description(1:tmplines(1)-1); end;
	        catch, end;
	        uilist   = { uilist{:}, ...
	         { 'Style', 'text', 'string', allfields{index} }, ...
	         { 'Style', 'pushbutton', 'string', description, 'callback', ...
			 [ 'tmpuserdata = get(gcf, ''userdata'');' ...
			   'tmpuserdata{' int2str(index) '} = pop_comments(tmpuserdata{' int2str(index) '}, ''Comments on event field: ' allfields{index} ''');' ...
			   'set(gcbo, ''string'', tmpuserdata{' int2str(index) '});' ...
			   'set(gcf, ''userdata'', tmpuserdata); clear tmpuserdata;' ] }, ...
	         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  allfields{index} }, ...
	         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', ['tagtest = ''' allfields{index} ''';' commandload ] }, ...
	         { }, fastif(strcmp(allfields{index}, 'epoch'), {}, { 'Style', 'checkbox', 'string', '    ' }),{ } };
	         listboxtext = [ listboxtext '|' allfields{index} ]; 
	    end;
	    geometry = { geometry{:} [1 1 1 1 0.45 0.35 0.45] [1] [1 2 0.6 1.3 1.5] };
	    uilist   = { uilist{:}, ...
	         { 'Style', 'edit', 'string', ''}, ...
	         { 'Style', 'edit', 'string', '' }, ...
	         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'newfield' }, ...
	         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', ['tagtest = ''newfield'';' commandload ] }, ...
	         { }, { },{ } };
         uilist   = { uilist{:}, ...
                { } ...
                { 'Style', 'text', 'string', 'Rename', 'fontweight', 'bold' }, ...
                { 'Style', 'listbox', 'string', listboxtext }, ...
                { 'Style', 'text', 'string', 'as', 'fontweight', 'bold' }, ...
                { 'Style', 'edit', 'string', '' } ...
                fastif(isunix,{ 'Style', 'text', 'string', 'Click on selected field!','fontweight','bold' },{ })};

        [results userdat ]= inputgui( geometry, uilist, 'pophelp(''pop_editeventfield'');', ...
									  'Edit event field(s) -- pop_editeventfield()', EEG.eventdescription );
        if length(results) == 0, return; end;

	    % decode top inputs
	    % -----------------
	    args = {};
	    if ~isempty( results{1} ), args = { args{:}, 'indices', results{1} }; end;
	    if results{2} == 0       , args = { args{:}, 'append', 'no' }; end;
	    
	    % dealing with existing fields
	    %-----------------------------
	    for index = 1:length(allfields) 
	        if results{index*2+2} == 1, args = { args{:}, allfields{index}, [] };
	        else 
				if ~isempty( results{index*2+1} )
	                if exist(results{index*2+1}) == 2,  args = { args{:}, allfields{index}, [ results{index*2+1} ] }; % file
	                else                                args = { args{:}, allfields{index}, results{index*2+1} }; end;
				end;
				try, 
					if ~strcmp( userdat{index}, EEG.eventdescription{index})
						args = { args{:}, [ allfields{index} 'info' ], userdat{index} }; 
					end;
				catch, end;
	        end;     
	    end;
	    
	    % dealing with the new field
	    %---------------------------
	    sub = 4;
	    if ~isempty( results{end-sub} )
	        if ~isempty( results{end-sub+2} )
	            if exist(results{end-sub+2}) == 2,  args = { args{:}, results{end-sub}, [ results{end-sub+2} ] }; % file
	            else                                args = { args{:}, results{end-sub}, results{end-sub+2} }; end;
	            if ~isempty( results{end-sub+1} ),  args = { args{:}, [ results{end-sub} 'info' ], [ results{end-sub+1} ] }; end;
	        else
	            disp(['The new field' results{end-sub} ' was ignored since no input data were given for it.' ]);
	        end;
	    end;  
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
        if isstr(args{index+1})                 args{index+1} = args{index+1}; % string 
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

% test the presence of variables
% ------------------------------
try, g.skipline;      catch, g.skipline = 0; end;
try, g.indices;  g.append = 'yes'; catch, g.indices = []; end;
try, g.append; 	      catch, g.append = 'yes'; end;
try, g.timeunit; 	  catch, g.timeunit = 1; end;
try, g.align; 	      catch, g.align = NaN; end;
try, g.delim; 	      catch, g.delim = char([9 32]); end;
g.align.val = g.align;
if isstr(g.indices), g.indices = eval([ '[' g.indices ']' ]); end;

tmpfields = fieldnames(g);
% scan all the fields of g
% ------------------------
for curfield = tmpfields'
    if ~isempty(EEG.event), allfields = fieldnames(EEG.event);
    else                    allfields = {}; end;
    switch lower(curfield{1})
       case {'append', 'fields', 'skipline', 'indices', 'timeunit', 'align', 'delim' }, ; % do nothing now
       case 'rename',
            if isempty( findstr('->',g.rename) ), disp('Set warning: bad syntax for rename'); end;
            oldname = g.rename(1:findstr('->',g.rename)-1);
            newname = g.rename(findstr('->',g.rename)+2:end);
            indexmatch = strmatch(oldname, allfields);
            if isempty(indexmatch), disp('Set warning: name not found for rename'); 
            else
                for index  = 1:length(EEG.event)
                     eval([ 'EEG.event(index).' newname '=EEG.event(index).' oldname ';']);  
                end;    
                EEG.event = rmfield(EEG.event, oldname);
            end;
       otherwise, % user defined field command
                  % --------------------------
            infofield = findstr(curfield{1}, 'info');
            if ~isempty(infofield) & infofield == length( curfield{1} )-3
                % description of a field
                % ----------------------     
                fieldname = curfield{1}(1:infofield-1);
                indexmatch = strmatch( fieldname, allfields);
                if isempty( indexmatch )
                    disp(['Setevent warning: Field ' fieldname ' not found to add description, ignoring']);
                else
                    EEG.eventdescription{indexmatch} = getfield(g, curfield{1});
                end;
            else              
                % not an field for description
                % ----------------------------      
	            if isempty( getfield(g, curfield{1}) ) % delete
	                 indexmatch = strmatch( curfield{1}, allfields);
                     if isempty( indexmatch )
                        disp(['Set warning: Field ''' curfield{1} ''' not found for deletion, ignoring']);
                     else
	                    EEG.event = rmfield(EEG.event, curfield{1}); 
	                    allfields(indexmatch) = [];
						try,
							EEG.eventdescription(indexmatch) = [];
						catch, end;
	                 end;    
	            else % interpret
		            switch g.append 
		                case 'no'
		                      EEG.event = load_file_or_array( getfield(g, curfield{1}), g.skipline, g.delim );
		                      allfields = { curfield{1} };
                              EEG.event = eeg_eventformat(EEG.event, 'struct', allfields);
                              EEG.event = recomputelatency( EEG.event, 1:length(EEG.event), EEG.srate, g.timeunit, g.align);
		                 case 'yes' % match existing fields
		                            % ---------------------
		                      tmparray = load_file_or_array( getfield(g, curfield{1}), g.skipline, g.delim );
		                      if isempty(g.indices) g.indices = [1:size(tmparray(:),1)] + length(EEG.event); end;
		                      
		                      indexmatch = strmatch(curfield{1}, allfields);
		                      if isempty(indexmatch) % no match
		                          disp(['Set: field ''' curfield{1} ''' not found, creating new field']);
		                      end;
		                      EEG.event = setstruct(EEG.event, curfield{1}, g.indices, cell2mat(tmparray));     
							  if strcmp(curfield{1}, 'latency')
								  EEG.event = recomputelatency( EEG.event, g.indices, EEG.srate, g.timeunit, g.align);
							  end;
		             end;
	            end;
	        end;    
      end;
end;

if isempty(EEG.event) % usefull 0xNB empty structure
    EEG.event = [];
end;
EEG = eeg_checkset(EEG, 'eventconsistency');

% generate the output command
% ---------------------------
com = sprintf('%s = pop_editeventfield( %s', inputname(1), inputname(1));
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if isstr( args{i+1} ) com = sprintf('%s, ''%s'', %s', com, args{i}, str2str(args{i+1}) );
        else    
			if ~iscell( args{i+1} ) com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
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
function str = str2str( array )
	str = '';
	for index = 1:size(array,1)
		str = [ str ', ''' array(index,:) '''' ];
	end;
	if size(array,1) > 1
		str = [ 'strvcat(' str(2:end) ')'];
	else
		str = str(2:end);
	end;	
return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline, delim );
	if exist(varname) == 2 % mean that it is a filename
                           % --------------------------
        array = loadtxt( varname );
        
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
