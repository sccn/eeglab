% pop_importepoch() - Export epoch and/or epoch event information to the event 
%                     structure array of an EEG dataset. If the dataset is 
%                     the only input, a window pops up to ask for the relevant 
%                     parameter values.
% Usage:
%   >> EEGOUT = pop_importepoch( EEG ); % pop-up window mode
%   >> EEGOUT = pop_importepoch( EEG, filename, fieldlist, 'key', 'val', ...);
%
% Graphic interface:
%  "Epoch file or array" - [edit box] enter epoch text file name. Use "Browse" 
%                      button to browse for a file. If a file with the given name
%                      can not be found, the function search for a variable with
%                      this name in the global workspace. Command line 
%                      equivalent: filename.
%  "File input field ..." - [edit box] enter a name for each of the column in the
%                      text file. If columns names are defined in the text file,
%                      they cannnot be used and you must copy their names
%                      in this edit box (and skip the rows). One column name
%                      for each column must be provided. The keywords "type" and
%                      "latency" should not be used. Columns names can be
%                      separated by comas, quoted or not. Command line 
%                      equivalent: fieldlist.
%  "Field name(s) containing event latencies" - [edit box] enter columns name(s)
%                      containing latency information. It is not necessary to 
%                      define a latency field for epoch information. All fields 
%                      that contain latencies will be imported as different event 
%                      types. For instance, if field 'RT' contains latencies, 
%                      events of type 'RT' will be created with latencies given 
%                      in the RT field. See notes. Command line 
%                      equivalent: 'latencyfields'.
%  "Field name(s) containing event durations" - [edit box] enter columns name(s)
%                      containing duration information. It is not necessary to 
%                      define a latency field for epoch information, but if you
%                      do, a duration field (or 0) must be entered for each 
%                      latency field you define. For instance if the latency fields
%                      are "'rt1' 'rt2'", then you must have duration fields
%                      such as "'dr1' 'dr2'". If duration is not defined for event 
%                      latency 'tr1', you may enter "0 'rt2'". Command line 
%                      equivalent: 'durationfields'.
%  "Field name containing time locking event type(s)" - [edit box] if one column
%                      contain the epoch type, its name must be defined in the 
%                      previous edit box and copied here. It is not necessary to 
%                      define a type field for the time-locking event (TLE). By 
%                      default it is defined as type ''TLE'' at time 0 for all 
%                      epochs. Command line equivalent: 'typefield'.
%  "Latency time unit rel. to seconds" - [edit box] specify the time unit for 
%                      latency columns defined above. Command line 
%                      equivalent: 'timeunit'.
%  "Number of header lines to ignore" - [edit box] for some text files, the first
%                      rows do not contain epoch information and have to be
%                      skipped. Command line equivalent: 'headerlines'.
%  "Remove old epoch and event info" - [checkbox] check this checkbox
%                      to remove any prior event or epoch information. Command
%                      line equivalent: 'clearevents'.
%
% Inputs:
%   EEG              - Input EEG dataset
%   filename         - Name of an ascii file with epoch and/or epoch event 
%                      information organised in columns. ELSE, name of a Matlab
%                      variable with the same information (either a Matlab array 
%                      or cell array). 
%   fieldlist        - {cell array} Label of each column (data field) in the file.
%
% Optional inputs:
%   'typefield'      - ['string'] Name of the field containing the type(s)
%                      of the epoch time-locking events (at time 0). 
%                      By default, all the time-locking events are assigned 
%                      type 'TLE' (for "time-locking event"). 
%   'latencyfields'  - {cell array} Field names that contain the latency 
%                      of an event. These fields are transferred into 
%                      events whose type will be the same as the name of
%                      the latency field. (Ex: field RT -> type 'RT' events).
%   'durationfields'  - {cell array} Field names that contain the duration 
%                      of an event. 
%   'timeunit'       - [float] Optional unit for latencies relative to seconds. 
%                      Ex: sec -> 1, msec -> 1e-3. Default: Assume latencies 
%                      are in time points (relative to the time-zero time point 
%                      in the epoch). 
%   'headerlines'    - [int] Number of header lines in the input file to ignore. 
%                      {Default 0}.
%   'clearevents'    - ['on'|'off'], 'on'-> clear the old event array. 
%                      {Default 'on'}
%
% Output:
%   EEGOUT - EEG dataset with modified event structure
%
% FAQ:
% 1) Why is this function so complex? This function can handle as many events
%    per epochs as needed, and the information is stored in terms of events
%    rather than epoch information, which requires some conversion.
% 2) Can I access epoch information later? The epoch information is stored in
%    "EEG.event" and the information is stored in terms of events only. For 
%    user convenience the "EEG.epoch" structure is generated automatically
%    from the event structure. See EEGLAB manual for more information.
%
% Authors: Arnaud Delorme & Scott Makeig, CNL/Salk Institute, 11 March 2002
%
% See also: eeglab()
 
% Copyright (C) 15 Feb 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% graphic interface INFOS
% 03/18/02 debugging variable passing - ad & lf
% 03/18/02 adding event updates and incremental calls -ad
% 03/25/02 adding default event description -ad
% 03/28/02 fixed latency calculation -ad

function [EEG, com] = pop_importepoch( EEG, filename, fieldlist, varargin);
    
com ='';
if nargin < 1
    help pop_importepoch
    return;
end;
if nargin < 2
    geometry    = { [ 1 1 1.86] [1] [1 0.66] [2.5 1 0.6] [2.5 1 0.6] [2.5 1 0.6] [1] [2 0.5 0.5] [2 0.5 0.5] [2 0.17 0.86]};
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    helpstrtype = ['It is not necessary to define a type field for the time-locking event.' 10 ...
			   'By default it is defined as type ''TLE'' at time 0 for all epochs'];
    helpstrlat  = ['It is not necessary to define a latency field for epoch information.' 10 ...
			   'All fields that contain latencies will be imported as different event types.' 10 ...
			   'For instance, if field ''RT'' contains latencies, events of type ''RT''' 10 ...
                           'will be created with latencies given in the RT field'];
    helpstrdur  = ['It is not necessary to define a duration for each event (default is 0).' 10 ...
			   'However if a duration is defined, a corresponding latency must be defined too' 10 ...
               '(in the edit box above). For each latency field, you have define a duration field.' 10 ...
               'if no duration field is defined for a specific event latency, enter ''0'' in place of the duration field' ];
	uilist = { ...
         { 'Style', 'text', 'string', 'Epoch file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''globfile'';' commandload ] }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { }...
         { 'Style', 'text', 'string', 'File input field (col.) names', 'fontweight', 'bold' }, { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', '           Field name(s) containing event latencies', 'horizontalalignment', 'right', ...
           'fontweight', 'bold', 'tooltipstring', helpstrlat },  ...
  		 { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', '(Ex: RT)', 'tooltipstring', helpstrlat }, ...
         { 'Style', 'text', 'string', '           Field name(s) containing event durations', 'horizontalalignment', 'right', ...
           'fontweight', 'bold', 'tooltipstring', helpstrdur },  ...
  		 { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', 'NOTE', 'tooltipstring', helpstrdur }, ...
         { 'Style', 'text', 'string', '           Field name containing time-locking event type(s)', 'horizontalalignment', 'right', ...
                                      'tooltipstring', helpstrtype },  ...
  		 { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', 'NOTE', 'tooltipstring', helpstrtype }, ...
         { } ...
         { 'Style', 'text', 'string', 'Latency time unit rel. to seconds. Ex: ms -> 1E-3 (NaN -> samples)', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '1' }, { } ...         
         { 'Style', 'text', 'string', 'Number of file header lines to ignore', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '0' }, { },...        
         { 'Style', 'text', 'string', 'Remove old epoch and event info (set = yes)', 'horizontalalignment', 'left' }, { 'Style', 'checkbox', 'value', isempty(EEG.event) }, { } };         
    result = inputgui( geometry, uilist, 'pophelp(''pop_importepoch'');', 'Import epoch info (data epochs only) -- pop_importepoch()');
    if length(result) == 0, return; end;

    filename    = result{1};
    fieldlist   = parsetxt( result{2} );
    options = {};
    if ~isempty( result{3}), options = { options{:} 'latencyfields' parsetxt( result{3} ) }; end; 
    if ~isempty( result{4}), options = { options{:} 'durationfields' parsetxt( result{4} ) }; end; 
    if ~isempty( result{5}), options = { options{:} 'typefield' result{5} }; end; 
    if ~isempty( result{6}), options = { options{:} 'timeunit' eval(result{6}) }; end; 
    if ~isempty( result{7}), options = { options{:} 'headerlines' eval(result{7}) }; end; 
    if ~result{8}, options = { options{:} 'clearevents' 'off'}; end; 
else 
    if ~isempty(varargin) & ~isstr(varargin{1})
        % old call compatibility
        options = { 'latencyfields' varargin{1} };
        if nargin > 4
            options = { options{:} 'timeunit' varargin{2} }; 
        end; 
        if nargin > 5
            options = { options{:} 'headerlines' varargin{3} }; 
        end; 
        if nargin > 6
            options = { options{:} 'clearevents' fastif(varargin{4}, 'on', 'off') }; 
        end; 
    else
        options = varargin;
    end;
end;

g = finputcheck( options, { 'typefield'      'string'   []       ''; ...
                            'latencyfields'  'cell'     []       {}; ...
                            'durationfields' 'cell'     []       {}; ...
                            'timeunit'       'real'     [0 Inf]  1/EEG.srate; ...
                            'headerlines'    'integer'  [0 Inf]  0; ...
                            'clearevents'    'string'   {'on';'off'}  'on'}, 'pop_importepoch');
if isstr(g), error(g); end;

% check duration field
% --------------------
if ~isempty(g.durationfields)
    if length(g.durationfields) ~= length(g.latencyfields) 
        error( [ 'If duration field(s) are defined, their must be as many duration' 10 ...
              'fields as there are latency fields (or enter 0 instead of a field for no duration' ]);
    end;
else
    for index = 1:length(g.latencyfields) 
        g.durationfields{index} = 0;
    end;
end;

% convert filename
% ----------------
fprintf('Pop_importepoch: Loading file or array...\n');
if isstr(filename)
	% check filename
	% --------------
	if exist(filename) == 2 & evalin('base', ['exist(''' filename ''')']) == 1
		disp('Pop_importepoch WARNING: FILE AND ARRAY WITH THE SAME NAME, LOADING FILE');
	end;
    values = load_file_or_array( filename, g.headerlines );
else
    values = filename;
    filename = inputname(2);
end;

% check parameters
% ----------------
if size(values,1) < size(values,2), values = values'; end;
if length(fieldlist) ~= size(values,2)
    values = values';
    if length(fieldlist) ~= size(values,2)
        error('There must be as many field names as there are columsn in the file/array');
    end;
end;
if ~iscell(fieldlist)
    otherfieldlist = { fieldlist };
    fieldlist = { fieldlist };
end;
otherfieldlist = setdiff( fieldlist, g.latencyfields);
otherfieldlist = setdiff( otherfieldlist, g.typefield);
for index = 1:length(g.durationfields)
    if isstr(g.durationfields{index})
        otherfieldlist = setdiff( otherfieldlist, g.durationfields{index});
    end;
end;
if size(values,1) ~= EEG.trials
    error( [ 'Pop_importepoch() error: the number of rows in the input file/array does' 10 ... 
             'not match the number of trials. Maybe you forgot to specify the file header length?' ]);
end;    

% create epoch array info
% -----------------------
if iscell( values )
    for indexfield = 1:length(fieldlist)
        for index=1:EEG.trials
            eval( ['EEG.epoch(index).' fieldlist{ indexfield } '=values{ index, indexfield };'] );
        end;
    end;    
else
    for indexfield = 1:length(fieldlist)
        for index=1:EEG.trials
            eval( ['EEG.epoch(index).' fieldlist{ indexfield } '=values( index, indexfield);'] );
        end;
    end;    
end;

if isempty( EEG.epoch )
    error('Pop_importepoch: cannot process empty epoch structure');
end;
epochfield = fieldnames( EEG.epoch );

% determine the name of the non latency fields
% --------------------------------------------
tmpfieldname = {};
for index = 1:length(otherfieldlist)
    if isempty(strmatch( otherfieldlist{index}, epochfield ))
         error(['Pop_importepoch: field ''' otherfieldlist{index} ''' not found']);
    end;
    switch otherfieldlist{index}
       case {'type' 'latency'}, tmpfieldname{index} = [ 'epoch' otherfieldlist{index} ];
       otherwise,               tmpfieldname{index} = otherfieldlist{index};
    end;   
end;

if ~isempty(EEG.event)
    if ~isfield(EEG.event, 'epoch')
        g.clearevents = 'on';
        disp('Pop_importepoch: cannot add events to a non-epoch event structure, erasing old epoch structure');
    end;
end;
if strcmpi(g.clearevents, 'on')
    if ~isempty(EEG.event)
        fprintf('Pop_importepoch: deleting old events if any\n');
    end;
    EEG.event = [];
else 
    fprintf('Pop_importepoch: appending new events to the existing event array\n');
end;
           
% add time locking event fields
% -----------------------------
if EEG.xmin <= 0
    fprintf('Pop_importepoch: adding automatically Time Locking Event (TLE) events\n');
    if ~isempty(g.typefield)
        if isempty(strmatch( g.typefield, epochfield )) 
            error(['Pop_importepoch: type field ''' g.typefield ''' not found']);
        end;
    end;
    for trial = 1:EEG.trials
        EEG.event(end+1).epoch = trial; 
        if ~isempty(g.typefield)
            eval( ['EEG.event(end).type = EEG.epoch(trial).' g.typefield ';'] );
        else 
            EEG.event(end).type = 'TLE';
        end;
        EEG.event(end).latency  = -EEG.xmin*EEG.srate+1+(trial-1)*EEG.pnts;
        EEG.event(end).duration = 0;
    end;
end;

% add latency fields
% ------------------
for index = 1:length(g.latencyfields)
    if isempty(strmatch( g.latencyfields{index}, epochfield )) 
         error(['Pop_importepoch: latency field ''' g.latencyfields{index} ''' not found']);
    end;
    for trials = 1:EEG.trials
        EEG.event(end+1).epoch  = trials; 
        EEG.event(end).type     = g.latencyfields{index};
        EEG.event(end).latency  = (getfield(EEG.epoch(trials), g.latencyfields{index})*g.timeunit-EEG.xmin)*EEG.srate+1+(trials-1)*EEG.pnts;
        if g.durationfields{index} ~= 0 & g.durationfields{index} ~= '0'
            EEG.event(end).duration = getfield(EEG.epoch(trials), g.durationfields{index})*g.timeunit*EEG.srate;
        else
            EEG.event(end).duration = 0;
        end;
    end;
end;

% add non latency fields
% ----------------------
if ~isfield(EEG.event, 'epoch') % no events added yet
    for trial = 1:EEG.trials
        EEG.event(end+1).epoch = trial;
    end;
end;
for indexevent = 1:length(EEG.event)
    if ~isempty( EEG.event(indexevent).epoch )
        for index2 = 1:length(tmpfieldname)
            eval( ['EEG.event(indexevent).' tmpfieldname{index2} ' = EEG.epoch(EEG.event(indexevent).epoch).' otherfieldlist{index2} ';' ] );
    	end;
    end;
end;

% adding desciption to the fields
% -------------------------------
if ~isfield(EEG, 'eventdescription' ) | isempty( EEG.eventdescription )
	allfields = fieldnames(EEG.event);
    EEG.eventdescription{strmatch('epoch', allfields, 'exact')} = 'Epoch number';
	if ~isempty(strmatch('type', allfields)), EEG.eventdescription{strmatch('type', allfields)} = 'Event type'; end;
	if ~isempty(strmatch('latency', allfields)), EEG.eventdescription{strmatch('latency', allfields)} = 'Event latency'; end;
	if ~isempty(strmatch('duration', allfields)), EEG.eventdescription{strmatch('duration', allfields)} = 'Event duration'; end;
end;

% checking and updating events
% ----------------------------
EEG = pop_editeventvals( EEG, 'sort', { 'epoch', 0 } ); % resort fields
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% generate the output command
% ---------------------------
if isempty(filename) & nargout == 2
    disp('Pop_importepoch: cannot generate command string'); return;
else 
	com = sprintf('%s = pop_importepoch( %s, ''%s'', %s);', inputname(1), inputname(1), ...
                  filename, vararg2str( { fieldlist options{:} }));
end;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline );

    if exist( varname ) == 2
        if exist(varname) ~= 2, error( [ 'Set error: no filename ' varname ] ); end;

		fid=fopen(varname,'r','ieee-le');
		if fid<0, error( ['Set error: file ''' varname ''' found but error while opening file'] ); end;  

		for index=1:skipline	fgetl(fid); end; % skip lines ---------
        inputline = fgetl(fid);
        linenb = 1;
        while inputline~=-1
            colnb = 1;
            while ~isempty(deblank(inputline))
                [tmp inputline] = strtok(inputline);
                tmp2 = str2num( tmp );
                if isempty( tmp2 ), array{linenb, colnb} = tmp;
                else                array{linenb, colnb} = tmp2;
                end;
                colnb = colnb+1;
            end;
            inputline = fgetl(fid);
            linenb = linenb +1;
        end;        
                
		fclose(fid);

    else % variable in the global workspace
         % --------------------------
         array = evalin('base', varname);
    end;     
return;

