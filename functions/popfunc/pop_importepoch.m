% pop_importepoch() - Export epoch and/or epoch event information to the event 
%                     structure array of an EEG dataset. If the dataset is 
%                     the only input, a window pops up to ask for the relevant 
%                     parameter values.
% Usage:
%   >> EEGOUT = pop_importepoch( EEG, filename, fieldlist, 'key', 'val', ...);
%
% Inputs:
%   EEG              - Input EEG dataset
%   filename         - Name of an ascii file with epoch and/or epoch event information 
%                      organised in columns. ELSE, name of a Matlab variable with the 
%                      same information (either a Matlab array or cell array). 
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
%   'timeunit'       - [float] Optional unit for latencies relative to seconds. 
%                      Ex: sec -> 1, msec -> 1e-3. Default: Assume latencies 
%                      are in time points (relative to the time-zero time point 
%                      in the epoch). 
%   'headerlines'    - [int] Number of header lines in the input file to ignore. 
%                      {Default 0}.
%   'clearevents'    - ['on'|'off'], 'on'-> clear the current event array. 
%                      {Default 'on'}
%
% Output:
%   EEGOUT - EEG dataset with modified event structure
%
% Authors: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 11 March 2002
%
% See also: eeglab()
 
%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.20  2002/11/15 00:51:21  arno
% debugging TLE latencies, add more feedback to the user
%
% Revision 1.19  2002/10/29 22:57:03  scott
% text
%
% Revision 1.18  2002/10/29 22:55:49  scott
% text
%
% Revision 1.17  2002/10/29 22:54:46  scott
% text .
%
% Revision 1.16  2002/10/29 17:19:18  arno
% rework box size
% s
%
% Revision 1.15  2002/10/28 23:43:46  scott
% time-lock event -> TLE ; edited help message and popup strings -sm
%
% Revision 1.14  2002/10/28 20:35:12  arno
% new version, different syntax, different text (added optional type)
%
% Revision 1.13  2002/08/22 00:01:40  arno
% adding error message
%
% Revision 1.12  2002/08/06 21:55:36  arno
% spelling
%
% Revision 1.11  2002/05/02 19:31:29  arno
% debugging strmatch (exact)
%
% Revision 1.10  2002/05/02 19:15:39  arno
% typo
%
% Revision 1.9  2002/05/02 19:14:18  arno
% updating command return
%
% Revision 1.8  2002/05/02 19:12:21  arno
% editing message
%
% Revision 1.7  2002/04/26 21:36:07  arno
% correcting bug for the comments
%
% Revision 1.6  2002/04/22 21:42:58  arno
% corrected returned command
%
% Revision 1.5  2002/04/18 18:25:02  arno
% typo can not2
%
% Revision 1.4  2002/04/18 18:24:41  arno
% typo can not
%
% Revision 1.3  2002/04/11 22:42:18  arno
% debuging empty latency fields input array
%
% Revision 1.2  2002/04/11 19:37:51  arno
% additional warning for file and array with the same name
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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
    geometry    = { [ 1 1 1.86] [1] [1 0.9] [2.5 1 0.6] [2.5 1 0.6] [1] [1.5 0.5 1] [1.5 0.5 1] [1.5 0.17 1.36]};
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
	uilist = { ...
         { 'Style', 'text', 'string', 'Epoch file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''globfile'';' commandload ] }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { }...
         { 'Style', 'text', 'string', 'File input field (col.) names (not type or latency)', 'fontweight', 'bold' }, { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', '           Field name containing time-locking event type(s)', 'horizontalalignment', 'right', ...
                                      'fontweight', 'bold', 'tooltipstring', helpstrtype },  ...
  		 { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', 'NOTE', 'tooltipstring', helpstrtype }, ...
         { 'Style', 'text', 'string', '           Field name(s) containing latencies', 'horizontalalignment', 'right', ...
           'fontweight', 'bold', 'tooltipstring', helpstrlat },  ...
  		 { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', '(Ex: RT)', 'tooltipstring', helpstrlat }, ...
         { } ...
         { 'Style', 'text', 'string', 'Latency time unit rel. to seconds. Ex: ms -> 1E-3', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '1' }, { } ...         
         { 'Style', 'text', 'string', 'Number of file header lines to ignore', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '0' }, { },...        
         { 'Style', 'text', 'string', 'Remove current epoch and event info (set = yes)', 'horizontalalignment', 'left' }, { 'Style', 'checkbox', 'value', isempty(EEG.event) }, { } };         
    result = inputgui( geometry, uilist, 'pophelp(''pop_importepoch'');', 'Import epoch info (data epochs only) -- pop_importepoch()');
    if length(result) == 0, return; end;

    filename    = result{1};
    fieldlist   = parsetxt( result{2} );
    options = {};
    if ~isempty( result{3}), options = { options{:} 'typefield' result{3} }; end; 
    if ~isempty( result{4}), options = { options{:} 'latencyfields' parsetxt( result{4} ) }; end; 
    if ~isempty( result{5}), options = { options{:} 'timeunit' eval(result{5}) }; end; 
    if ~isempty( result{6}), options = { options{:} 'headerlines' eval(result{6}) }; end; 
    if ~result{7}, options = { options{:} 'clearevents' 'off'}; end; 
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

g = finputcheck( options, { 'typefield'     'string'   []       ''; ...
                            'latencyfields' 'cell'     []       {}; ...
                            'timeunit'      'real'     [0 Inf]  1/EEG.srate; ...
                            'headerlines'   'integer'  [0 Inf]  0; ...
                            'clearevents'   'string'   {'on' 'off'}  'on'}, 'pop_importepoch');
if isstr(g), error(g); end;

%g.insertepoch = 1;
%if strcmpi(g.cleanevent, 'on')
%    g.insertepoch == 0;
%end;
%if ~isfield( EEG.event, 'epoch')
%    g.insertepoch = 1;
%end;    

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
	error('There must be as many field names as there are columsn in the file/array');
end;
if ~iscell(fieldlist)
    otherfieldlist = { fieldlist };
    fieldlist = { fieldlist };
end;
if ~iscell(g.latencyfields)
    g.latencyfields = { g.latencyfields };
end;
otherfieldlist = setdiff( fieldlist, g.latencyfields);
otherfieldlist = setdiff( otherfieldlist, g.typefield);
if size(values,1) ~= EEG.trials
    error('Pop_importepoch() error: the number of rows in the input file/array does not match the number of trials');
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
    fprintf('Pop_importepoch: deleting old events if any\n');
    EEG.event = [];
else 
    fprintf('Pop_importepoch: appending new events to the existing event array\n');
end;
           
% add time locking event fields
% -----------------------------
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
    eval( ['EEG.event(end).latency = -EEG.xmin*EEG.srate+1+(trial-1)*EEG.pnts;' ] );
end;

% add latency fields
% ------------------
for index = 1:length(g.latencyfields)
    if isempty(strmatch( g.latencyfields{index}, epochfield )) 
         error(['Pop_importepoch: latency field ''' g.latencyfields{index} ''' not found']);
    end;
    for trials = 1:EEG.trials
        EEG.event(end+1).epoch = trials; 
        eval( ['EEG.event(end).type = '''  g.latencyfields{index} ''';'] );
        eval( ['EEG.event(end).latency = (EEG.epoch(trials).' g.latencyfields{index} '*g.timeunit-EEG.xmin)*EEG.srate+1+(trials-1)*EEG.pnts;' ] );
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
if isempty( EEG.eventdescription )
	allfields = fieldnames(EEG.event);
    EEG.eventdescription{strmatch('epoch', allfields, 'exact')} = 'Epoch number';
	if ~isempty(strmatch('type', allfields)), EEG.eventdescription{strmatch('type', allfields)} = 'Event type'; end;
	if ~isempty(strmatch('latency', allfields)), EEG.eventdescription{strmatch('latency', allfields)} = 'Event latency'; end;
end;

% checking and updating events
% ----------------------------
EEG = pop_editeventvals( EEG, 'sort', { 'epoch', 0 } ); % resort fields
EEG = eeg_checkset(EEG, 'eventconsistency');

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

