% pop_importepoch() - Export epoch information to the event structure array
%                  of an EEG dataset. If the dataset is the only input, 
%                  a window pops up to ask for the relevant parameter values.
%
% Usage:
%   >> EEGOUT = pop_importepoch( EEG, filename, fieldlist, latencyfieldlist, ...
%                                  timeunit, headerlines, cleanevents);
%
% Inputs:
%   EEG              - input dataset
%   filename         - name of a file or an array with epoch information
%                      organised in columns. It can be either an array or
%                      a cell array of values.
%   fieldlist        - cell array of name of each column in the file.
%   latencyfieldlist - cell array of fields (defined previously) that 
%                      contain the latency of an event. These fields are 
%                      transfered into the event structure of the input dataset.
%   timeunit         - optional latency units in second. Default is that
%                      latencies are expressed in time points unit.
%   headerlines      - number of line in the file to skip. Default is 0.
%   cleanevents      - [0|1], 1=clean the event array. Default is 1.
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

function [EEG, com] = pop_importepoch( EEG, filename, fieldlist, latencyfieldlist, timeunit, headerlines, cleanevent);
com ='';
if nargin < 1
    help pop_importepoch
    return;
end;
if nargin < 2
    geometry    = { [ 1 1 1] [1.2 0.8] [0.8 1 1.2] [1.5 0.5 1.5] [1.5 0.5 1.5] [1.5 0.2 1.8]};
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    helpstr = ['It is not necessary to define a latency field for epoch information.' 10 ...
			   'Actually all fields that contain latencies will be imported as different event types.' 10 ...
			   'For instance, if entering ''rt'', events of type ''rt'' will be created and these events will have a latency'];
	uilist = { ...
         { 'Style', 'text', 'string', 'Epoch file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''globfile'';' commandload ] }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         ...
         { 'Style', 'text', 'string', 'All input field (column) names (ex: type latency)', 'fontweight', 'bold' }, { 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', '           Out of which', 'horizontalalignment', 'right', 'fontweight', 'bold', 'tooltipstring', helpstr },  ...
		{ 'Style', 'edit', 'string', '' }, ...
         { 'Style', 'text', 'string', 'contain latencies (ex: rt)', 'fontweight', 'bold', 'tooltipstring', helpstr }, ...
         ...
         { 'Style', 'text', 'string', 'File header lines', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '0' }, { },...        
         { 'Style', 'text', 'string', 'Latency time unit (sec) Ex: 1E-3 = ms', 'horizontalalignment', 'left' }, { 'Style', 'edit', 'string', '1' }, { } ...         
         { 'Style', 'text', 'string', 'Clean previous info (set=yes)', 'horizontalalignment', 'left' }, { 'Style', 'checkbox', 'value', isempty(EEG.event) }, { } };         
    result = inputgui( geometry, uilist, 'pophelp(''pop_importepoch'');', 'Import epoch info (data epochs only) -- pop_importepoch()');
    if length(result) == 0, return; end;

    filename    = result{1};
    fieldlist   = parsetxt( result{2} );
    latencyfieldlist = parsetxt( result{3} );
    headerlines = eval( result{4} );
    timeunit    = eval( result{5} );
    cleanevent    = result{6};
end;

if exist('timeunit') ~= 1
    timeunit = 1/EEG.srate;
end;    
if exist('headerlines') ~= 1
    headerlines = 1;
end;    
if exist('cleanevent') ~= 1
    cleanevent = 1;
end;
insertepoch = 1;
if cleanevent == 0
    insertepoch == 0;
end;
if ~isfield( EEG.event, 'epoch')
    insertepoch = 1;
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
    values = load_file_or_array( filename, headerlines );
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
if ~iscell(latencyfieldlist)
    latencyfieldlist = { latencyfieldlist };
end;
otherfieldlist = setdiff( fieldlist, latencyfieldlist);
if size(values,1) ~= EEG.trials
    error('Pop_importepoch error: the number of row of the input file/array does not match the number of trials');
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

if isempty( EEG.epoch)
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
        cleanevent = 1;
        disp('Pop_importepoch: cannot add events to a non-epoch event structure, erasing old epoch structure');
    end;
end;
if cleanevent
    EEG.event = [];
end;
           
% add latency fields
% ------------------
for index = 1:length(latencyfieldlist)
    if isempty(strmatch( latencyfieldlist{index}, epochfield )) 
         error(['Pop_importepoch: field ''' latencyfieldlist{index} ''' not found']);
    end;
    for trial = 1:EEG.trials
       if insertepoch, 
            EEG.event(end+1).epoch = trial; 
            eval( ['EEG.event(end).type = '''  latencyfieldlist{index} ''';'] );
       else eval( ['EEG.event(end+1).type = '''  latencyfieldlist{index} ''';'] );
       end;
       eval( ['EEG.event(end).latency = (EEG.epoch(trial).' latencyfieldlist{index} '*timeunit-EEG.xmin)*EEG.srate+1+(trial-1)*EEG.pnts;' ] );
    end;
end;

% add non latency fields
% ----------------------
if ~isfield(EEG.event, 'epoch')
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
EEG = eeg_checkset(EEG, 'eventconsistency');

% generate the output command
% ---------------------------
if isempty(filename) & nargout == 2
    disp('Pop_importepoch: cannot generate command string'); return;
else 
	com = sprintf('%s = pop_importepoch( %s, ''%s'', { ', inputname(1), inputname(1), filename);
	for i=1:length(fieldlist)
		com = sprintf('%s ''%s'',', com, fieldlist{i} );
	end;    
	com = [ com(1:end-1) '} , { ' ];
	for i=1:length(latencyfieldlist)
		com = sprintf('%s ''%s'',', com, latencyfieldlist{i} );
	end;    
	com = [ com(1:end-1) '}'];
	com = sprintf('%s, %d, %d, %d);', com, timeunit, headerlines, cleanevent);
	return;
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

