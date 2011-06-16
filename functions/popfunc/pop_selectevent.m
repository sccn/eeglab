% pop_selectevent() - Find events in an EEG dataset. If the dataset
%                     is the only input, a window pops up to
%                     ask for the relevant parameter values.
%
% Usage: >> [EEGOUT,event_indices] = pop_selectevent( EEG, 'key1', value1, ...
%                                                   'key2', value2, ... );
% Input:
%   EEG  - EEG dataset
%
% Optional inputs:
%   'latency'     - [latency_range] latency range of events to include
%                      Ex: 'latency','400 <= 700' Include all events with
%                           latnecy in the range [400,700]
%   'omitlatency' - [latency_range] latency range of events to exclude
%   'type'        - [type_range] event type(s) to include
%                      Ex: 'type',3  or [2 3 5] or 1:10
%   'omittype'    - [type_range], type(s) of events to exclude
%   'event'       - [event_range], indices of events to include
%   'omitevent'   - [event_range], indices of events to exclude
%   'USER_VAR'    - [VAR_range], 'USER_VAR' is any user-defined field in
%                   the event structure. Includes events with values of
%                   field 'USER_VAR' in the specified range. Use [vector]
%                   format for integers, 'min<max' format for real numbers.
%   'omitUSER_VAR' - [VAR_range], 'USER_VAR' range of events to exclude
%   'select'       - ['normal'|'inverse'] invert the selection of events. 
%                    {Default is 'normal'}
%   'deleteepochs' - ['on'|'off'] 'on' = Delete ALL epochs that do not include
%                   any of the specified events {Default = 'on'}.
%                   This option is relevant only for epoched datasets derived
%                   from continuous datasets.
%   'invertepochs' - ['on'|'off'] 'on' = Invert epoch selection. {Default = 'off'}.
%   'deleteevents' - ['on'|'off'] 'on' = Delete ALL events except
%                   the selected events. {Default = 'off'}.
%   'renametype'   - [string] rename the type of selected events with the
%                  string given as parameter. {Default is [], do not rename
%                  field}.
%   'oldtypefield' - [string] in conjunction with the previous parameter, 
%                  create a new field (whose 'name' is provided as parameter)
%                  to store the (old) type of the event whose type has been 
%                  renamed. {Default is [], do not create field}.
%
% Outputs:
%   EEGOUT - EEG dataset with the selected events only
%   event_indices - indexes of the selected events
%
%   Ex:  [EEGTARGETS,target_indices] = pop_selectevent(EEG,'type',[1 6 11 16 21]);
%
%        % Returns ONLY THOSE epochs containing any of the 5 specified
%          types of target events.
%
% Note: By default, if several optional inputs are given, the function
%       performs their conjunction (&).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002-
%
% See also: eeg_eventformat(), pop_importevent()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002, arno@salk.edu
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

%02/01/2002 added inputgui and finalize function - ad
%02/06/2002 work on the header and event format - sm & ad
%02/09/2002 modify function according to new structure - ad
%02/12/2002 add getepoch compatibility - ad
%02/19/2002 add event indices & correct event selection - ad

function [EEG, Ievent, com] = pop_selectevent(EEG, varargin);

if nargin < 1
   help pop_selectevent;
   return;
end;
Ievent = [];
com ='';
event_indices = [];

% note that this function is also used for epochs
% -----------------------------------------------
I = [];
if isempty(EEG.event)
    disp('Getevent: cannot deal with empty event structure');
    return;
end;   

% remove the event field if present
% ---------------------------------
allfields = fieldnames(EEG.event);
indexmatch = strmatch('urevent', allfields);
if ~isempty(indexmatch)
    allfields = { allfields{1:indexmatch-1} allfields{indexmatch+1:end} };
end;

if nargin<2
    geometry = { [0.6 2.1 1.2 0.8 ] };
    uilist = { ...
         { 'Style', 'text', 'string', 'Field', 'horizontalalignment', 'center', 'fontweight', 'bold'  }, ...
         {} ...
         { 'Style', 'text', 'string', 'Selection', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Set=NOT THESE', 'fontweight', 'bold'  } };
					   
    % add all fields to graphic interface
    % -----------------------------------
    ind1 = strmatch('type', allfields, 'exact');
    ind2 = strmatch('latency' , allfields, 'exact');
    ind3 = strmatch('duration', allfields, 'exact');
    neworder = [ ind2 ind3 ind1 setdiff(1:length(allfields), [ind1 ind2 ind3]) ];
    allfields = { allfields{neworder} }; 
    
    for index = 1:length(allfields)
        % format the description to fit a help box
        % ----------------------------------------
        if index <= length( EEG.eventdescription )
             tmptext = EEG.eventdescription{ neworder(index) };
			 if ~isempty(tmptext)
				 if size(tmptext,1) > 15,    stringtext = [ tmptext(1,1:15) '...' ]; 
				 else                        stringtext = tmptext(1,:); 
				 end;
			 else stringtext = 'no description'; tmptext = 'no description (use menu Edit > Event Field)';
			 end;
        else stringtext = 'no description'; tmptext = 'no description (use menu Edit > Event Field)';
        end;

		descrip = { 'string', stringtext, 'callback', ['questdlg2(' vararg2str(tmptext) ...
					',''Description of field ''''' allfields{index} ''''''', ''OK'', ''OK'');' ] };

        % create the gui for this field
        % -----------------------------
        textfield = allfields{index};
        if strcmp(textfield, 'latency') | strcmp(textfield, 'duration')
            if EEG.trials > 1, textfield = [ textfield ' (ms)' ];
            else textfield = [ textfield ' (s)' ];
            end;
            middletxt  = { { 'Style', 'text', 'string', 'min' } { 'Style', 'edit', 'string', '' 'tag' [ 'min' allfields{index} ] } ...
                           { 'Style', 'text', 'string', 'max' } { 'Style', 'edit', 'string', '' 'tag' [ 'max' allfields{index} ] } };
            middlegeom = [ 0.3 0.35 0.3 0.35 ];
        elseif strcmp(textfield, 'type')
            commandtype = [ 'if ~isfield(EEG.event, ''type'')' ...
                           '   errordlg2(''No type field'');' ...
                           'else' ...
                           '   tmpevent = EEG.event;' ...
                           '   if isnumeric(EEG.event(1).type),' ...
                           '        [tmps,tmpstr] = pop_chansel(unique([ tmpevent.type ]));' ...
                           '   else,' ...
                           '        [tmps,tmpstr] = pop_chansel(unique({ tmpevent.type }));' ...
                           '   end;' ...
                           '   if ~isempty(tmps)' ...
                           '       set(findobj(''parent'', gcbf, ''tag'', ''type''), ''string'', tmpstr);' ...
                           '   end;' ...
                           'end;' ...
                           'clear tmps tmpv tmpevent tmpstr tmpfieldnames;' ];
            middletxt  = { { 'Style', 'edit', 'string', '' 'tag' 'type' } { 'Style', 'pushbutton', 'string', '...' 'callback' commandtype } };
            middlegeom = [ 0.95 0.35 ];
        else
            middletxt  = { { 'Style', 'edit', 'string', '' 'tag' textfield } };
            middlegeom = 1.3;
        end;
        geometry = { geometry{:} [0.55 0.65 middlegeom 0.1 0.22 0.1] };
        uilist   = { uilist{:}, ...
         { 'Style', 'text', 'string', textfield }, ...
         { 'Style', 'pushbutton', descrip{:}, 'horizontalalignment', 'left' }, ...
         middletxt{:}, ...
         { }, { 'Style', 'checkbox', 'string', '    ' 'tag' [ 'not' allfields{index} ] },{ } };
    end;
    
    % event indices
    % -------------
    uilist = { uilist{:} ...
            { 'Style', 'text', 'string', 'Event indices' }, ...
            { }, ...
            { 'Style', 'edit', 'string', '' 'tag' 'indices' }, ...
            { }, { 'Style', 'checkbox', 'string', '    ' 'tag' 'notindices' },{ } };
    geometry = { geometry{:} [0.55 0.65 1.3 0.1 0.22 0.1] };
    
    % rename/keep events
    % ------------------
    geometry = { geometry{:} [1] [1] [.1 2 .3 .2] [.1 1.5 0.5 0.5]  [.1 1 0.5 1] [.1 1 0.5 1] };
    uilist = { uilist{:} { } ...
                { 'Style', 'text', 'string','Event selection', 'fontweight', 'bold' } ...
                {} { 'Style', 'checkbox', 'string','Select all events NOT selected above (Set this button and "all BUT" buttons (above) for logical OR)' 'tag' 'invertevent' } { } { } ...
                {} { 'Style', 'checkbox', 'string','Keep only selected events and remove all other events', ...
                'value', fastif(EEG.trials>1, 0, 1) 'tag' 'rmevents' } { } { } ...
                {} { 'Style', 'text', 'string', 'Rename selected event type(s) as type:' } ...
                   { 'Style', 'edit', 'string', ''  'tag' 'rename' } { } ...
                {} { 'Style', 'text', 'string', 'Retain old event type name(s) in (new) field named:' } ...
                   { 'Style', 'edit', 'string', '' 'tag' 'retainfield'  } { }  };
    
    % epoch selections
    % ----------------
    if EEG.trials > 1
        geometry = { geometry{:} [1] [0.1 2 0.5 0.5] [0.1 2 0.5 0.5]};
        uilist   = { uilist{:} ...
                     { 'Style', 'text', 'string','Epoch selection', 'fontweight', 'bold' } ...
                     { } { 'Style', 'checkbox', 'string','Remove epochs not referenced by any selected event', ...
                       'value', 1  'tag' 'rmepochs' } { } { } ...
                     { } { 'Style', 'checkbox', 'string','Invert epoch selection', ...
                       'value', 0 'tag' 'invertepoch' } { } { } };
    end;
    
	[results tmp2 tmp3 res] = inputgui( geometry, uilist, 'pophelp(''pop_selectevent'')', 'Select events -- pop_selectevent()');
    if length(results) == 0, return; end;
   
    % decode inputs
    % -------------
    args = {};
    if ~res.notindices, args = { args{:},     'event', eval( [ '[' res.indices ']' ]) };
    else                args = { args{:}, 'omitevent', eval( [ '[' res.indices ']' ]) }; 
    end;
    for index = 1:length(allfields)
        textfield = allfields{index};
        tmpflag = getfield(res, [ 'not' textfield ]);
        if strcmp(textfield, 'duration') | strcmp(textfield, 'latency') 
            tmpres = [];
            minlat = getfield(res, [ 'min' textfield ]);
            maxlat = getfield(res, [ 'max' textfield ]);
            if ~isempty(minlat) & ~isempty(maxlat)
                tmpres = [ minlat '<=' maxlat ];
            end;
        else
            tmpres  = getfield(res, textfield);
            try, tmpres2 = eval( [ '[' tmpres ']' ] );
                if ~isnumeric(tmpres2),
                    if tmpres(1) == ''''
                        tmpres = eval( [ '{' tmpres '}' ] );
                    else
                        tmpres = parsetxt( tmpres ); 
                    end;
                else
                    tmpres = tmpres2;
                end;
            catch, tmpres = parsetxt( tmpres ); end;
        end
        if ~isempty(tmpres)
            if ~tmpflag, args = { args{:}, textfield, tmpres };
            else         args = { args{:}, [ 'omit' textfield], tmpres }; 
            end;
        end;
    end;
    if res.invertevent,  args = { args{:}, 'select', 'inverse' }; end;
    if ~isempty(res.rename),       args = { args{:}, 'renametype', res.rename }; end;
    if ~isempty(res.retainfield),  args = { args{:}, 'oldtypefield', res.retainfield }; end;
    args = { args{:}, 'deleteevents', fastif(res.rmevents,     'on', 'off') };
    if EEG.trials > 1
        args = { args{:}, 'deleteepochs', fastif(res.rmepochs    , 'on', 'off') };        
        args = { args{:}, 'invertepochs', fastif(res.invertepoch , 'on', 'off') };
    end;
else % no interactive inputs
    args = varargin;
end;

% setting default for the structure
% ---------------------------------
fieldlist = { 'event'         'integer'     []                                       [1:length(EEG.event)] ;
			  'omitevent'     'integer'     []                                       [] ;
			  'deleteepochs'  'string'      { 'yes','no','on','off' }                'on' ;
			  'invertepochs'  'string'      { 'on','off' }                           'off' ;
			  'deleteevents'  'string'      { 'yes','no','on','off' }                'off';
			  'renametype'    'string'      []                                       '';
			  'oldtypefield'  'string'      []                                       '';
			  'select'        'string'      { 'normal','inverse','remove','keep' }   'normal' };
for index = 1:length(allfields) 
	fieldlist{end+1, 1} = allfields{index};
	fieldlist{end  , 2} = '';
	fieldlist{end+1, 1} = [ 'omit' allfields{index} ];
	fieldlist{end  , 2} = '';
end;
g = finputcheck( args, fieldlist, 'pop_selectevent');
if isstr(g), error(g); end;
if isempty(g.event), g.event = [1:length(EEG.event)]; end;
if strcmpi(g.select, 'remove'), g.select = 'inverse'; end;
if strcmpi(g.select, 'keep'  ), g.select = 'normal'; end;
if strcmpi(g.deleteepochs, 'yes'  ), g.deleteepochs = 'on'; end;
if strcmpi(g.deleteepochs, 'no'  ),  g.deleteepochs = 'off'; end;
if ~isempty(g.oldtypefield) & isempty(g.renametype)
    error('A name for the new type must be defined');
end;

% select the events to keep
% -------------------------
Ievent    = g.event;
Ieventrem = g.omitevent;

for index = 1:length(allfields)

    % convert the value if the field is a string field
    % ------------------------------------------------
	tmpvar = getfield(g, {1}, allfields{index});
    
	if ~isempty(tmpvar)
		if isnumeric(tmpvar)
			if isstr(getfield( EEG.event, {1}, allfields{index}))
				for tmpind = 1:length(tmpvar)
					tmpvartmp{tmpind} = num2str(tmpvar(tmpind));
				end;
				tmpvar = tmpvartmp;
			end;
		elseif isstr(tmpvar) & isempty( findstr(tmpvar, '<='))
			if isnumeric(getfield( EEG.event, {1}, allfields{index}))
				error(['numerical values must be entered for field ''' allfields{index} '''']);
			end;
		end;
	end;
		
	if isstr(tmpvar) & isempty( findstr(tmpvar, '<='))
		tmpvar = { tmpvar };
	end;

    % scan each field of EEG.event
    % ----------------------------
    if ~isempty( tmpvar )
        if  iscell( tmpvar ) % strings
            eval( [ 'tmpevent = EEG.event; tmpvarvalue = {tmpevent(:).' allfields{index} '};'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                tmpindex = strmatch( tmpvar{index2}, tmpvarvalue, 'exact');
                if isempty( tmpindex ),
                    fprintf('Warning: ''%s'' field value ''%s'' not found\n', allfields{index}, tmpvar{index2});
                end;
                Ieventtmp = unique( [ Ieventtmp; tmpindex ]);
            end;
            Ievent = intersect( Ievent, Ieventtmp );
        elseif isstr( tmpvar ) % real range
            tmpevent = EEG.event;
            eval( [ 'tmpvarvalue = [ tmpevent(:).' allfields{index} ' ];'] );
            min = eval(tmpvar(1:findstr(tmpvar, '<=')-1));
            max = eval(tmpvar(findstr(tmpvar, '<=')+2:end));
			if strcmp(allfields{index}, 'latency')
				if EEG.trials > 1
                    tmpevent = EEG.event;
					tmpvarvalue = eeg_point2lat(tmpvarvalue, {tmpevent.epoch}, EEG.srate, ...
											[EEG.xmin EEG.xmax]*1000, 1E-3);
				else
					tmpvarvalue = eeg_point2lat(tmpvarvalue, ones(1,length(EEG.event)), EEG.srate, ...
											[EEG.xmin EEG.xmax], 1);
				end;
			end;
			if strcmp(allfields{index}, 'duration')
				if EEG.trials > 1, tmpvarvalue = tmpvarvalue/EEG.srate*1000;
                else               tmpvarvalue = tmpvarvalue/EEG.srate;
                end;
            end;
			Ieventlow  = find( tmpvarvalue >= min);
			Ieventhigh = find( tmpvarvalue <= max);
			Ievent = intersect( Ievent, intersect( Ieventlow, Ieventhigh ) );
        else
			if strcmp(allfields{index}, 'latency')
				fprintf(['pop_selectevent warning: latencies are continuous values\n' ...
						 'so you may use the ''a<=b'' notation to select these values\n']);
			end;
            eval( [ 'tmpevent = EEG.event; tmpvarvalue = [ tmpevent(:).' allfields{index} ' ];'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                Ieventtmp = unique( [ Ieventtmp find(tmpvarvalue == tmpvar(index2)) ] );
            end;
			Ievent = intersect( Ievent, Ieventtmp );
        end;
     end;
        
    % scan each field of EEG.event (omit)
    % -----------------------------------
    tmpvar = eval(['g.omit' allfields{index} ]);
	if eval(['isstr(EEG.event(1).' allfields{index} ')' ]) & isnumeric(tmpvar) & ~isempty(tmpvar)
		for tmpind = 1:length(tmpvar) 
			tmpvartmp{tmpind} = num2str(tmpvar(tmpind));
		end;
		tmpvar = tmpvartmp;
	end;
	if isstr(tmpvar) & isempty( findstr(tmpvar, '<='))
		tmpvar = { tmpvar };
	end;
    if ~isempty( tmpvar )
        if  iscell( tmpvar )
            eval( [ 'tmpevent = EEG.event; tmpvarvalue = {tmpevent(:).' allfields{index} '};'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                tmpindex = strmatch( tmpvar{index2}, tmpvarvalue, 'exact');
                if isempty( tmpindex ),
                    fprintf('Warning: ''%s'' field value ''%s'' not found\n', allfields{index}, tmpvar{index2});
                end;
                Ieventtmp = unique( [ Ieventtmp; tmpindex ]);
            end;
            Ieventrem = union( Ieventrem, Ieventtmp );
         elseif isstr( tmpvar )
            tmpevent = EEG.event;
            eval( [ 'tmpvarvalue = [ tmpevent(:).' allfields{index} ' ];'] );
            min = eval(tmpvar(1:findstr(tmpvar, '<=')-1));
            max = eval(tmpvar(findstr(tmpvar, '<=')+2:end));
			if strcmp(allfields{index}, 'latency')
				if EEG.trials > 1
                    tmpevent = EEG.event;
					tmpvarvalue = eeg_point2lat(tmpvarvalue, {tmpevent.epoch}, EEG.srate, ...
											[EEG.xmin EEG.xmax]*1000, 1E-3);
				else
					tmpvarvalue = eeg_point2lat(tmpvarvalue, ones(1,length(EEG.event)), EEG.srate, ...
											[EEG.xmin EEG.xmax], 1);
				end;
			end;
			if strcmp(allfields{index}, 'duration')
				if EEG.trials > 1, tmpvarvalue = tmpvarvalue/EEG.srate*1000;
                else               tmpvarvalue = tmpvarvalue/EEG.srate;
                end;
            end;
            Ieventlow  = find( tmpvarvalue > min);
            Ieventhigh = find( tmpvarvalue < max);
            Ieventrem = union( Ieventrem, intersect( Ieventlow, Ieventhigh ) );
        else
			if strcmp(allfields{index}, 'latency')
				fprintf(['pop_selectevent warning: latencies are continuous values\n' ...
						 'so you may use the ''a<=b'' notation to select these values\n']);
			end;
            tmpevent = EEG.event;
            eval( [ 'tmpvarvalue = [ tmpevent(:).' allfields{index} ' ];'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                Ieventtmp = unique( [ Ieventtmp find( tmpvarvalue ==tmpvar(index2)) ] );
            end;
            Ieventrem = union( Ieventrem, Ieventtmp );
        end;
	end;
end;

Ievent = setdiff( Ievent, Ieventrem);
if strcmp(g.select, 'inverse')
	Ievent = setdiff( [1:length(EEG.event)], Ievent );
end;

% checking if trying to remove boundary events (in continuous data)
if isfield(EEG.event, 'type')
    if isstr(EEG.event(1).type) & EEG.trials == 1 
        Ieventrem = setdiff([1:length(EEG.event)], Ievent );
        tmpevent  = EEG.event;
        boundaryindex = strmatch('boundary', { tmpevent(Ieventrem).type });
        if ~isempty(boundaryindex)
            boundaryindex = Ieventrem(boundaryindex);
            Ievent = [ Ievent boundaryindex ];
        end;
        Ievent = sort(Ievent);
    else boundaryindex = [];
    end;
else boundaryindex = [];
end;

% rename events if necessary
% --------------------------
if ~isempty(g.renametype)
    fprintf('Pop_selectevent: renaming %d selected events (out of %d)\n', length(Ievent), length(EEG.event));
    if ~isempty(g.oldtypefield)
        for index = setdiff(Ievent, boundaryindex)
            eval([ 'EEG.event(index).' g.oldtypefield '= EEG.event(index).type;']);
            EEG.event(index).type = g.renametype;
        end;
    else
        for index = setdiff(Ievent, boundaryindex)
            EEG.event(index).type = g.renametype;
        end;
    end;
end;

% Events: delete epochs
% ---------------------
if strcmp( lower(g.deleteepochs), 'on') & EEG.trials > 1
	% ask for confirmation
	% --------------------
	Iepoch = ones(1, EEG.trials);
	for index = 1:length(Ievent)
		Iepoch(EEG.event(Ievent(index)).epoch) = 0;
	end;
    if strcmpi(g.invertepochs, 'on')
        Iepoch = ~Iepoch;
    end;
	Iepoch = find(Iepoch == 0);
	if length(Iepoch) == 0,
		error('Empty dataset: all epochs have been removed');
	end;
	if nargin < 2 
		ButtonName=questdlg2(strvcat([ 'Warning: delete ' num2str(EEG.trials-length(Iepoch)) ...
                            ' (out of ' int2str(EEG.trials) ') un-referenced epochs ?' ]), ...
							'Confirmation', ...
							 'Cancel', 'Ok','Ok');
	else ButtonName = 'ok'; end;
	
	switch lower(ButtonName),
	 case 'cancel', return; 
	 case 'ok',
	  if strcmpi(g.deleteevents, 'on')
          EEG.event = EEG.event(Ievent);
      end;
      EEG = pop_select(EEG, 'trial', Iepoch);
	end % switch
else 
    % delete events if necessary
    % --------------------------
    if strcmpi(g.deleteevents, 'on')
        EEG.event = EEG.event(Ievent);
    end;
end;

EEG = eeg_checkset(EEG, 'eventconsistency');

% generate the output command
% ---------------------------
argsout = {};
for index =1:2:length(args)
	if ~isempty(args{index+1})
		argsout = { argsout{:} args{index}  args{index+1}};
	end;
end;
com = sprintf('EEG = pop_selectevent( %s, %s);', inputname(1), vararg2str(argsout));

% chop the text so that it fits into the description window
% ---------------------------------------------------------
function  chopedtext = choptext( tmptext )
    chopedtext = '';
    while length(tmptext) > 30
          blanks = findstr( tmptext, ' ');
          [tmp I] = min( abs(blanks - 30) );
          chopedtext = [ chopedtext ''' 10 ''' tmptext(1:blanks(I)) ];
          tmptext  = tmptext(blanks(I)+1:end);
    end;    
    chopedtext = [ chopedtext ''' 10 ''' tmptext];
    chopedtext = chopedtext(7:end);
return;
