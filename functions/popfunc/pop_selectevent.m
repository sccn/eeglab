% pop_selectevent() - Find events in a EEG dataset. If the dataset
%              is the only input, a window pops up to
%              ask for the relevant parameter values.
%
% Usage: >> [EEGOUT,event_indices] = pop_selectevent( EEG, 'key1', value1, ...
%                                                   'key2', value2, ... );
% Input:
%   EEG  - EEG dataset
%
% Optional inputs:
%   'type'        - [type_range] event type(s) to include
%                      Ex: 'type',3  or [2 3 5] or 1:10
%   'omittype'    - [type_range], type(s) of events to exclude
%   'latency'     - [latency_range] latency range of the events to include
%                      Ex: 'latency','400 <= 700' Include all events with
%                           latnecy in the range [400,700]
%   'omitlatency' - [latency_range] latency range of the events to exclude
%   'event'       - [event_range], indices of the events to include
%   'omitevent'   - [event_range], indices of the events to exclude
%   'USER_VAR'    - [VAR_range], 'USER_VAR' is any user-defined field in
%                   the event structure. Includes events with values of
%                   field 'USER_VAR' in the specified range. Use [vector]
%                   format for integers, 'min<max' format for real numbers.
%   'omitUSER_VAR' - [VAR_range], 'USER_VAR' range of events to exclude
%   'select'       - ['normal'|'inverse'] invert the selection of events. Default
%                    is 'normal'.
%   'deleteepochs' - ['on'|'off'] 'on' = Delete ALL epochs that do not include
%                   the specified events. {NOTE Default = 'on'}.
%                   This option is relevant only for epoched datasets derived
%                   from continuous datasets.
%   'deleteevents' - ['on'|'off'] 'on' = Delete ALL events except
%                   the selected events. {NOTE Default = 'off'}.
%   'renametype'   - [string] rename the type of selected events with the
%                  string given as parameter. Default is [], do not rename
%                  field.
%   'oldtypefield' - [string] in conjunction with the previous parameter, 
%                  create a new field (which name is provided as parameter)
%                  to store the (old) type of the event which type have been 
%                  renamed. Default is [], do not create field.
%
% Outputs:
%   EEGOUT - EEG dataset with the selected events only
%   event_indices - indexes of the selected events
%
%   Ex:  [EEGTARGETS,target_indices] = getevent(EEG,'type',[1 6 11 16 21]);
%
%        % Returns ONLY THOSE epochs containing any of the 5 specified
%          types of target events.
%
% Note: By default, if several optional inputs are given, the function
%       performs their conjunction (&).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002
%
% See also: eeg_eventformat(), pop_importevent()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.30  2002/10/29 22:30:55  scott
% text
%
% Revision 1.29  2002/10/29 22:30:08  scott
% text
%
% Revision 1.28  2002/10/29 22:25:59  scott
% text
%
% Revision 1.27  2002/10/29 17:27:17  arno
% change location of time unit message
%
% Revision 1.26  2002/10/29 17:07:07  arno
% text editing
%
% Revision 1.25  2002/10/29 16:40:35  arno
% text
%
% Revision 1.24  2002/10/29 16:14:18  arno
% debug continuous data call
%
% Revision 1.23  2002/10/29 16:07:23  arno
% new version with renaming etc ...
%
% Revision 1.22  2002/10/27 17:11:06  julie
% help msg
%
% Revision 1.21  2002/09/05 16:33:28  arno
% remove warning
%
% Revision 1.20  2002/08/28 00:57:19  arno
% [Amodifying error message
%
% Revision 1.19  2002/08/19 21:57:55  arno
% debug for MAC
%
% Revision 1.18  2002/08/19 19:06:52  arno
% debugging
%
% Revision 1.17  2002/08/13 16:15:29  scott
% text
%
% Revision 1.16  2002/08/12 18:36:10  arno
% questdlg2
%
% Revision 1.15  2002/08/08 02:09:19  arno
% contraining boundary event to be preserved for continuous data
%
% Revision 1.14  2002/08/01 22:21:41  arno
% debugging
%
% Revision 1.13  2002/07/10 02:16:41  arno
% adding the 'select' input, the function check
%
% Revision 1.12  2002/05/03 02:49:53  arno
% updating interface
%
% Revision 1.11  2002/05/03 02:42:11  arno
% editing interface
%
% Revision 1.10  2002/05/03 02:29:17  arno
% allow to select latencies
%
% Revision 1.9  2002/04/25 17:17:15  scott
% editting msg -sm
%
% Revision 1.8  2002/04/25 02:14:43  arno
% debugging multi-line event field description
%
% Revision 1.7  2002/04/22 23:37:29  arno
% temporary modification removed
%
% Revision 1.6  2002/04/22 23:36:53  arno
% temporary modif
%
% Revision 1.5  2002/04/18 18:25:51  arno
% typo can not
%
% Revision 1.4  2002/04/10 22:39:47  arno
% removing debuging comment
%
% Revision 1.3  2002/04/10 00:37:47  arno
% debuging string event selection
%
% Revision 1.2  2002/04/09 03:02:38  arno
% debuging epoch selection
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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
if isfield(EEG, 'tmpevent') & strmatch('event', allfields)
    indexmatch = strmatch('event', allfields);
	allfields = { allfields{1:indexmatch-1} allfields{indexmatch+1:end}};
end;   
 
if nargin<2
    geometry = { [0.8 1 2.3 0.6 ] [0.8 1.1 1.8 1 ] [0.65 0.85 1.3 0.45 0.25 0.1] };
    uilist = { ...
         { 'Style', 'text', 'string', 'Selection', 'horizontalalignment', 'center', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Field Description', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Range (value list or real range "min <= max")', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'If set,', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', '  Field', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'To edit: Edit/Event info'  }, ...
         { 'Style', 'text', 'string', 'Ex: 2:4,5  OR  ''COND1''  OR  4.5 <= 13'  }, ...
         { 'Style', 'text', 'string', 'select all but these', 'fontweight', 'bold'  }, ...
         ...
         { 'Style', 'text', 'string', 'Event indices' }, ...
         { }, ...
         { 'Style', 'edit', 'string', '' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ' },{ } };
					   
    % add all fields to graphic interface
    % -----------------------------------
    for index = 1:length(allfields)
        % format the description to fit a help box
        % ----------------------------------------
        if index <= length( EEG.eventdescription )
             tmptext = EEG.eventdescription{ index };
			 if ~isempty(tmptext)
				 if size(tmptext,1) > 15,    stringtext = [ tmptext(1,1:15) '...' ]; 
				 else                        stringtext = tmptext(1,:); 
				 end;
			 else stringtext = 'no-description'; tmptext = 'no-description';
			 end;
        else stringtext = 'no-description'; tmptext = 'no-description';
        end;

		descrip = { 'string', stringtext, 'callback', ['questdlg2(' vararg2str(tmptext) ...
					',''Description of field ''''' allfields{index} ''''''', ''OK'', ''OK'');' ] };

        % create the gui for this field
        % -----------------------------
        geometry = { geometry{:} [0.65 0.85 1.3 0.45 0.25 0.1] };
        uilist   = { uilist{:}, ...
         { 'Style', 'text', 'string', [allfields{index} '(s)'] }, ...
         { 'Style', 'pushbutton', descrip{:}, 'horizontalalignment', 'left' }, ...
         { 'Style', 'edit', 'string', '' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ' },{ } };
        if strcmpi(allfields{index}, 'latency')
            if EEG.trials > 1
                uilist{end-2} = { 'Style', 'text', 'string', '(ms)     ' };
            else
                uilist{end-2} = { 'Style', 'text', 'string', '(s)' };	   
            end;
        end;
    end;

    geometry = { geometry{:} [1] [1.3 2] };
    uilist   = { uilist{:} ...
        { }, ...
        { 'Style', 'checkbox', 'string','Select all events NOT selected above',} { } ...
        };

    geometry = { geometry{:} [1] [0.9 0.5 1 0.5] [2 1] };
    uilist = { uilist{:} { } ...
                { 'Style', 'text', 'string', 'Rename selected event(s) types as type:' } ...
                { 'Style', 'edit', 'string', '' } ...
                { 'Style', 'text', 'string', 'Retain old event type name(s) in (new) field named:' } ...
                { 'Style', 'edit', 'string', '' } ...
                { 'Style', 'checkbox', 'string','Keep only selected events and remove all other events', ...
                'value', fastif(EEG.trials>1, 0, 1) } { } };

    if EEG.trials > 1
        geometry = { geometry{:} [2 1] };
        uilist   = { uilist{:} ...
                     { 'Style', 'checkbox', 'string','Remove epochs not referenced by any selected event', ...
                       'fontweight', 'bold', 'value', 1  } { } };
    end;
    
	results = inputgui( geometry, uilist, 'pophelp(''pop_selectevent'')', 'Select events -- pop_selectevent()');
    if length(results) == 0, return; end;
   
    % decode inputs
    % -------------
    args = {};
    if ~results{2}, args = { args{:},     'event', eval( [ '[' results{1} ']' ]) };
    else            args = { args{:}, 'omitevent', eval( [ '[' results{1} ']' ]) }; 
    end;
    for index = 1:length(allfields) 
        tmpres = results{2*index+1};
        if isempty(findstr(tmpres, '<=')), 
            try, tmpres = eval( [ '[' tmpres ']' ] ); 
            catch, tmpres = parsetxt( tmpres ); end;
        end;
        if ~results{2*index+2}, args = { args{:}, allfields{index}, tmpres };
        else                    args = { args{:}, [ 'omit' allfields{index}], tmpres }; 
        end;
    end;
    if EEG.trials > 1
        if results{end-4},  args = { args{:}, 'select', 'inverse' }; end;
        if ~isempty(results{end-3}),  args = { args{:}, 'renametype', results{end-3} }; end;
        if ~isempty(results{end-2}),  args = { args{:}, 'oldtypefield', results{end-2} }; end;
        args = { args{:}, 'deleteevents', fastif(results{end-1}, 'on', 'off') };
        args = { args{:}, 'deleteepochs', fastif(results{end}, 'on', 'off') };
    else
        if results{end-3},  args = { args{:}, 'select', 'inverse' }; end;
        if ~isempty(results{end-2}),  args = { args{:}, 'renametype', results{end-2} }; end;
        if ~isempty(results{end-1}),  args = { args{:}, 'oldtypefield', results{end-1} }; end;
        args = { args{:}, 'deleteevents', fastif(results{end}, 'on', 'off') };
    end;
else % no interactive inputs
    args = varargin;
    for i=1:length(varargin)
        if iscell(args{i}), args{i} = { args{i} }; end; % double nested 
    end;    
end;

% setting default for the structure
% ---------------------------------
fieldlist = { 'event'         'integer'     []                                       [1:length(EEG.event)] ;
			  'omitevent'     'integer'     []                                       [] ;
			  'deleteepochs'  'string'      { 'yes' 'no' 'on' 'off' }                'on' ;
			  'deleteevents'  'string'      { 'yes' 'no' 'on' 'off' }               'off';
			  'renametype'    'string'      []                                       '';
			  'oldtypefield'  'string'      []                                       '';
			  'select'        'string'      { 'normal' 'inverse' 'remove' 'keep' }   'normal' };
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
Ievent = g.event;
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
	
	if isstr(tmpvar) & isempty( findstr(tmpvar, '<='))
		tmpvar = { tmpvar };
	end;

    % scan each field of EEG.event
    % ----------------------------
    if ~isempty( tmpvar )
        if  iscell( tmpvar ) % strings
            eval( [ 'tmpvarvalue = {EEG.event(:).' allfields{index} '};'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                Ieventtmp = unique( [ Ieventtmp; strmatch( tmpvar{index2}, tmpvarvalue, 'exact') ]);
            end;
            Ievent = intersect( Ievent, Ieventtmp );
        elseif isstr( tmpvar ) % real range
            eval( [ 'tmpvarvalue = cell2mat( {EEG.event(:).' allfields{index} '});'] );
            min = eval(tmpvar(1:findstr(tmpvar, '<=')-1));
            max = eval(tmpvar(findstr(tmpvar, '<=')+2:end));
			if strcmp(allfields{index}, 'latency')
				if EEG.trials > 1
					tmpvarvalue = eeg_point2lat(tmpvarvalue, {EEG.event.epoch}, EEG.srate, ...
											[EEG.xmin EEG.xmax]*1000, 1E-3);
				else
					tmpvarvalue = eeg_point2lat(tmpvarvalue, ones(1,length(EEG.event)), EEG.srate, ...
											[EEG.xmin EEG.xmax], 1);
				end;
			end;
			Ieventlow  = find( tmpvarvalue > min);
			Ieventhigh = find( tmpvarvalue < max);
			Ievent = intersect( Ievent, intersect( Ieventlow, Ieventhigh ) );
        else
			if strcmp(allfields{index}, 'latency')
				fprintf(['pop_selectevent warning: latencies are continuous values\n' ...
						 'so you may use the ''a<=b'' notation to select these values\n']);
			end;
            eval( [ 'tmpvarvalue = cell2mat( {EEG.event(:).' allfields{index} '});'] );
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
            eval( [ 'tmpvarvalue = {EEG.event(:).' allfields{index} '};'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                Ieventtmp = unique( [ Ieventtmp; strmatch( tmpvar{index2}, tmpvarvalue, 'exact') ]);
            end;
            Ieventrem = union( Ieventrem, Ieventtmp );
         elseif isstr( tmpvar )
            eval( [ 'tmpvarvalue = cell2mat( {EEG.event(:).' allfields{index} '});'] );
            min = eval(tmpvar(1:findstr(tmpvar, '<=')-1));
            max = eval(tmpvar(findstr(tmpvar, '<=')+2:end));
			if strcmp(allfields{index}, 'latency')
				if EEG.trials > 1
					tmpvarvalue = eeg_point2lat(tmpvarvalue, {EEG.event.epoch}, EEG.srate, ...
											[EEG.xmin EEG.xmax]*1000, 1E-3);
				else
					tmpvarvalue = eeg_point2lat(tmpvarvalue, ones(1,length(EEG.event)), EEG.srate, ...
											[EEG.xmin EEG.xmax], 1);
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
            eval( [ 'tmpvarvalue = cell2mat( {EEG.event(:).' allfields{index} '});'] );
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
if isfield(EEG.event, 'type') & isstr(EEG.event(1).type) & EEG.trials == 1 
	Ieventrem = setdiff([1:length(EEG.event)], Ievent );
	boundaryindex = strmatch('boundary', { EEG.event(Ieventrem).type });
	if ~isempty(boundaryindex)
		Ievent = [ Ievent Ieventrem(boundaryindex)];
	end;
	Ievent = sort(Ievent);
end;

% rename events if necessary
% --------------------------
if ~isempty(g.renametype)
    fprintf('Pop_selectevent: renaming %d selected events (out of %d)\n', length(Ievent), length(EEG.event));
    if ~isempty(g.oldtypefield)
        for index = Ievent
            eval([ 'EEG.event(index).' g.oldtypefield '= EEG.event(index).type;']);
            EEG.event(index).type = g.renametype;
        end;
    else
        for index = Ievent
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
	Iepoch = find(Iepoch == 0);
	if length(Iepoch) == 0,
		error('Empty dataset: all epochs have been removed');
	end;
	if nargin < 2 
		ButtonName=questdlg2(strvcat([ 'Warning: delete ' num2str(EEG.trials-length(Iepoch)) ...
                            ' (out of ' int2str(EEG.trials) ') un-referenced epochs ?' ]), ...
							'Confirmation', ...
							 'Cancel', 'Ok','Ok');
	else ButtonName = 'Yes'; end;
	
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
