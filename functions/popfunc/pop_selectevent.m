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
%                      Ex: 'latency','400 < 700' Include all events with
%                           latnecy in the range [400,700]
%   'omitlatency' - [latency_range] latency range of the events to exclude
%   'events'      - [events_range], indices of the events to include
%   'omitevents'  - [events_range], indices of the events to exclude
%   'USER_VAR'    - [VAR_range], 'USER_VAR' is any user-defined field in
%                   the event structure. Includes events with values of
%                   field 'USER_VAR' in the specified range. Use [vector]
%                   format for integers, 'min<max' format for real numbers.
%   'omitUSER_VAR' - [VAR_range], 'USER_VAR' range of events to exclude
%   'deleteothers' - ['yes'|'no'] 'yes' = Delete ALL epochs that do not include
%                   the specified events. {NOTE Default = 'yes'}.
%                   This option is relevant only for epoched datasets derived
%                   from continuous datasets.
% Outputs:
%   EEGOUT - EEG dataset with the selected events only
%   event_indices - indexes of the selected events
%
%   Ex:  [EEGTARGETS,target_indices] = getevent(EEG,'type',[1 6 11 16 21]);
%
%        % Returns ONLY THOSE epochs containing the 5 specified
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
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

%02/01/2002 added inputgui and finalize function - ad
%02/06/2002 work on the header and event format - sm & ad
%02/09/2002 modify function according to new structure - ad
%02/12/2002 add getepoch compatibility - ad
%02/19/2002 add event indices & correct event selection - ad

function [EEG, event_indices, com] = pop_selectevent(EEG, varargin);

com ='';
if nargin < 1
   help pop_selectevent;
   return;
end;	

% note that this function is also used for epochs
% -----------------------------------------------
if isfield(EEG, 'tmpevent'), txtgui = 'Epoch'; helpfunction = 'pophelp(''getepoch'');';
                             txtdel = 'Remove non selected epochs';
else                         txtgui = 'Event'; helpfunction = 'pophelp(''getevent'');';
                             txtdel = 'Remove epochs not referenced by any selected event';
end;

I = [];
if isempty(EEG.event)
    disp('Getevent: can not deal with empty event structure');
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
    geometry = { [0.8 1 2.3 0.6 ] [0.8 1.1 1.8 1 ] [0.65 0.85 1.5 0.25 0.25 0.1] };
    uilist = { ...
         { 'Style', 'text', 'string', 'Selection', 'horizontalalignment', 'center', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Field Description', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Range  (int vector or real range "min < max")', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'If set,', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', '  Field', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', '(See Edit event info)'  }, ...
         { 'Style', 'text', 'string', '         Ex: 2:4,5,9:133 OR 4.5 < 13'  }, ...
         { 'Style', 'text', 'string', '    Remove range', 'fontweight', 'bold'  }, ...
         ...
         { 'Style', 'text', 'string', [ txtgui '(s)' ]}, ...
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
             if length(tmptext) > 15,    stringtext = [ tmptext(1:15) '...' ]; 
             else                        stringtext = tmptext; 
             end;
        else stringtext = 'no-description'; tmptext = 'no-description';
        end;
        descrip = { 'string', stringtext, 'callback', ['questdlg([''' choptext( tmptext ) '''],''Description of field ' allfields{index} ''', ''OK'', ''OK'');' ] };

        % create the gui for this field
        % -----------------------------
        geometry = { geometry{:} [0.65 0.85 1.5 0.25 0.25 0.1] };
        uilist   = { uilist{:}, ...
         { 'Style', 'text', 'string', [allfields{index} '(s)'] }, ...
         { 'Style', 'pushbutton', descrip{:}, 'horizontalalignment', 'left' }, ...
         { 'Style', 'edit', 'string', '' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ' },{ } }; 
    end;

    geometry = { geometry{:} [1]  [2 1] };
    uilist   = { uilist{:} ...
        { }, ...
        { 'Style', 'checkbox', 'string', txtdel, 'tag', '10', 'fontweight', 'bold', 'value', 1  } { } };
    results = inputgui( geometry, uilist, helpfunction, 'Select events -- pop_selectevent()');
    if length(results) == 0, return; end;
   
    % decode inputs
    % -------------
    args = {};
    if ~results{2}, args = { args{:},     'event', eval( [ '[' results{1} ']' ]) };
    else            args = { args{:}, 'omitevent', eval( [ '[' results{1} ']' ]) }; 
    end;
    for index = 1:length(allfields) 
        tmpres = results{2*index+1};
        if isempty(findstr(tmpres, '<')), 
            try, tmpres = eval( [ '[' tmpres ']' ] ); 
            catch, tmpres = parsetxt( tmpres ); end;
        end;
        if ~results{2*index+2}, args = { args{:}, allfields{index}, tmpres };
        else                    args = { args{:}, [ 'omit' allfields{index}], tmpres }; 
        end;
    end;    
    if results{end},  args = { args{:}, 'deleteothers', 'yes'}; end;
else % no interactive inputs
    args = varargin;
    for i=1:length(varargin)
        if iscell(args{i}), args{i} = { args{i} }; end; % double nested 
    end;    
end;

% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, error('Wrong syntax in function arguments'); end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try, if isempty(g.event) g.event = [1:length(EEG.event)]; end; catch, g.event = [1:length(EEG.event)]; end;
try, g.omitevent; 	 	 catch, g.omitevent = []; end;
try, g.deleteothers; 	 catch, g.deleteothers = 'no'; end;
for index = 1:length(allfields) 
   try, eval(['g.' allfields{index} ';']); catch, eval(['g.' allfields{index} '=[];' ]); end; 
   try, eval(['g.omit' allfields{index} ';']); catch, eval(['g.omit' allfields{index} '=[];' ]); end; 
end;
g

% select the events to keep
% -------------------------
Ievent = g.event;
Ieventrem = g.omitevent;

for index = 1:length(allfields)

    % scan each field of EEG.event
    % ----------------------------
    tmpvar = eval(['g.' allfields{index} ]);
    if ~isempty( tmpvar )
        if  iscell( tmpvar )
            eval( [ 'tmpvarvalue = {EEG.event(:).' allfields{index} '};'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                Ieventtmp = unique( [ Ieventtmp strmatch( tmpvar{index2}, tmpvarvalue, 'exact') ]);
            end;
            Ievent = intersect( Ievent, Ieventtmp );
        elseif isstr( tmpvar )
            eval( [ 'tmpvarvalue = cell2mat( {EEG.event(:).' allfields{index} '});'] );
            min = eval(tmpvar(1:findstr(tmpvar, '<')-1));
            max = eval(tmpvar(findstr(tmpvar, '<')+1:end));
            Ieventlow  = find( tmpvarvalue >= min);
            Ieventhigh = find( tmpvarvalue <= max);
            Ievent = intersect( Ievent, intersect( Ieventlow, Ieventhigh ) )
        else
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
    if ~isempty( tmpvar )
        if  iscell( tmpvar )
            eval( [ 'tmpvarvalue = {EEG.event(:).' allfields{index} '};'] );
            Ieventtmp = [];
            for index2 = 1:length( tmpvar )
                Ieventtmp = unique( [ Ieventtmp strmatch( tmpvar{index2}, tmpvarvalue, 'exact') ]);
            end;
            Ieventrem = union( Ieventrem, Ieventtmp );
         elseif isstr( tmpvar )
            eval( [ 'tmpvarvalue = cell2mat( {EEG.event(:).' allfields{index} '});'] );
            min = eval(tmpvar(1:findstr(tmpvar, '<')-1));
            max = eval(tmpvar(findstr(tmpvar, '<')+1:end));
            Ieventlow  = find( tmpvarvalue >= min);
            Ieventhigh = find( tmpvarvalue <= max);
            Ieventrem = union( Ieventrem, intersect( Ieventlow, Ieventhigh ) );
        else
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

% Events: delete epochs
% ---------------------
switch lower(g.deleteothers)
 case 'yes', % ask for confirmation
			 % --------------------
			 Iepoch = ones(1, EEG.trials);
			 for index = 1:length(Ievent)
				 Iepoch(EEG.event(Ievent(index)).epoch) = 0;
			 end;
			 Iepoch = find(Iepoch == 0);
			 if nargin < 2 
				 ButtonName=questdlg([ 'Warning: ' num2str(length(Ievent)) ' events selected' 10 ...
					'Delete '  num2str(EEG.trials-length(Iepoch)) ' un-referenced epochs ?' ], ...
									 'Confirmation', ...
									 'No', 'Cancel', 'Yes','Yes');
			 else ButtonName = 'Yes'; end;
			 
			 switch lower(ButtonName),
			  case 'cancel', return; 
			  case 'yes',
			   length(Ievent)
			   EEG.event = EEG.event(Ievent);
			   EEG = pop_select(EEG, 'trial', Iepoch);
			  case 'no',
			   EEG.event = EEG.event(Ievent);
			 end % switch
			 
 otherwise, EEG.event = EEG.event(Ievent);
end;

event_indices = Ievent;

% generate the output command
% ---------------------------
com = sprintf('EEG = getevent( %s', inputname(1));
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if isstr( args{i+1} )   com = sprintf('%s, ''%s'', ''%s''', com, args{i}, args{i+1} );
        else                    com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
        end;
    end;    
end;
com = [com ');'];

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
