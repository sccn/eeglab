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
%   'event'       - [event_range], indices of the events to include
%   'omitevent'   - [event_range], indices of the events to exclude
%   'USER_VAR'    - [VAR_range], 'USER_VAR' is any user-defined field in
%                   the event structure. Includes events with values of
%                   field 'USER_VAR' in the specified range. Use [vector]
%                   format for integers, 'min<max' format for real numbers.
%   'omitUSER_VAR' - [VAR_range], 'USER_VAR' range of events to exclude
%   'select'       - ['keep'|'remove'] keep or remove selected events. Default
%                    is 'keep'.
%   'deleteepochs' - ['yes'|'no'] 'yes' = Delete ALL epochs that do not include
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

function [EEG, event_indices, com] = pop_selectevent(EEG, varargin);

com ='';
if nargin < 1
   help pop_selectevent;
   return;
end;	

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
    geometry = { [0.8 1 2.3 0.6 ] [0.8 1.1 1.8 1 ] [0.65 0.85 1.5 0.25 0.25 0.1] };
    uilist = { ...
         { 'Style', 'text', 'string', 'Selection', 'horizontalalignment', 'center', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Field Description', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Range (value list or real range "min < max")', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'If set,', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', '  Field', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'To edit: Edit/Event info'  }, ...
         { 'Style', 'text', 'string', 'Ex: 2:4,5  OR  ''COND1''  OR  4.5 < 13'  }, ...
         { 'Style', 'text', 'string', '          non-range', 'fontweight', 'bold'  }, ...
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

		descrip = { 'string', stringtext, 'callback', ['questdlg(' vararg2str(tmptext) ...
					',''Description of field ''''' allfields{index} ''''''', ''OK'', ''OK'');' ] };

        % create the gui for this field
        % -----------------------------
        geometry = { geometry{:} [0.65 0.85 1.5 0.25 0.25 0.1] };
        uilist   = { uilist{:}, ...
         { 'Style', 'text', 'string', [allfields{index} '(s)'] }, ...
         { 'Style', 'pushbutton', descrip{:}, 'horizontalalignment', 'left' }, ...
         { 'Style', 'edit', 'string', '' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ' },{ } }; 
    end;

    geometry = { geometry{:} [1]  [1 2] [2 1] };
    uilist   = { uilist{:} ...
        { }, ...
        { 'Style', 'checkbox', 'string','Remove selected events',} { } ...
        { 'Style', 'checkbox', 'string','Remove epochs not referenced by any selected event', ...
		  'tag', '10', 'fontweight', 'bold', 'value', 1  } { } };

	if isfield(EEG.event, 'latency')
		if EEG.trials > 1
			uilist{end} = { 'Style', 'text', 'string', 'Note: latency unit is millisecond' };
		else
			uilist{end} = { 'Style', 'text', 'string', 'Note: latency unit is second' };	   
		end;
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
        if isempty(findstr(tmpres, '<')), 
            try, tmpres = eval( [ '[' tmpres ']' ] ); 
            catch, tmpres = parsetxt( tmpres ); end;
        end;
        if ~results{2*index+2}, args = { args{:}, allfields{index}, tmpres };
        else                    args = { args{:}, [ 'omit' allfields{index}], tmpres }; 
        end;
    end;    
    if results{end-1},  args = { args{:}, 'select', 'remove' }; end;
    if results{end},    args = { args{:}, 'deleteepochs', 'yes'}; end;
else % no interactive inputs
    args = varargin;
    for i=1:length(varargin)
        if iscell(args{i}), args{i} = { args{i} }; end; % double nested 
    end;    
end;

% setting default for the structure
% ---------------------------------
fieldlist = { 'event'         'integer'     []                       [1:length(EEG.event)] ;
			  'omitevent'     'integer'     []                       [] ;
			  'deleteepochs'  'string'      { 'yes' 'no' }           'no' ;
			  'select'        'string'      { 'keep' 'remove' }      'keep' };
for index = 1:length(allfields) 
	fieldlist{end+1, 1} = allfields{index};
	fieldlist{end  , 2} = '';
	fieldlist{end+1, 1} = [ 'omit' allfields{index} ];
	fieldlist{end  , 2} = '';
end;
g = finputcheck( args, fieldlist, 'pop_selectevent');
if isstr(g), error(g); end;
if isempty(g.event), g.event = [1:length(EEG.event)]; end;
	
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
		elseif isstr(tmpvar) & isempty( findstr(tmpvar, '<'))
			if isnumeric(getfield( EEG.event, {1}, allfields{index}))
				error(['numerical values must be entered for field ''' allfields{index} '''']);
			end;
		end;
	end;
		
	if isstr(tmpvar) & isempty( findstr(tmpvar, '<'))
		tmpvar = { tmpvar };
	end;
	
	if isstr(tmpvar) & isempty( findstr(tmpvar, '<'))
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
            min = eval(tmpvar(1:findstr(tmpvar, '<')-1));
            max = eval(tmpvar(findstr(tmpvar, '<')+1:end));
			if strcmp(allfields{index}, 'latency')
				if EEG.trials > 1
					tmpvarvalue = eeg_point2lat(tmpvarvalue, {EEG.event.epoch}, EEG.srate, ...
											[EEG.xmin EEG.xmax]*1000, 1E-3);
				else
					tmpvarvalue = eeg_point2lat(tmpvarvalue, ones(1,length(EEG.event)), EEG.srate, ...
											[EEG.xmin EEG.xmax], 1);
				end;
			end;
			Ieventlow  = find( tmpvarvalue >= min);
			Ieventhigh = find( tmpvarvalue <= max);
			Ievent = intersect( Ievent, intersect( Ieventlow, Ieventhigh ) );
        else
			if strcmp(allfields{index}, 'latency')
				fprintf(['pop_selectevent warning: latencies are continuous values\n' ...
						 'so you may use the ''a<b'' notation to select these values\n']);
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
	if isstr(eval(['EEG.event(1).' allfields{index} ';' ])) & isnumeric(tmpvar) & ~isempty(tmpvar)
		for tmpind = 1:length(tmpvar) 
			tmpvartmp{tmpind} = num2str(tmpvar(tmpind));
		end;
		tmpvar = tmpvartmp;
	end;
	if isstr(tmpvar) & isempty( findstr(tmpvar, '<'))
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
            min = eval(tmpvar(1:findstr(tmpvar, '<')-1));
            max = eval(tmpvar(findstr(tmpvar, '<')+1:end));
			if strcmp(allfields{index}, 'latency')
				if EEG.trials > 1
					tmpvarvalue = eeg_point2lat(tmpvarvalue, {EEG.event.epoch}, EEG.srate, ...
											[EEG.xmin EEG.xmax]*1000, 1E-3);
				else
					tmpvarvalue = eeg_point2lat(tmpvarvalue, ones(1,length(EEG.event)), EEG.srate, ...
											[EEG.xmin EEG.xmax], 1);
				end;
			end;
            Ieventlow  = find( tmpvarvalue >= min);
            Ieventhigh = find( tmpvarvalue <= max);
            Ieventrem = union( Ieventrem, intersect( Ieventlow, Ieventhigh ) );
        else
			if strcmp(allfields{index}, 'latency')
				fprintf(['pop_selectevent warning: latencies are continuous values\n' ...
						 'so you may use the ''a<b'' notation to select these values\n']);
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
if strcmp(g.select, 'remove')
	Ievent = setdiff( [1:length(EEG.event)], Ievent );
end;

% Events: delete epochs
% ---------------------
if strcmp( lower(g.deleteepochs), 'yes') & EEG.trials > 1
	% ask for confirmation
	% --------------------
	Iepoch = ones(1, EEG.trials);
	for index = 1:length(Ievent)
		Iepoch(EEG.event(Ievent(index)).epoch) = 0;
	end;
	Iepoch = find(Iepoch == 0);
	if nargin < 2 
		ButtonName=questdlg([ 'Warning: keeping ' num2str(length(Ievent)) ' events' 10 ...
					'Delete '  num2str(EEG.trials-length(Iepoch)) ' un-referenced epochs ?' ], ...
							'Confirmation', ...
							'No', 'Cancel', 'Yes','Yes');
	else ButtonName = 'Yes'; end;
	
	switch lower(ButtonName),
	 case 'cancel', return; 
	 case 'yes',
	  EEG.event = EEG.event(Ievent);
	  EEG = pop_select(EEG, 'trial', Iepoch);
	 case 'no',
	  EEG.event = EEG.event(Ievent);
	end % switch
else
	if EEG.trials == 1
		disp('Pop_selectevent: delete trials option ignored since the data is continuous');
	end;
	EEG.event = EEG.event(Ievent);
end;
event_indices = Ievent;

% generate the output command
% ---------------------------
argsout = {};
for index =1:2:length(args)
	if ~isempty(args{index+1})
		argsout = { argsout{:} args{index}  args{index+1}};
	end;
end;
com = sprintf('EEG = pop_selectevent( %s, %s)', inputname(1), vararg2str(argsout));

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
