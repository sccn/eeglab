% pop_editeventvals() - Edit events contained in an EEG dataset structure. 
%               If the dataset is the only input, a window pops up 
%               allowing the user to insert the relevant parameter values.
%
% Usage: >> EEGOUT = pop_editeventvals( EEG, 'key1', value1, ...
%                                    'key2', value2, ... );
% Input:
%   EEG  - EEG dataset
%
% Optional inputs:
%   'sort'        - { field1 dir1 field2 dir2 } Sort events based on field1
%                   then on optional field2. Arg dir1 indicates the sort 
%                   direction (0 = increasing, 1 = decreasing).
%   'changefield' - {num field value} Insert the given value into the specified 
%                   field in event num. (Ex: {34 'latency' 320.4})
%   'changeevent' - {num value1 value2 value3 ...} Change the values of
%                   all fields in event num.
%   'add'         - {num value1 value2 value3 ...} Insert event before
%                   event num having the specified values.
%   'delete'      - vector of indices of events to delete
%
% Outputs:
%   EEGOUT        - EEG dataset with the selected events only
%
% Ex:  EEG = pop_editeventvals(EEG,'changefield', { 1 'type' 'target'});
%        % set field type of event number 1 to 'target'
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2002
%
% See also: pop_selectevent(), pop_importevent()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 15 March 2002, arno@salk.edu
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
% Revision 1.16  2002/12/06 03:43:23  arno
% debuging event sorting
%
% Revision 1.15  2002/08/12 18:31:03  arno
% questdlg2
%
% Revision 1.14  2002/06/28 02:32:55  arno
% disabling ori fields
%
% Revision 1.13  2002/06/25 13:58:09  arno
% typo
%
% Revision 1.12  2002/05/21 20:45:23  scott
% removed ; from evalin() calls -sm
%
% Revision 1.11  2002/05/03 02:35:15  arno
% allowing sorting on latency
%
% Revision 1.10  2002/05/03 01:41:57  arno
% updating call for modifying latency
%
% Revision 1.9  2002/04/25 02:14:30  arno
% adding event field description
%
% Revision 1.8  2002/04/22 23:47:57  arno
% debugging 2 variable sorting
%
% Revision 1.7  2002/04/19 20:38:55  arno
% debuging sorting for integer arrays
%
% Revision 1.6  2002/04/18 18:23:39  arno
% typo can not
%
% Revision 1.5  2002/04/18 15:34:07  scott
% editted help msg -sm
%
% Revision 1.4  2002/04/18 15:29:23  scott
% [same] -sm
%
% Revision 1.3  2002/04/18 15:26:41  scott
% added number of events to title -sm
%
% Revision 1.2  2002/04/09 20:54:55  arno
% debuging latency display for latency in continuous data
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 03-16-02 text interface editing -sm & ad 
% 03-18-02 automatic latency switching display (epoch/continuous) - ad & sm 
% 03-18-02 debug soring order - ad
% 03-18-02 put latencies in ms - ad, lf & sm
% 03-29-02 debug latencies in ms - ad & sm
% 04-02-02 debuging test - ad & sm

function [EEG, com] = pop_editeventvals(EEG, varargin);

com ='';
if nargin < 1
   help pop_editeventvals;
   return;
end;	

if isempty(EEG.event)
    disp('Getevent: cannot deal with empty event structure');
    return;
end;   

allfields = fieldnames(EEG.event);
if nargin<2
    % transfer events to global workspace
    evalin('base', [ 'eventtmp = ' inputname(1) '.event;' ]);

    % add field values
    % ----------------
    geometry = { 1 };
    tmpstr = sprintf('Edit event field values (currently %d events)',length(EEG.event));
    uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } };
    for index = 1:length(allfields) 
        geometry = { geometry{:} [1 1 1 1] };
        if strcmp( allfields{index}, 'latency')
            if EEG.trials > 1
               inputstr =  [ allfields{index} ' (ms)'];
               valuestr = num2str(eeg_point2lat( getfield(EEG.event,{1}, allfields{index}), EEG.event(1).epoch,EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3));
               strassign = [ 'eval([ ''eventtmp(valnum).' allfields{index} '= eeg_lat2point(editval,eventtmp(valnum).epoch,EEG.srate,[EEG.xmin EEG.xmax]*1000, 1E-3);'']);'];
            else
               inputstr =  [ allfields{index} ' (sec)'];
               valuestr = num2str((getfield(EEG.event,{1}, allfields{index})-1)/EEG.srate+EEG.xmin);
               strassign = [ 'eval([ ''eventtmp(valnum).' allfields{index} '= (editval- EEG.xmin)*EEG.srate+1;'' ]);'];
            end;   
		else inputstr =  allfields{index};
            valuestr = num2str(getfield(EEG.event,{1}, allfields{index}));
            strassign = [ 'eval([ ''eventtmp(valnum).' allfields{index} '= editval;'']);'];
		end;
		% callback for displaying help
		% ----------------------------
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
		callbackpushbutton = ['questdlg2(' vararg2str(tmptext) ...
					',''Description of field ''''' allfields{index} ''''''', ''OK'', ''OK'');' ];
		if strcmp(allfields{index}, 'ori_time') | strcmp(allfields{index}, 'ori_index'), enabling = 'off';
		else enabling = 'on'; end;
		uilist   = { uilist{:}, ...
					 { }, ...
					 { 'Style', 'pushbutton', 'string', inputstr, 'callback',callbackpushbutton  }, ...
					 { 'Style', 'edit', 'tag', allfields{index}, 'string', valuestr, 'enable', enabling, 'callback', ...
					   [ 'valnum   = str2num(get(findobj(''parent'', gcbf, ''tag'', ''numval''), ''string''));' ...
						 'editval = get(gcbo, ''string'');' ...
						 'if ~isempty(str2num(editval)), editval =str2num(editval);  end;' ...
						 strassign ...
						 'olduserdat = get( gcbf, ''userdata'');' ...
						 'if isempty(olduserdat), olduserdat = {}; end;' ...
						 'set( gcbf, ''userdata'', { olduserdat{:} ''changefield'' { valnum ''' allfields{index} ''' editval }});' ...
						 'clear editval valnum olduserdat;' ] } { } };
    end;

    % add buttons
    % -----------
    geometry = { geometry{:} [1] [1 0.7 0.7 1 0.7 0.7 1] [1 0.7 0.7 1 0.7 0.7 1] };
    callpart1 = [ 'valnum   = str2num(get(findobj(''parent'', gcbf, ''tag'', ''numval''), ''string''));' ];
    callpart2 = [ 'set(findobj(''parent'', gcbf, ''tag'', ''numval''), ''string'', num2str(valnum));' ];
    for index = 1:length(allfields) 
        if strcmp( allfields{index}, 'latency')
            if EEG.trials > 1
	             callpart2 = [ callpart2 'set(findobj(''parent'', gcbf, ''tag'', ''' allfields{index} ...
							   '''), ''string'', num2str(eeg_point2lat(eventtmp(valnum).' allfields{index} ...
							   ',eventtmp(valnum).epoch, EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3)));' ];
            else callpart2 = [ callpart2 'set(findobj(''parent'', gcbf, ''tag'', ''' allfields{index} ...
							           '''), ''string'', num2str((eventtmp(valnum).' allfields{index} '-1)/EEG.srate+EEG.xmin));' ]; 
			end;
		else callpart2 = [ callpart2 'set(findobj(''parent'', gcbf, ''tag'', ''' allfields{index} ...
						   '''), ''string'', num2str(eventtmp(valnum).' allfields{index} '));'  ];
      end;
    end;
    callpart2 = [ callpart2 'clear valnum;' ];
    
    listboxtext = 'No field selected';  
    for index = 1:length(allfields) 
         listboxtext = [ listboxtext '|' allfields{index} ]; 
    end;

    uilist   = { uilist{:}, ...
          { }, ...
          { },{ },{ }, {'Style', 'text', 'string', 'Event Num', 'fontweight', 'bold' }, { },{ },{ }, ...
          { 'Style', 'pushbutton', 'string', 'Delete event',  'callback', [callpart1 'eventtmp(valnum) = []; valnum = min(valnum,length(eventtmp));' ...
                'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
                'set( gcbf, ''userdata'', { olduserdat{:} ''delete'', valnum }); clear olduserdat' callpart2 ] }, ...
          { 'Style', 'pushbutton', 'string', '<<', 'callback', [callpart1 'valnum = max(valnum-10,1);' callpart2 ] }, ...
          { 'Style', 'pushbutton', 'string', '<',  'callback', [callpart1 'valnum = max(valnum-1,1);' callpart2 ] }, ...
          { 'Style', 'edit', 'string', '1', 'tag', 'numval', 'callback', [callpart1 'valnum = min(str2num(get(gcbo, ''string'')),length(eventtmp));' callpart2 ] }, ...
          { 'Style', 'pushbutton', 'string', '>',  'callback', [callpart1 'valnum = min(valnum+1,length(eventtmp));' callpart2 ] }, ...
          { 'Style', 'pushbutton', 'string', '>>', 'callback', [callpart1 'valnum = min(valnum+10,length(eventtmp));' callpart2 ] }, ...
          { 'Style', 'pushbutton', 'string', 'Insert event',  'callback', [callpart1 ...
                'eventtmp(end+3) = eventtmp(end);' ...
                'eventtmp(valnum+1:end-2) = eventtmp(valnum:end-3);' ...
                'eventtmp(valnum) = eventtmp(end-1);' ...
                'if isfield(eventtmp, ''epoch''), eventtmp(valnum).epoch = eventtmp(valnum+1).epoch; end;' ...
                'eventtmp = eventtmp(1:end-2);' ...
                'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
                'tmpcell = cell(1,1+length(fieldnames(eventtmp))); tmpcell{1} =valnum;' ...
                'set( gcbf, ''userdata'', { olduserdat{:} ''add'', tmpcell }); clear tmpcell olduserdat' callpart2 ] }, ...
          };

    % add sorting options
    % -------------------
    geometry = { geometry{:} [1] [1] [1 1 1] [1 1 1] [1 1.5 0.5] };
    uilist = {  uilist{:},...
         {}, { 'Style', 'text', 'string', 'Re-order events (for review only)', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Main sorting field:'  }, ...
         { 'Style', 'listbox', 'string', listboxtext }, ...
         { 'Style', 'checkbox', 'string', 'Click for decreasing order' } ...
         { 'Style', 'text', 'string', 'Secondary sorting field:'  }, ...
         { 'Style', 'listbox', 'string', listboxtext }, ...
         { 'Style', 'checkbox', 'string', 'Click for decreasing order' }, ...
         { 'Style', 'pushbutton', 'string', 'Re-sort', 'callback', 'set(findobj(''parent'', gcf, ''tag'', ''ok''), ''userdata'', ''retuninginputui'')' }, ...
         { 'Style', 'text', 'string', 'NB: after re-sorting call back this window' }, { }};
   
    [results userdat] = inputgui( geometry, uilist, 'pophelp(''pop_editeventvals'');', 'Edit event values -- pop_editeventvals()' );
    if length(results) == 0, return; end;

    % transfer events back from global workspace
    eventtmp = evalin('base', 'eventtmp');
    evalin('base', 'clear eventtmp');
    EEG.event = eventtmp;

    % handle sorting
    % --------------
    args = {};
    if results{end-3} ~= 1
        sortval = { allfields{ results{end-3}-1 } results{end-2} };
        if results{end-1} ~= 1
            sortval = { sortval{:} allfields{ results{end-1}-1 } results{end} };
        end;  
        args = { args{:}, 'sort', sortval }; 
    end;  
    
else % no interactive inputs
    args = varargin;
end;

% scan all the fields of g
% ------------------------
for curfield = 1:2:length(args)
    switch lower(args{curfield})
        case 'sort', 
            tmparg = args{ curfield+1 };
            if length(tmparg) < 2, dir1 = 0;
            else                   dir1 = tmparg{2}; 
            end;
            if length(tmparg) > 2
	            if length(tmparg) < 4, dir2 = 0;
	            else                   dir2 = tmparg{4}; 
	            end;
	            try, eval(['tmparray = cell2mat( { EEG.event.' tmparg{3} ' } );']);
	            catch, eval(['tmparray = { EEG.event.' tmparg{3} ' };']);
	            end;
				if strcmp( tmparg{3}, 'latency') & EEG.trials > 1
					tmparray = eeg_point2lat(tmparray, {EEG.event.epoch}, EEG.srate, [EEG.xmin EEG.xmax], 1);
				end;
	            [X I] = sort( tmparray );
	            if dir2 == 1, I = I(end:-1:1); end;
	            events = EEG.event(I);
	        else
	            events = EEG.event;
	        end;       
            try, eval(['tmparray = cell2mat( { events.' tmparg{1} ' } );']);
            catch, eval(['tmparray = { events.' tmparg{1} ' };']);
	        end;
			if strcmp( tmparg{1}, 'latency') & EEG.trials > 1
				tmparray = eeg_point2lat(tmparray, {events.epoch}, EEG.srate, [EEG.xmin EEG.xmax], 1);
			end;
	        [X I] = sort( tmparray );
	        if dir1 == 1, I = I(end:-1:1); end;
	        EEG.event = events(I);
	   case 'delete'
	        EEG.event(args{ curfield+1 })=[];
	   case 'changefield'
            tmpargs = args{ curfield+1 };
            if length( tmpargs ) < 3
                error('Pop_editeventvals: not enough arguments to change field value');
            end;
            valstr = reformat(tmpargs{3}, strcmp(tmpargs{2}, 'latency'), EEG.trials > 1, tmpargs{1} );
            eval([ 'EEG.event(' int2str(tmpargs{1}) ').'  tmpargs{2} '=' fastif(isempty(valstr), '[]', valstr) ';' ]);
	   case 'add'
            tmpargs = args{ curfield+1 };
            allfields = fieldnames(EEG.event);
            if length( tmpargs ) < length(allfields)+1
                error('Pop_editeventvals: not enough arguments to change all field values');
            end;
            num = tmpargs{1};
            EEG.event(end+1) = EEG.event(end);
            EEG.event(num+1:end) = EEG.event(num:end-1);
            for index = 1:length( allfields )
                valstr = reformat(tmpargs{index+1}, strcmp(allfields{index}, 'latency'), EEG.trials > 1, num );
                eval([ 'EEG.event(' int2str(num) ').' allfields{index} '=' fastif(isempty(valstr), '[]', valstr) ';' ]);
	        end;
	   case 'changeevent'
            tmpargs = args{ curfield+1 };
            num = tmpargs{1};
            allfields = fieldnames(EEG.event);
            if length( tmpargs ) < length(allfields)+1
                error('Pop_editeventvals: not enough arguments to change all field values');
            end;
            for index = 1:length( allfields )
                valstr = reformat(tmpargs{index+1}, strcmp(allfields{index}, 'latency'), EEG.trials > 1, num );
                eval([ 'EEG.event(' int2str(num) ').' allfields{index} '=' fastif(isempty(valstr), '[]', valstr) ';' ]);
	        end;
	end;
end;

% generate the output command
% ---------------------------
if exist('userdat') == 1
    if ~isempty(userdat)
        args = { args{:} userdat{:} };
    end;
end; 
com = sprintf('EEG = pop_editeventvals( %s', inputname(1));
for i=1:2:length(args)
    if iscell(args{i+1})
        com = sprintf('%s, ''%s'', {', com, args{i} );
        tmpcell = args{i+1};
        for j=1:length(tmpcell);
            if isstr( tmpcell{j} )   com = sprintf('%s ''%s'',', com, tmpcell{j} );
            else                     com = sprintf('%s [%s],',   com, num2str(tmpcell{j}) );
            end;
        end;
        com = sprintf('%s } ', com(1:end-1));     
    else
        com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
    end;       
end;
com = [com ');'];

return;

% format the output field
% -----------------------
function strval = reformat( val, latencycondition, trialcondition, eventnum)
    if latencycondition
        if trialcondition > 1
            strval = ['eeg_point2lat(' num2str(val) ', EEG.event(' int2str(eventnum) ').epoch, EEG.srate,[EEG.xmin EEG.xmax]*1000, 1E-3);' ];
        else    
            strval = [ '(' num2str(val) '-EEG.xmin)*EEG.srate+1;' ]; 
        end;
    else
        if isstr(val), strval = [ '''' val '''' ];
        else           strval = num2str(val);
        end;
    end;
