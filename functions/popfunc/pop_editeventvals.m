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
% Revision 1.27  2004/06/12 01:50:09  arno
% still debuging
%
% Revision 1.26  2004/06/11 22:52:07  arno
% remove debug msg
%
% Revision 1.25  2004/05/27 02:01:42  arno
% update setfield for urevent
%
% Revision 1.24  2004/05/26 23:24:29  arno
% same
%
% Revision 1.23  2004/05/26 23:23:34  arno
% adding event consistency check
%
% Revision 1.22  2004/05/22 00:49:30  arno
% enable duration
%
% Revision 1.21  2004/05/21 21:21:49  arno
% debug history if no modification
%
% Revision 1.20  2003/08/29 19:05:19  arno
% first shot at urevent inserting ...
%
% Revision 1.19  2003/06/27 23:30:53  arno
% adding contextual help
%
% Revision 1.18  2003/06/27 23:12:56  arno
% new implementation including append and urevents
%
% Revision 1.17  2003/02/04 21:33:18  arno
% debugging command line call with empty values
%
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

if nargin >= 2 | isstr(EEG) % interpreting command from GUI or command line
    
    if isstr(EEG) % GUI
        gui    = 1;
        varargin = { EEG varargin{:} };
        
        % user data
        % ---------
        userdata  = get(gcf, 'userdata');
        EEG       = userdata{1};
        eventtmp  = userdata{2};
        oldcom    = userdata{3};
        allfields = fieldnames(eventtmp);
        tmpind    = strmatch('urevent', allfields);
        allfields(tmpind) = [];
        
        % current event
        % -------------
        objevent  = findobj('parent', gcf, 'tag', 'numval');
        valnum    = str2num(get(objevent, 'string'));
    
    else % command line    
        
        gui = 0;
        if isempty(EEG.event)
            disp('Getevent: cannot deal with empty event structure');
            return;
        end;   
        
        eventtmp  = EEG.event;
        allfields = fieldnames(EEG.event);
        tmpind = strmatch('urevent', allfields);
        allfields(tmpind) = [];
        
    end;
    
    % scan inputs
    % -----------
    for indfield = 1:2:length(varargin)
        if length(varargin) >= indfield+1
            tmparg = varargin{ indfield+1 };
        end;
        
        switch lower(varargin{indfield})
        
     case 'goto', % ******************** GUI ONLY ***********************
      
      % shift time
      % ----------
      shift     = tmparg;
      valnum    = valnum + shift;
      if valnum < 1,                valnum = 1;                end;
      if valnum > length(eventtmp), valnum = length(eventtmp); end;
      set(objevent, 'string', num2str(valnum,3));

      % update fields
      % -------------
      for index = 1:length(allfields) 
          
          enable = 'on';
          if strcmp( allfields{index}, 'latency') & ~isempty(eventtmp(valnum).latency)
              if isfield(eventtmp, 'type') & strcmpi(eventtmp(valnum).type, 'boundary'), enable = 'off'; end;
              if isfield(eventtmp, 'epoch')
                   value = eeg_point2lat( eventtmp(valnum).latency, eventtmp(valnum).epoch, ...
                                          EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
              else value = (eventtmp(valnum).latency-1)/EEG.srate+EEG.xmin;
              end;
          elseif strcmp( allfields{index}, 'duration') & ~isempty(eventtmp(valnum).duration)
              if isfield(eventtmp, 'type') & strcmpi(eventtmp(valnum).type, 'boundary'), enable = 'off'; end;
              if isfield(eventtmp, 'epoch')
                   value = eventtmp(valnum).duration/EEG.srate*1000; % milliseconds
              else value = eventtmp(valnum).duration/EEG.srate;      % seconds
              end;
          else
              value = getfield( eventtmp(valnum), allfields{index});
          end;
          
          % update interface
          % ----------------
          tmpobj = findobj('parent', gcf, 'tag', allfields{index});
          set(tmpobj, 'string', num2str(value,3), 'enable', enable);
      end;
      
      % update original
      % --------------- 
      tmpobj = findobj('parent', gcf, 'tag', 'original');
      if isfield(eventtmp, 'urevent') & eventtmp(valnum).urevent ~= valnum
           set(tmpobj, 'string', [ 'originally ' int2str(eventtmp(valnum).urevent)], ...
                       'horizontalalignment', 'center');
      else set(tmpobj, 'string', ' '); 
      end;
          
     case { 'append' 'insert' 'add' }, % **********************************************************

      if gui
          shift     = tmparg; % shift is for adding before or after the event
      else
          if strcmpi(lower(varargin{indfield}), 'insert')
               shift = 1;
          else shift = 0;
          end;
          valnum = tmparg{1};
      end;
      
      % find ur index
      % -------------
      if isfield(eventtmp, 'epoch'), curepoch = eventtmp(valnum).epoch; end;
      if isfield(EEG, 'urevent') & isfield(eventtmp, 'urevent')
          urvalnum = eventtmp(valnum).urevent;
          if isfield(EEG.urevent, 'epoch'), urcurepoch = EEG.urevent(urvalnum).epoch; end;
          urvalnum = urvalnum + shift;
      end;
      valnum    = valnum   + shift;
      
      % update urevents
      % ---------------
      if isfield(EEG, 'urevent') & isfield(eventtmp, 'urevent')
          EEG.urevent(end+3)            = EEG.urevent(end);
          EEG.urevent(urvalnum+1:end-2) = EEG.urevent(urvalnum:end-3);
          EEG.urevent(urvalnum)         = EEG.urevent(end-1);
          EEG.urevent                   = EEG.urevent(1:end-2);
          if isfield(EEG.urevent, 'epoch'),  EEG.urevent(urvalnum).epoch = urcurepoch; end;
      end;
      
      % update events
      % -------------
      eventtmp(end+3)          = eventtmp(end);
      eventtmp(valnum+1:end-2) = eventtmp(valnum:end-3);
      eventtmp(valnum)         = eventtmp(end-1);
      eventtmp                 = eventtmp(1:end-2);
      if isfield(eventtmp, 'epoch'), eventtmp(valnum).epoch = curepoch; end;      
      if isfield(EEG, 'urevent') & isfield(eventtmp, 'urevent')
          eventtmp(valnum).urevent = urvalnum;
          for index = valnum+1:length(eventtmp)
              eventtmp(index).urevent = eventtmp(index).urevent+1;
          end;
      end;
      
      % update type field
      % -----------------
      for tmpind = 1:length(allfields)
          eventtmp = checkconsistency(eventtmp, valnum, allfields{tmpind});
      end;
      
      if gui
          % update gui
          % ----------
          userdata{2} = eventtmp;
          set(gcf, 'userdata', userdata);
          pop_editeventvals('goto', shift);
          
          % update commands
          % ---------------
          tmpcell    = cell(1,1+length(fieldnames(eventtmp))); 
          tmpcell{1} = valnum;
          oldcom     = { oldcom{:} 'append', tmpcell };
      else
          % update field values if command line
          % -----------------------------------
          if any(~cellfun('isempty', tmparg(2:end)))
              EEG.event = eventtmp;
              for ind=2:length(tmparg)
                  if ind-1 <= length(allfields) & ~strcmpi(allfields{ind-1}, 'urevent') % do not include urevent 
                      EEG = pop_editeventvals(EEG, 'changefield', { valnum allfields{ind-1} tmparg{ind} });
                  end;
              end;
              eventtmp = EEG.event;
          end;
      end;
      
     case 'delete', % **********************************************************
      
      eventtmp(valnum) = []; 
      
      if gui, 
          
          valnum           = min(valnum,length(eventtmp));
          set(objevent, 'string', num2str(valnum)); 

          % update gui
          % ----------
          userdata{2} = eventtmp;
          set(gcf, 'userdata', userdata);
          pop_editeventvals('goto', 0);
          
          set(gcf, 'userdata', userdata);
          pop_editeventvals('goto', 0);
          
          % update commands
          % ---------------
          oldcom = { oldcom{:} 'delete', valnum };
      
      end;
    
     case { 'assign' 'changefield' }, % **********************************************************
      
      if gui, % GUI case
          
          field    = tmparg;
          objfield = findobj('parent', gcf, 'tag', field);
          editval     = get(objfield, 'string');
          if ~isempty(editval) & ~isempty(str2num(editval)), editval = str2num(editval); end;
      
      else % command line case
          
          valnum  = tmparg{1};
          field   = tmparg{2};
          editval = tmparg{3};
          
      end;
          
      % latency and duration case
      % -------------------------
      editvalori  = editval;
      if strcmp( field, 'latency') & ~isempty(editval)
          if isfield(eventtmp, 'epoch')
               editval = eeg_lat2point( editval, eventtmp(valnum).epoch, ...
                                       EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
          else editval = (editval- EEG.xmin)*EEG.srate+1;
          end;
      end;
      if strcmp( field, 'duration') & ~isempty(editval)
          if isfield(eventtmp, 'epoch')
               editval = editval/1000*EEG.srate; % milliseconds
          else editval = editval*EEG.srate;      % seconds
          end;
      end;
      
      % adapt to other formats
      % ----------------------
      eventtmp(valnum) = setfield(eventtmp(valnum), field, editval);
      eventtmp = checkconsistency(eventtmp, valnum, field);
      
      % update urevents
      % ---------------
      if isfield(EEG, 'urevent') & isfield(eventtmp, 'urevent')
          urvalnum  = eventtmp(valnum).urevent;
          
          % latency case
          % ------------
          if strcmp( field, 'latency') & ~isempty(editval)
              if isfield(EEG.urevent, 'epoch')
                  urepoch = EEG.urevent(urvalnum).epoch;
                  
                  
                  % find closest event latency
                  % --------------------------
                  if valnum<length(eventtmp)
                      if eventtmp(valnum+1).epoch == urepoch
                          urlatency = EEG.urevent(eventtmp(valnum+1).urevent).latency;
                          latency   = eventtmp(valnum+1).latency;
                      end;
                  end;
                  if valnum>1
                      if eventtmp(valnum-1).epoch == urepoch
                          urlatency = EEG.urevent(eventtmp(valnum-1).urevent).latency;
                          latency   = eventtmp(valnum-1).latency;
                      end;
                  end;
                  
                  % update event
                  % ------------
                  if exist('urlatency') ~=1
                      disp('Urevent not updated: could not find other event in the epoch');
                  else
                      editval = urlatency - ( latency - editval ); % new latency value
                  end;
              else 
                  editval = eeg_urlatency(EEG.event, eventtmp(valnum).latency);
              end;
          elseif strcmp( field, 'latency') % empty editval
              eventtmp(valnum).latency = NaN;
          end;
          
          % duration case
          % ------------
          if strcmp( field, 'duration') & ~isempty(editval)
              if isfield(eventtmp, 'epoch')
                   editval = editval/1000*EEG.srate; % milliseconds -> point
              else editval = editval*EEG.srate;      % seconds -> point
              end;
          end;
          
          EEG.urevent = setfield(EEG.urevent, {urvalnum}, field, editval);
      end;

      if gui
          % update history
          % --------------
          oldcom = { oldcom{:} 'changefield' { valnum field editvalori }};
      end;
      
     case 'sort', % **********************************************************
      
      if gui % retrieve data
          field1 = get(findobj('parent', gcf, 'tag', 'listbox1'), 'value');
          field2 = get(findobj('parent', gcf, 'tag', 'listbox2'), 'value');
          dir1   = get(findobj('parent', gcf, 'tag', 'order1'),   'value');
          dir2   = get(findobj('parent', gcf, 'tag', 'order2'),   'value');
          
          if field1 > 1, field1 = allfields{field1-1}; else return; end;
          if field2 > 1, field1 = allfields{field2-1}; else field2 = []; end;
          oldevents = EEG.event;
          EEG.event = eventtmp;
      else % command line
          field1 = tmparg{1};
          if length(tmparg) < 2, dir1 = 0;
          else                   dir1 = tmparg{2}; 
          end;
          if length(tmparg) < 3, field2 = [];
          else                   field2 = tmparg{3}; 
          end;
          if length(tmparg) < 4, dir2 = 0;
          else                   dir2 = tmparg{4}; 
          end;
      end;
      
      if ~isempty(field2)
          try, eval(['tmparray = cell2mat( { EEG.event.' field2 ' } );']);
          catch, eval(['tmparray = { EEG.event.' field2 ' };']);
          end;
          if strcmp(field2, 'latency') & EEG.trials > 1
              tmparray = eeg_point2lat(tmparray, {EEG.event.epoch}, EEG.srate, [EEG.xmin EEG.xmax], 1);
          end;
          [X I] = mysort( tmparray );
          if dir2 == 1, I = I(end:-1:1); end;
          events = EEG.event(I);
      else
          events = EEG.event;
      end;  
      try,   eval(['tmparray = cell2mat( { events.' field1 ' } );']);
      catch, eval(['tmparray = { events.' field1 ' };']);
      end;
      if strcmp( field1, 'latency') & EEG.trials > 1
          tmparray = eeg_point2lat(tmparray, {events.epoch}, EEG.srate, [EEG.xmin EEG.xmax], 1);
      end;
      [X I] = mysort( tmparray );
      if dir1 == 1, I = I(end:-1:1); end;
      EEG.event = events(I);
      
      if gui
          eventtmp  = EEG.event;
          EEG.event = oldevents;
          
          % update gui
          % ----------
          userdata{2} = eventtmp;
          set(gcf, 'userdata', userdata);
          pop_editeventvals('goto', 0);
          
          
          % update history
          % --------------
          oldcom = { oldcom{:} 'sort' { field1 dir1 field2 dir2 } };
          
          % warn user
          % ---------
          warndlg2('Sorting done');
      else 
          eventtmp = EEG.event;
      end;
      
    end; % end switch
    end; % end loop

    % save userdata
    % -------------
    if gui
        userdata{1} = EEG;
        userdata{2} = eventtmp;
        userdata{3} = oldcom;
        set(gcf, 'userdata', userdata);
    else
        EEG.event = eventtmp;
    end;
    return;
end;

% ----------------------
% graphic interface part
% ----------------------

if isempty(EEG.event)
    disp('Getevent: cannot deal with empty event structure');
    return;
end;   

allfields = fieldnames(EEG.event);
tmpind = strmatch('urevent', allfields);
allfields(tmpind) = [];

if nargin<2
    % transfer events to global workspace
    evalin('base', [ 'eventtmp = ' inputname(1) '.event;' ]);

    % add field values
    % ----------------
    geometry = { [2 0.5] };
    tmpstr = sprintf('Edit event field values (currently %d events)',length(EEG.event));
    uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } ...
               { 'Style', 'pushbutton', 'string', 'Delete event', 'callback', 'pop_editeventvals(''delete'');'  }};

    for index = 1:length(allfields) 

        geometry = { geometry{:} [1 1 1 1] };
        
        % input string
        % ------------
        if strcmp( allfields{index}, 'latency') | strcmp( allfields{index}, 'duration') 
            if EEG.trials > 1
                 inputstr =  [ allfields{index} ' (ms)'];
            else inputstr =  [ allfields{index} ' (sec)'];
            end;   
		else inputstr =  allfields{index};
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
		cbbutton = ['questdlg2(' vararg2str(tmptext) ...
					',''Description of field ''''' allfields{index} ''''''', ''OK'', ''OK'');' ];

        % create control
        % --------------
        cbedit = [ 'pop_editeventvals(''assign'', ''' allfields{index} ''');' ]; 
		uilist   = { uilist{:}, { }, ...
					 { 'Style', 'pushbutton', 'string', inputstr, 'callback',cbbutton  }, ...
					 { 'Style', 'edit', 'tag', allfields{index}, 'string', '', 'callback', cbedit } ...
                     { } };
    end;

    % add buttons
    % -----------
    geometry = { geometry{:} [1] [1.2 0.6 0.6 1 0.6 0.6 1.2] [1.2 0.6 0.6 1 0.6 0.6 1.2] [2 1 2] };
    
    tpappend = 'Append event before the current event';
    tpinsert = 'Insert event after the current event';
    tporigin = 'Original index of the event (in EEG.urevent table)';
    uilist   = { uilist{:}, ...
          { }, ...
          { },{ },{ }, {'Style', 'text', 'string', 'Event Num', 'fontweight', 'bold' }, { },{ },{ }, ...
          { 'Style', 'pushbutton', 'string', 'Append event',  'callback', 'pop_editeventvals(''append'', 0);', 'tooltipstring', tpappend }, ...
          { 'Style', 'pushbutton', 'string', '<<',            'callback', 'pop_editeventvals(''goto'', -10);' }, ...
          { 'Style', 'pushbutton', 'string', '<',             'callback', 'pop_editeventvals(''goto'', -1);' }, ...
          { 'Style', 'edit',       'string', '1',             'callback', 'pop_editeventvals(''goto'', 0);', 'tag', 'numval' }, ...
          { 'Style', 'pushbutton', 'string', '>',             'callback', 'pop_editeventvals(''goto'', 1);' }, ...
          { 'Style', 'pushbutton', 'string', '>>',            'callback', 'pop_editeventvals(''goto'', 10);' }, ...
          { 'Style', 'pushbutton', 'string', 'Insert event',  'callback', 'pop_editeventvals(''append'', 1);', 'tooltipstring', tpinsert }, ...
          { }, { 'Style', 'text',  'string', ' ', 'tag', 'original' 'horizontalalignment' 'center' 'tooltipstring' tporigin } { } };

    % add sorting options
    % -------------------
    listboxtext = 'No field selected';  
    for index = 1:length(allfields) 
         listboxtext = [ listboxtext '|' allfields{index} ]; 
    end;
    geometry = { geometry{:} [1] [1 1 1] [1 1 1] [1 1.5 0.5] };
    uilist = {  uilist{:}, ...
         { 'Style', 'text',       'string', 'Re-order events (for review only)', 'fontweight', 'bold'  }, ...
         { 'Style', 'text',       'string', 'Main sorting field:'  }, ...
         { 'Style', 'listbox',    'string', listboxtext, 'tag', 'listbox1' }, ...
         { 'Style', 'checkbox',   'string', 'Click for decreasing order', 'tag', 'order1' } ...
         { 'Style', 'text',       'string', 'Secondary sorting field:'  }, ...
         { 'Style', 'listbox',    'string', listboxtext, 'tag', 'listbox2' }, ...
         { 'Style', 'checkbox',   'string', 'Click for decreasing order', 'tag', 'order2' }, ...
         { 'Style', 'pushbutton', 'string', 'Re-sort', 'callback', 'pop_editeventvals(''sort'');' }, ...
         { }, { }};
   
    userdata = { EEG EEG.event {} };
    inputgui( geometry, uilist, 'pophelp(''pop_editeventvals'');', ...
                                  'Edit event values -- pop_editeventvals()', userdata, 'plot');
    pop_editeventvals('goto', 0);
    
    % wait for figure
    % ---------------
    fig = gcf;
    waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');
    try, userdata = get(fig, 'userdata'); close(fig); % figure still exist ?
    catch, return; end;
    
    % transfer events
    % ---------------
    if ~isempty(userdata{3})
        com = sprintf('%s = pop_editeventvals(%s,%s);', inputname(1), inputname(1), vararg2str(userdata{3}));
    end;
    if isempty(findstr('''sort''', com))
        if ~isempty(userdata{3}) % some modification have been done
            EEG       = userdata{1};
            EEG.event = userdata{2};
            disp('Checking event consistency...');
            TMPEEG = pop_editeventvals(EEG, 'sort', { 'latency' });
            if ~isequal(TMPEEG.event, EEG.event)
                EEG = TMPEEG;
                disp('Event resorted by increasing latencies. Some event indices have changed.');
            end;
            EEG = eeg_checkset(EEG, 'eventconsistency');
        end;
    else 
        com = '';
        disp('WARNING: all edits discarded because of event resorting. The EEGLAB event structure');
        disp('            must contain events sorted by latency (you may obtain an event structure');
        disp('            with resorted event by calling this function from the command line).');
    end;
    return;
    
end;
return;

% scan all the fields of g
% ------------------------
for curfield = 1:2:length(args)
    switch lower(args{curfield})
	   case 'delete'
	        EEG.event(args{ curfield+1 })=[];
	   case 'changefield'
            tmpargs = args{ curfield+1 };
            if length( tmpargs ) < 3
                error('Pop_editeventvals: not enough arguments to change field value');
            end;
            valstr = reformat(tmpargs{3}, strcmp(tmpargs{2}, 'latency'), EEG.trials > 1, tmpargs{1} );
            if strcmp(tmpargs{2}, 'duration'), 
                if EEG.trials > 1
                     valstr = num2str( tmpargs{3}/1000*EEG.srate ); % millisecond
                else valstr = num2str( tmpargs{3}*EEG.srate );      % second
                end;
            end;
            eval([ 'EEG.event(' int2str(tmpargs{1}) ').'  tmpargs{2} '=' fastif(isempty(valstr), '[]', valstr) ';' ]);
	   case { 'add' 'append' }
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

return;

% format the output field
% -----------------------
function strval = reformat( val, latencycondition, trialcondition, eventnum)
    if latencycondition
        if trialcondition
            strval = ['eeg_lat2point(' num2str(val) ', EEG.event(' int2str(eventnum) ').epoch, EEG.srate,[EEG.xmin EEG.xmax]*1000, 1E-3);' ];
        else    
            strval = [ '(' num2str(val) '-EEG.xmin)*EEG.srate+1;' ]; 
        end;
    else
        if isstr(val), strval = [ '''' val '''' ];
        else           strval = num2str(val);
        end;
    end;

% sort also empty values
% ----------------------
function [X, I] = mysort(tmparray);
    if iscell(tmparray)
        if all(cellfun('isreal', tmparray))
            tmpempty = cellfun('isempty', tmparray);
            tmparray(tmpempty) = { 0 };
            tmparray = cell2mat(tmparray);
        end;
    end;
    try, 
        [X I] = sort(tmparray);
    catch,
        sadf
        disp('Sorting failed. Check that selected fields contain uniform value format.');
        X = tmparray;
        I = 1:length(X);
    end;
    
% checkconsistency of new event
% -----------------------------
function eventtmp = checkconsistency(eventtmp, valnum, field)
    
    otherval = mod(valnum+1, length(eventtmp)+1)+1;
    
    if isstr(getfield(eventtmp(valnum), field)) & ~isstr(getfield(eventtmp(otherval), field))
        eventtmp(valnum) = setfield(eventtmp(valnum), field, str2num(getfield(eventtmp(valnum), field)));
    end;
    if ~isstr(getfield(eventtmp(valnum), field)) & isstr(getfield(eventtmp(otherval), field))
        eventtmp(valnum) = setfield(eventtmp(valnum), field, num2str(getfield(eventtmp(valnum), field)));
    end;
    if strcmpi(field, 'latency') & isempty(getfield(eventtmp(valnum), field))
        eventtmp(valnum).latency = NaN;
    end;
