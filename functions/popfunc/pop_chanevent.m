% pop_chanevent() - import event latencies from 'edge' values of a specified EEG.data channel 
%
% Usage:
%   >> OUTEEG = pop_chanevent( INEEG ); % select parameters via a pop-up window
%   >> OUTEEG = pop_chanevent( INEEG, chanindices, 'key', 'val' ... ); % no pop-up
%
% Graphic interface:
%   "Event channel(s)" - [edit box] indices of event channel(s) to import.
%                  Command line equivalent: chanindices.
%   "Preprocessing transform" - [edit box] apply this preprocessing
%                  formula or function to the selected data channel(s) X, 
%                  transforming X into the command output before edge
%                  extraction. Command line equivalent 'oper'.
%   "Transition to extract" - [list box] extract events when the event
%                  channel values go up ('leading'), down ('trailing')
%                  or both ('both'). Command line equivalent: 'edge'.
%   "Transition length" - [edit box] Increase this number to avoid having 
%                  events very close to each other due to a not perfectly 
%                  straight edge. Command line equivalent: 'edgelen'.
%   "Assign duration to events?" - [checkbox] . Assign duration to each 
%                  extracted event.  This option can only be used when 
%                  extracting events on leading edges. Event will last
%                  until next trailing edge (down) event. Command line 
%                  equivalent: 'duration'.
%   "Delete event channel(s)" - [checkbox] check to delete the event channel
%                  after events have been extracted from it.
%                  Command line equivalent: 'delchan'.
%   "Delete old events if any" - [checkbox] check this checkbox to 
%                  remove any prior events in the dataset. Otherwise 
%                  imported events are appended to old events. Command
%                  line equivalent: 'delevent'.
%   "Only one event type" - [checkbox] check this checkbox to assign
%                  all transitions in the event channel to one event 
%                  type. Else, one type is assigned for each non-zero
%                  channel value. Command line equivalent: 'nbtype'.
% Inputs:
%   INEEG          - input dataset structure
%   chanindices    - indices of an event channel(s)
%
% Optional inputs:
%   'edge'         - ['leading'|'trailing'|'both'] extract events when values
%                    in the event channel go up ('leading'), down ('trailing')
%                    or both ('both'). {Default is 'both'}.
%   'edgelen'      - [integer] maximum edge length (for some data edge do not
%                    take whole value and it takes a few sample points for
%                    signal to rise. Default is 1 (perfect edges).
%   'oper'         - [string] prior to extracting edges, preprocess data
%                    channel(s) using the string command argument.
%                    In this command, the data channel(s) are designated by
%                    (capital) X. For example, 'X>3' will test the value of X 
%                    at each time point (returning 1 if the data channel value 
%                    is larger than 3, and 0 otherwise). You may also use 
%                    any function (Ex: 'myfunction(X)').
%   'duration'     - ['on'|'off'] extract event duration. This option can only be
%                    used when extracting events on leading edges. Event will last
%                    until next trailing-edge (down) event { 'off' }.
%   'delchan'      - ['on'|'off'] delete channel from data { 'on' }.
%   'delevent'     - ['on'|'off'] delete old events if any { 'on' }.
%   'nbtype'       - [1|NaN] setting this to 1 will force the program to 
%                    consider all events to have the same type. {Default is NaN}.
%                    If set (1), all transitions are considered the same event type
%			         If unset (NaN), each (transformed) event channel value following a
%                    transition determines an event type (Ex: Detecting leading-edge
%                    transitions of 0 0 1 0 2 0 ...  produces event types 1 and 2).
%   'typename'     - [string] event type name. Only relevant if 'nbtype' is 1
%                    or if there is only one event type in the event channel.
%                    {Default is 'chanX', X being the index of
%                    the selected event channel}.
%
% Outputs:
%   OUTEEG         - EEGLAB output data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 29 July 2002
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.36  2004/06/16 16:09:05  arno
% new option edge length
%
% Revision 1.35  2004/05/24 17:17:27  arno
% extracting event duration
%
% Revision 1.34  2004/05/13 22:27:10  arno
% debug operation
%
% Revision 1.33  2004/04/16 15:40:39  arno
% nothing
%
% Revision 1.32  2004/03/05 23:25:54  arno
% debug gui order
%
% Revision 1.31  2004/03/04 18:21:26  arno
% more header
%
% Revision 1.30  2004/03/04 18:19:51  arno
% editing text
%
% Revision 1.29  2004/03/04 17:08:22  arno
% programming new field oper, suppressing threshold
%
% Revision 1.28  2004/03/04 03:09:58  arno
% implementing threshold, debug trials
%
% Revision 1.27  2004/01/15 17:55:46  scott
% same
%
% Revision 1.26  2004/01/15 17:53:34  scott
% same
%
% Revision 1.25  2004/01/15 17:52:38  scott
% same
%
% Revision 1.24  2004/01/15 17:50:23  scott
% pop-up window text
%
% Revision 1.23  2004/01/15 17:47:27  scott
% edit direction checkbox
%
% Revision 1.22  2004/01/15 17:46:03  scott
% edit help msg
%
% Revision 1.21  2004/01/15 17:37:00  scott
% edited pop-window text
%
% Revision 1.20  2003/12/11 20:24:53  arno
% nothing
%
% Revision 1.19  2003/11/18 22:50:22  arno
% do not know
%
% Revision 1.18  2003/09/22 23:30:30  arno
% [6~[6~nothing
%
% Revision 1.17  2003/09/22 23:28:03  arno
% nothing
%
% Revision 1.16  2003/06/19 16:09:28  arno
% make ur
%
% Revision 1.15  2003/04/10 17:30:33  arno
% header edit
%
% Revision 1.14  2002/12/06 03:49:35  arno
% resorint events
%
% Revision 1.13  2002/12/06 03:20:37  arno
% correcting typo
%
% Revision 1.12  2002/12/06 03:12:10  arno
% removing debuging message
%
% Revision 1.11  2002/12/06 02:52:12  arno
% inserting epoch number
%
% Revision 1.10  2002/12/06 02:44:41  arno
% debugging last
%
% Revision 1.9  2002/12/06 02:36:08  arno
% updating header and help
%
% Revision 1.8  2002/12/06 02:32:25  arno
% adding type name
%
% Revision 1.7  2002/10/09 22:28:47  arno
% update text
%
% Revision 1.6  2002/10/09 22:25:31  arno
% debugging
%
% Revision 1.5  2002/10/02 23:02:09  arno
% debug delevent and 'both' options
%
% Revision 1.4  2002/08/22 21:13:36  arno
% debug
%
% Revision 1.3  2002/08/06 21:39:11  arno
% spelling
%
% Revision 1.2  2002/07/29 17:57:12  arno
% debugging
%
% Revision 1.1  2002/07/29 17:53:03  arno
% Initial revision
%

function [EEG, command] = pop_chanevent(EEG, chan, varargin); 
command = '';

if nargin < 1
    help pop_chanevent;
    return;
end;

if nargin < 2
	geometry = { [1.5 1 1] [1] [1.5 1 1] [1.5 1 1] [1.5 1 1] [1.5 0.2 0.36 0.84] ...
                 [1] [1.5 0.21 1] [1.5 0.21 1] [1.5 0.21 1] };
    
    % callback from listbox to disable duration checkbox (if leading event is not selected)
    % --------------------------------------------------
    cb_list = [ 'if get(gcbo, ''value'') == 1,' ...
                '  set(findobj(gcbf, ''tag'', ''dur''), ''enable'', ''on'');' ...
                'else,' ...
                '  set(findobj(gcbf, ''tag'', ''dur''), ''enable'', ''off'', ''value'', 0);' ...
                'end;' ];
  
	strgui = { { 'style' 'text' 'string' 'Event channel(s)' 'tooltipstring' 'indexes of event channels' } ...
			   { 'style' 'edit' 'string' '' } { } ...
               {} ...
			   { 'style' 'text' 'string' 'Preprocessing transform (data=''X'')' 'tooltipstring'  ...
                 [ 'For example, ''X>3'' will test the value of X' 10 ...
                  'at each time point (returning 1 if the data channel value' 10 ...
                  'is larger than 3, and 0 otherwise).' ] } ...
			   { 'style' 'edit' 'string' '' } { 'style' 'text' 'string' 'Optional. Ex: X>3' } ...
			   { 'style' 'text' 'string' 'Transitions to extract? (up|down)' 'tooltipstring' ...
				 [ 'Extract events whenever values in the (transformed) event channel(s) shift up' 10 ...
				   '(''leading''), down (''trailing'') or either (''both'').' 10 ...
				   'AFTER SCROLLING CLICK TO SELECT' ] } ...
			   { 'style' 'listbox' 'string' 'up (leading)|both|down (trailing)' 'value' 1 'callback' cb_list } ...
               { 'style' 'text' 'string' '(click to select)'} ...
			   { 'style' 'text' 'string' 'Transition length (1=perfect edges)' 'tooltipstring'  ...
                 [ 'Increase this number to avoid having events very close to each other due.' 10 ...
                   'to a not perfectly straight edge' ] } ...
			   { 'style' 'edit' 'string' '1' } { } ...
			   { 'style' 'text' 'string' 'Assign duration to each events?' 'tag' 'dur' 'tooltipstring' ...
				 [ 'You may assign an event duration to each event if you select to detect' 10 ...
				   'event on the leading edge above. Event will last as long as the signal is non-0.' ] } ...
			   { 'style' 'checkbox' 'string' '' 'value' 0 'tag' 'dur'} { } ...
               { 'style' 'text' 'string' '(set=yes)' } ...
               {} ...
               { 'style' 'text' 'string' 'Delete event channel(s)? ' } ...
			   { 'style' 'checkbox' 'value' 1 } { 'style' 'text' 'string' '        (set = yes)'} ...
			   { 'style' 'text' 'string' 'Delete old events if any? ' } ...
			   { 'style' 'checkbox' 'value' 1 } { } ...
			   { 'style' 'text' 'string' 'All events of same type? ' 'tooltipstring' ...
			    ['If set, all transitions are considered the same event type,' 10 ...
				 'If unset, each (transformed) event channel value following a transition' 10 ...
                 'determines an event type (Ex: Detecting leading-edge transitions of' 10 ... 
                 '0 0 1 0 2 0 ...  produces event types 1 and 2).' ] } ...
			   { 'style' 'checkbox' 'value' 0 } { } };
	result       = inputgui( geometry, strgui, 'pophelp(''pop_chanevent'');', 'Extract event from channel(s) - pop_chanevent()');
	
	if length(result) == 0 return; end;
	chan   = eval( [ '[' result{1} ']' ] );
	if ~isempty(result{2}), g.oper = result{2}; else g.oper = ''; end;
	switch result{3},
		case 1, g.edge = 'leading';
		case 2, g.edge = 'both';
		case 3, g.edge = 'trailing';
	end;
    g.edgelen = eval( [ '[' result{4} ']' ] );
    if result{5}, g.duration = 'on'; else g.duration = 'off'; end;
	if result{6}, g.delchan  = 'on'; else g.delchan  = 'off'; end;
	if result{7}, g.delevent = 'on'; else g.delevent = 'off'; end;
	if result{8}, g.nbtype   = 1;    else g.nbtype   = NaN; end;
    g.typename =  [ 'chan' int2str(chan) ];
else 
	listcheck = { 'edge'      'string'     { 'both' 'leading' 'trailing'}     'both';
                  'edgelen'   'integer'    [1 Inf]                            1;
				  'delchan'   'string'     { 'on' 'off' }                     'on';
				  'oper'      'string'     []                                 '';
				  'delevent'  'string'     { 'on' 'off' }                     'on';
				  'duration'  'string'     { 'on' 'off' }                     'off';
                  'typename'  'string'     []                                 [ 'chan' int2str(chan) ];
				  'nbtype'    'integer'    [1 NaN]                             NaN };
	g = finputcheck( varargin, listcheck, 'pop_chanedit');
	if isstr(g), error(g); end;
end;

% check inut consistency
% ----------------------
if strcmpi(g.duration, 'on') & ~strcmpi(g.edge, 'leading')
    error('Must detect leading edge to extract event duration');
end;

% process events
% --------------
fprintf('pop_chanevent: importing events from data channel %d ...\n', chan);
counte = 1; % event counter
events(10000).latency = 0;
if isnan(g.nbtype)
    if length(unique(EEG.data(chan, :))) == 2, g.nbtype = 1; end;
end;

for ci = chan
    X = EEG.data(ci, :);

    % apply preprocessing
    % -------------------
    if ~isempty(g.oper)
        try, eval( [ 'X = ' g.oper ';' ]);
        catch, error('pop_chanevent: error executing preprocessing string');
        end;
    end;    
    
    % extract edges
    % -------------
    tmpdiff =  diff(X);
    switch g.edge
     case 'both'    , tmpevent1 = find( tmpdiff > 0)-1; tmpevent2 = find( tmpdiff < 0);
     case 'trailing', tmpevent2 = find( tmpdiff < 0);
     case 'leading' , tmpevent1 = find( tmpdiff > 0)-1; tmpdur = find( tmpdiff < 0);
    end;
    
    % fuse close events if necessary
    % ------------------------------
    if exist('tmpevent1')
        tmpclose = find( tmpevent1(2:end)-tmpevent1(1:end-1) < g.edgelen)+1;
        tmpevent1(tmpclose) = [];
        tmpevent    = tmpevent1+1;
        tmpeventval = tmpevent1+2;
    end;
    if exist('tmpevent2')
        tmpclose = find( tmpevent2(2:end)-tmpevent2(1:end-1) < g.edgelen); % not +1
        tmpevent2(tmpclose) = [];
        tmpevent    = tmpevent2+1;
        tmpeventval = tmpevent2;
    end;
    if exist('tmpevent1') & exist('tmpevent2')
        tmpevent    = sort([ tmpevent1+1 tmpevent2+1]);
        tmpeventval = sort([ tmpevent1+2 tmpevent2]);
    end;
    
    % adjust edges for duration  if necessary
    % ---------------------------------------
    if strcmpi(g.duration, 'on')
        tmpclose = find( tmpdur(2:end)-tmpdur(1:end-1) < g.edgelen); % not +1 (take out the first)
        tmpdur(tmpclose) = [];
        if tmpdur(1) < tmpevent(1), tmpdur(1) = []; end;
        if length(tmpevent) > length(tmpdur), tmpdur(end+1) = EEG.pnts; end;
        if length(tmpevent) ~= length(tmpdur)
            error([ 'Error while attempting to extract event durations' 10 ...
                    'Maybe edges are not perfectly defined, try increasing edge length' ]);
        end;
    end;
    
    if isempty(tmpevent), 
        fprintf('No event found for channel %d\n', ci);
    else
        for tmpi = 1:length(tmpevent)
            if ~isnan(g.nbtype)
                events(counte).type    = g.typename;
            else
                events(counte).type    = X(tmpeventval(tmpi));
            end;
            events(counte).latency = tmpevent(tmpi);
            if strcmpi(g.duration, 'on')
                events(counte).duration = tmpdur(tmpi) - tmpevent(tmpi);
            end;
            counte = counte+1;
        end;
    end;
    events = events(1:counte-1);
end;

% resort events
% --------------
if strcmp(g.delevent, 'on')
	EEG.event = events;
    if EEG.trials > 1
        for index = 1:length(events)
            EEG.event(index).epoch = 1+floor((EEG.event(index).latency-1) / EEG.pnts);
        end;
    end;
else
	for index = 1:length(events)
		EEG.event(end+1).type  = events(index).type;
		EEG.event(end).latency = events(index).latency;
        if EEG.trials > 1 | isfield(EEG.event, 'epoch');
            EEG.event(end).epoch = 1+floor((EEG.event(end).latency-1) / EEG.pnts);
        end;
	end;
    if EEG.trials > 1
        EEG = pop_editeventvals( EEG, 'sort', {  'epoch' 0 'latency', [0] } );
    else
        EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
    end;
end;
if isfield(EEG.event, 'urevent'), EEG.event = rmfield(EEG.event, 'urevent'); end;
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% delete channels
% ---------------
if strcmp(g.delchan, 'on')
	EEG = pop_select(EEG, 'nochannel', chan);
end;

command = sprintf('%s = pop_chanevent(%s, %s);', inputname(1), inputname(1), ...
				  vararg2str({ chan 'edge', g.edge 'nbtype' g.nbtype 'delchan' g.delchan 'delevent' g.delevent})); 
return;
