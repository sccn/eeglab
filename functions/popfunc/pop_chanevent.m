% pop_chanevent() - import event latencies from 'edge' values of a specified EEG.data channel 
%
% Usage:
%   >> OUTEEG = pop_chanevent( INEEG ); % select parameters via a pop-up window
%   >> OUTEEG = pop_chanevent( INEEG, chanindices, 'key', 'val' ... ); % no pop-up
%
% Graphic interface:
%   "Event channel(s)" - [edit box] indices of event channel(s) to import.
%                  Command line equivalent: chanindices.
%   "Edge type to extract" - [list box] extract events when the event
%                  channel values go up ('leading'), down ('trailing')
%                  or both ('both'). Command line equivalent: 'edge'.
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
%   'delchan'      - ['on'|'off'] delete channel from data { 'on' }.
%   'delevent'     - ['on'|'off'] delete old events if any { 'on' }.
%   'nbtype'       - [1|NaN] setting this to 1 will force the program to 
%                    consider all events to have the same type. {Default is NaN}.
%   'typename'     - [string] event type name. Only relevant if 'nbtype' is 1
%                    or if there is only one event type in the event channel.
%                     {Default is 'chanX', X being the index of
%                     the selected event channel}.
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
	geometry = { [1.465 2.05] [1.5 1 1] [1.5 0.21 1] [1.5 0.21 1] [1.5 0.21 1] };
	strgui = { { 'style' 'text' 'string' 'Event channel(s)' 'tooltipstring' 'indexes of event channels' } ...
			   { 'style' 'edit' 'string' '' } ...
			   { 'style' 'text' 'string' 'Transitions to extract? (up|down)' 'tooltipstring' ...
				 [ 'Extract events when values in the event channel(s) shift up' 10 ...
				   '(''leading''), down (''trailing'') or either (''both'').' 10 ...
				   'IN UNIX, AFTER SCROLLING CLICK TO SELECT' ] } ...
			   { 'style' 'listbox' 'string' 'both|up (leading)|down (trailing)' 'value' 1 } { 'style' 'text' 'string' '(click to select)'} ...
			   { 'style' 'text' 'string' 'Delete event channel(s)? ' } ...
			   { 'style' 'checkbox' 'value' 1 } { 'style' 'text' 'string' '        (set = yes)'} ...
			   { 'style' 'text' 'string' 'Delete old events if any? ' } ...
			   { 'style' 'checkbox' 'value' 1 } { } ...
			   { 'style' 'text' 'string' 'All events of same type? ' 'tooltipstring' ...
			    ['If set, all transitions are considered the same,' 10 ...
				 'If unset, each signal value is assigned a different type'] } ...
			   { 'style' 'checkbox' 'value' 0 } { } };
	result       = inputgui( geometry, strgui, 'pophelp(''pop_chanevent'');', 'Extract event from channel(s) - pop_chanevent()');
	
	if length(result) == 0 return; end;
	chan   = eval( [ '[' result{1} ']' ] );
	switch result{2},
		case 1, g.edge = 'both';
		case 2, g.edge = 'leading';
		case 3, g.edge = 'trailing';
	end;
	if result{3}, g.delchan = 'on'; else g.delchan  = 'off'; end;
	if result{4}, g.delevent= 'on'; else g.delevent = 'off'; end;
	if result{5}, g.nbtype  = 1;     else g.nbtype   = NaN; end;
    g.typename =  [ 'chan' int2str(chan) ];
else 
	listcheck = { 'edge'     'string'     { 'both' 'leading' 'trailing'}     'both';
				  'delchan'  'string'     { 'on' 'off' }                     'on';
				  'delevent' 'string'     { 'on' 'off' }                     'on';
                  'typename' 'string'     []                                 [ 'chan' int2str(chan) ];
				  'nbtype'   'integer'    [1 NaN]                             NaN };
	g = finputcheck( varargin, listcheck, 'pop_chanedit');
	if isstr(g), error(g); end;
end;
if length(chan) ~= 1
	error('One (single) channel must be selected');
end;

% process events
% --------------
fprintf('pop_chanevent: importing events from data channel %d ...\n', chan);
counte = 1; % event counter
events(10000).latency = 0;
if isnan(g.nbtype)
    if length(unique(EEG.data(chan, :))) == 2, g.nbtype = 1; end;
end;

counttrial = 1;
switch g.edge
 case 'both'    , tmpevent = find( diff(EEG.data(chan, :)) ~= 0);
 case 'trailing', tmpevent = find( diff(EEG.data(chan, :)) < 0);
 case 'leading' , tmpevent = find( diff(EEG.data(chan, :)) > 0);
end;
if isempty(tmpevent), disp('No event found'); return; end;
tmpevent = tmpevent+1;
for tmpi = tmpevent
    if ~isnan(g.nbtype)
        events(counte).type    = g.typename;
    else
        events(counte).type    = EEG.data(chan, tmpi);
    end;
    events(counte).latency = tmpi;
    counte = counte+1;
end;
events = events(1:counte-1);

% resort events
% --------------
if strcmp(g.delevent, 'on')
	EEG.event = events;
else
	for index = 1:length(events)
		EEG.event(end+1).type  = events(index).type;
		EEG.event(end).latency = events(index).latency;
        if EEG.trials > 1
            EEG.event(end).epoch = 1+floor(EEG.event(end).latency / EEG.pnts);
        end;
	end;
    if EEG.trials > 1
        EEG = pop_editeventvals( EEG, 'sort', {  'epoch' 0 'latency', [0] } );
    else
        EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
    end;
end;
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
