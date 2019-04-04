% pop_chanevent() - import event latencies from the rising and/or falling 'edge' 
%                   latencies of a specified event-marker channel in EEG.data 
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
%   chanindices    - index|(indices) of the event channel(s)
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
%                    any function (Ex: 'myfunction(X)'). If an equal character
%                    is not present, this function preprend 'X=' to your 
%                    expression before evaluating it. Otherwise it just 
%                    evaluate the expression. For example, one may use
%                    'if X(1)>100,X(1)=0;end; X(find(X>100))=X(find(X>100)-1);'
%                    to remove spike artifacts in the data.
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
% Outputs:
%   OUTEEG         - EEGLAB output data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 29 July 2002
%
% See also: eeglab()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [EEG, command] = pop_chanevent(EEG, chan, varargin); 
command = '';

if nargin < 1
    help pop_chanevent;
    return;
end

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
			   { 'style' 'edit' 'string' '0' } { } ...
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
	
	if length(result) == 0 return; end
	chan   = eval( [ '[' result{1} ']' ] );
    options = {};
	if ~isempty(result{2}), options = { options{:} 'oper' result{2} }; end
	switch result{3},
		case 1, options = { options{:} 'edge' 'leading' };
		case 2, options = { options{:} 'edge' 'both' };
		case 3, options = { options{:} 'edge' 'trailing' };
	end; 
    options = { options{:} 'edgelen' eval( [ '[' result{4} ']' ] ) };
    if result{5},  options = { options{:} 'duration' 'on' }; end
	if ~result{6}, options = { options{:} 'delchan'  'off'}; end
	if ~result{7}, options = { options{:} 'delevent' 'off'}; end
	if result{8},  options = { options{:} 'nbtype'  1}; end
else
    options = varargin;
end
listcheck = { 'edge'      'string'     { 'both';'leading';'trailing'}     'both';
              'edgelen'   'integer'    [1 Inf]                            1;
              'delchan'   'string'     { 'on';'off' }                     'on';
              'oper'      'string'     []                                 '';
              'delevent'  'string'     { 'on';'off' }                     'on';
              'duration'  'string'     { 'on';'off' }                     'off';
              'typename'  'string'     []                                 [ 'chan' int2str(chan) ];
              'nbtype'    'integer'    [1 NaN]                             NaN };
g = finputcheck( options, listcheck, 'pop_chanedit');
if ischar(g), error(g); end

% check inut consistency
% ----------------------
if strcmpi(g.duration, 'on') && ~strcmpi(g.edge, 'leading')
    error('Must detect leading edge to extract event duration');
end

% process events
% --------------
fprintf('pop_chanevent: importing events from data channel %d ...\n', chan);
counte = 1; % event counter
events(10000).latency = 0;
if isnan(g.nbtype)
    if length(unique(EEG.data(chan, :))) == 2, g.nbtype = 1; end
end

for ci = chan
    X = EEG.data(ci, :);

    % apply preprocessing
    % -------------------
    if ~isempty(g.oper)
        if ~any( g.oper == '=' )
             g.oper = [ 'X = ' g.oper ';' ];
        else g.oper = [ g.oper ';' ];
        end
        try, eval( g.oper );
        catch, error('pop_chanevent: error executing preprocessing string');
        end
    end;    
    
    % extract edges
    % -------------
    tmpdiff =  diff(abs([ X X(end) ]));
    switch g.edge
     case 'both'    , tmpevent1 = find( tmpdiff > 0)-1; tmpevent2 = find( tmpdiff < 0);
     case 'trailing', tmpevent2 = find( tmpdiff < 0);
     case 'leading' , tmpevent1 = find( tmpdiff > 0)-1; tmpdur = find( tmpdiff < 0);
    end
    
    % fuse close events if necessary
    % ------------------------------
    if exist('tmpevent1')
        tmpclose = find( tmpevent1(2:end)-tmpevent1(1:end-1) < g.edgelen)+1;
        tmpevent1(tmpclose) = [];
        tmpevent    = tmpevent1+1;
        tmpeventval = tmpevent1+2;
    end
    if exist('tmpevent2')
        tmpclose = find( tmpevent2(2:end)-tmpevent2(1:end-1) < g.edgelen); % not +1
        tmpevent2(tmpclose) = [];
        tmpevent    = tmpevent2+1;
        tmpeventval = tmpevent2;
    end
    if exist('tmpevent1') && exist('tmpevent2')
        tmpevent    = sort([ tmpevent1+1 tmpevent2+1]);
        tmpeventval = sort([ tmpevent1+2 tmpevent2]);
    end
    
    % adjust edges for duration  if necessary
    % ---------------------------------------
    if strcmpi(g.duration, 'on')
        tmpclose = find( tmpdur(2:end)-tmpdur(1:end-1) < g.edgelen); % not +1 (take out the first)
        tmpdur(tmpclose) = [];
        if tmpdur(1) < tmpevent(1), tmpdur(1) = []; end
        if length(tmpevent) > length(tmpdur), tmpdur(end+1) = EEG.pnts; end
        if length(tmpevent) ~= length(tmpdur)
            error([ 'Error while attempting to extract event durations' 10 ...
                    'Maybe edges are not perfectly defined, try increasing edge length' ]);
        end
    end

    if isempty(tmpevent), 
        fprintf('No event found for channel %d\n', ci);
    else
        for tmpi = 1:length(tmpevent)
            if ~isnan(g.nbtype)
                events(counte).type    = g.typename;
            else
                events(counte).type    = X(tmpeventval(tmpi));
            end
            events(counte).latency = tmpevent(tmpi);
            if strcmpi(g.duration, 'on')
                events(counte).duration = tmpdur(tmpi) - tmpevent(tmpi);
            end
            counte = counte+1;
        end
    end
    events = events(1:counte-1);
end

% resort events
% --------------
if strcmp(g.delevent, 'on')
	EEG.event = events;
    if EEG.trials > 1
        for index = 1:length(events)
            EEG.event(index).epoch = 1+floor((EEG.event(index).latency-1) / EEG.pnts);
        end
    end
else
	for index = 1:length(events)
		EEG.event(end+1).type  = events(index).type;
		EEG.event(end).latency = events(index).latency;
        if EEG.trials > 1 || isfield(EEG.event, 'epoch');
            EEG.event(end).epoch = 1+floor((EEG.event(end).latency-1) / EEG.pnts);
        end
	end
    if EEG.trials > 1
        EEG = pop_editeventvals( EEG, 'sort', {  'epoch' 0 'latency', [0] } );
    else
        EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
    end
end
if isfield(EEG.event, 'urevent'), EEG.event = rmfield(EEG.event, 'urevent'); end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% delete channels
% ---------------
if strcmp(g.delchan, 'on')
	EEG = pop_select(EEG, 'nochannel', chan);
end

if nargin < 2
    command = sprintf('EEG = pop_chanevent(EEG, %s);', vararg2str({ chan options{:} })); 
end
return;
