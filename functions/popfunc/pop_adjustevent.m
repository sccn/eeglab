% pop_adjustevent() - Adjust event offset of all (or specified) events
% Usage:
% >> EEG = pop_adjustevent(EEG); % launch a GUI
% >> EEG = pop_adjustevent(EEG,20);
% >> EEG = pop_adjustevent(EEG,20,'eventtype',{'rt'});
% >> EEG = pop_adjustevent(EEG,20,'eventtype',{'rt','square'});
% >> EEG = pop_adjustevent(EEG,20,'eventnums', [4 5 6 7 ]); % adjust events with these numbers
% >> EEG = pop_adjustevent(EEG,20,'eventrange', [5 20]);    % adjust events with numbers in this range
%
% Graphic interface:
%    "Offset Value (ms)"        - [edit box] See input 'offset'.
%    "Event type(s) ({}=all )"  - [edit box] See Optional input 'eventtype'
%
% Inputs:
%   EEG           - Input dataset
%   offset        - Positive or negative floating point number specifying the latency
%                    shift in [msec] to apply to the latencies of all events of the 
%                    specified type(s)   {default: 0}
% Optional Inputs
%   eventtype     - Cell array of one or more strings giving the event type(s) of
%                   the events to apply the latency shift to.
%                   {default: {} --> apply the latency shift to all events}
%   eventnums     - Integer array of specific event number(s) to which to apply
%                    the latency shift {default: [] --> all}
%   eventrange    - Integer vector [min, max] giving  a range of event numbers to
%                   which to apply the specified latency shift {default: [] --> all}
% Outputs:
%   EEG           - Input dataset with latencies shifted
%   com           - Command for EEGLAB history
%
%
% Authors: Ramon Martinez-Cancino, Arnaud Delorme, Scott Makeig
%
% See also: 
%
% Copyright (C) 2017 Ramon Martinez-Cancino, SCCN ramon@sccn.ucsd.edu 
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

function [EEG,com] = pop_adjustevent(EEG,offset,varargin);
com = [];

% Call help
if nargin <1 
    help pop_adjustevent;
    return;
end

% Options
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('pop_adjustevent() error: calling convention {''key'', value, ... } error'); return;
end;
try g.eventtype;       catch, g.eventtype   = {}';        end;
try g.eventnums;       catch, g.eventnums   = [];        end;
try g.eventrange;      catch, g.eventrange  = [];        end;

if exist('offset','var') && isnumeric(offset)
    g.offset = offset;
end

% GUI if inputs not provided
if nargin < 2
%     callback_forward  = ['set(findobj(''tag'',''checkbox_backward''), ''value'',~get(gcbo,''value''))'];
%     callback_backward = ['set(findobj(''tag'',''checkbox_forward''),  ''value'',~get(gcbo,''value''))'];
    srate = EEG.srate;
    cback_time   = ['set(findobj(''tag'', ''edit_samples''),    ''string'', num2str(str2num(get(gcbo,''string''))*' num2str(srate) '))'];
    cback_sample = ['set(findobj(''tag'', ''edit_time''),       ''string'', num2str(str2num(get(gcbo,''string''))/' num2str(srate) '))'];
    cbevent = ['if ~isfield(EEG.event, ''type'')' ...
        '   errordlg2(''No type field'');' ...
        'else' ...
        '   tmpevent = EEG.event;' ...
        '   if isnumeric(EEG.event(1).type),' ...
        '        [tmps,tmpstr] = pop_chansel(unique([ tmpevent.type ]));' ...
        '   else,' ...
        '        [tmps,tmpstr] = pop_chansel(unique({ tmpevent.type }));' ...
        '   end;' ...
        '   if ~isempty(tmps)' ...
        '       set(findobj(''parent'', gcbf, ''tag'', ''events''), ''string'', tmpstr);' ...
        '   end;' ...
        'end;' ...
        'clear tmps tmpevent tmpv tmpstr tmpfieldnames;' ];
    
    uilist = { { 'style' 'text' 'string' 'Offset Value : '} ...
        { 'style' 'edit' 'string'  ' ' 'tag' 'edit_time' 'callback' cback_time} ...
        { 'style' 'text' 'string' 'ms'}...
        { 'style' 'edit' 'string'  ' ' 'tag' 'edit_samples' 'callback' cback_sample}...
        { 'style' 'text' 'string' 'samples '} ...
        { 'style' 'text' 'string' 'Event type(s) ( {}=all ): '} ...
        { 'style' 'edit' 'string'  '{}'   'tag' 'events' } ...
        { 'style' 'pushbutton' 'string' '...' 'callback' cbevent } ...
        { }...
        { }...
        };
    
    addlistenerval = ['addlistener(findobj(''tag'',''edit_time''),''string'',''PostSet'',set(findobj(''tag'', ''edit_samples''),    ''string'', num2str(str2num(get(findobj(''tag'',''edit_time''),''string''))*' num2str(srate) ')))'];
    
    uigeom = { [1 0.7 0.5 0.7 0.5] [1 0.7 0.5 0.7 0.5]};
    result = inputgui( 'uilist', uilist, 'geometry', uigeom, 'title', 'Adjust event offset - pop_adjustevent()');
    if length(result) == 0 return; end;
    
    % Collecting inputs
    g.offset    = str2double(result{1});
    g.eventtype  = result{2} ;
    
    % --- From pop_epoch
   if strcmpi(result{3}, '{}'), result{3} = ''; end;
   if ~isempty(result{3})
       if strcmpi(result{3}(2),'''') 
            g.eventtype = eval( [ '{' result{3} '}' ] );
       else
           g.eventtype = parsetxt( result{3});
       end;
   else
       g.eventtype = {};
   end
   %---
end

% Starting here is the same for GUI or command line call

% Direction of offset
if sign(g.offset)== 1, g.direction = 'forward'; else g.direction = 'backward'; end;

% Finding events indices
indx2shift = [];
if isempty(g.eventtype) ||  any(strcmpi(g.eventtype, '[]'))      % all events
    indx2shift = 1:length({EEG.event.type});
    
elseif length(g.eventtype) == 1  % Only one event
    indx2shift =  find(ismember({EEG.event.type}, g.eventtype));
    
elseif length(g.eventtype) > 1   % Multiple events
    for i=1:length(g.eventtype)
        
        indxtmp = find(ismember({EEG.event.type}, g.eventtype{i}));
        
        if isempty(indxtmp)
            warning(['Event type ''' g.eventtype{i} ''' not found.']);
        end
        indx2shift = cat(2,indx2shift,indxtmp);
    end
end

% Checking eventnums and eventrange
flagevtnumrange = 0;
if ~isempty(g.eventnums)
    indx2shift = find(ismember(indx2shift, g.eventnums));
    flagevtnumrange = 1;
elseif ~isempty(g.eventrange) &&  length(g.eventrange)==2
    indx2shift = find(indx2shift >= g.eventrange(1) & indx2shift <= g.eventrange(2));
    flagevtnumrange = 2;
end

if isempty(indx2shift) || any(isnan(indx2shift))
    error(['Unable to proceed. Event type(s) requested not found']);
end

% Shifting events
switch g.direction
    case 'forward'
       tmpevents = num2cell([EEG.event(indx2shift).latency]+g.offset);
    case 'backward'
        tmpevents = num2cell([EEG.event(indx2shift).latency]-g.offset);
end
[EEG.event(indx2shift).latency] = tmpevents{:};

% Checking consistency of changes
EEG = eeg_checkset(EEG, 'eventconsistency');

% com output
if nargout > 1
    if isempty(g.eventtype)
        eventstring = '{';
    else
        eventstring = ['{''' g.eventtype{1} '''' ];
    end
    
    if length(g.eventtype)>1
        for ievents = 2:length(g.eventtype)
            eventstring = [eventstring ' , '''  g.eventtype{ievents} '''' ];
        end
    end
    
    if flagevtnumrange == 0
        evtnumstring   = '[]';        evtrangestring = '[]';
    elseif flagevtnumrange == 1
        evtnumstring   = ['[' num2str(g.eventnums) ']']; evtrangestring = '[]';
    elseif flagevtnumrange == 2
        evtnumstring   = '[]';        evtrangestring = ['[' num2str(g.eventrange) ']'];
    end
        
    com = ['[EEG,com] = pop_adjustevent(EEG,' num2str(g.offset) ',''eventtype'',' eventstring '},''eventnums'','  evtnumstring ',''eventrange'',' evtrangestring ')'];
end