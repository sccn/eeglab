% pop_adjustevent() - Adjust event offset of all (or specified) events
%
% Usage:
% >> EEG = pop_adjustevent(EEG); % launch a GUI
% >> EEG = pop_adjustevent(EEG,'key', val);
%
% Inputs:
%   EEG           - Input dataset
%
% Optional Inputs
%   'addms'       - [float] Positive or negative floating point number specifying
%                   the number of millisecond to add to each event latency
%   'addsamples'  - [float] Positive or negative floating point number specifying
%                   the number of samples to add to each event latency
%   'eventtype'   - [cell] Cell array of one or more strings giving the list 
%                   of event type(s) to apply the latency shift to. Default
%                   is all events.
%
% Outputs:
%   EEG           - Input dataset with latencies shifted
%   com           - Command for EEGLAB history
%
% Authors: Ramon Martinez-Cancino, Arnaud Delorme, Scott Makeig

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

function [EEG,com] = pop_adjustevent(EEG,varargin)
com = [];

% Call help
if nargin <1 
    help pop_adjustevent;
    return;
end

% GUI if inputs not provided
if nargin < 2
    srate = EEG.srate;
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
    
    uilist = { ...
        { 'style' 'text' 'string' 'Event type(s) to adjust (all by default): '} ...
        { 'style' 'edit' 'string'  ''   'tag' 'events' } ...
        { 'style' 'pushbutton' 'string' '...' 'callback' cbevent } ...
        { 'style' 'text' 'string' 'Add in milliseconds (can be negative)'} ...
        { 'style' 'edit' 'string'  '' 'tag' 'edit_time'} ...
        { }...
        { 'style' 'text' 'string' 'Or add in samples'} ...
        { 'style' 'edit' 'string'  '' 'tag' 'edit_samples'}...
        { }...
        };
        
    uigeom = { [1 0.7 0.5] [1 0.7 0.5] [1 0.7 0.5]};
    result = inputgui( 'uilist', uilist, 'geometry', uigeom, 'title', 'Adjust event latencies - pop_adjustevent()');
    if isempty(result) return; end
    
    % Collecting inputs
    options = { };
    if ~isempty(result{2})
        options = { options{:} 'addms' str2num(result{2}) };
    elseif ~isempty(result{3})
        options = { options{:} 'addsamples' str2num(result{3}) };
    end
    
    % Parse events
    eventTypes = result{1};
    if ~isempty(eventTypes)
        if strcmpi(eventTypes, '{}'), eventTypes = ''; end
        if ~isempty(eventTypes)
            if strcmpi(eventTypes(2),'''')
                g.eventtype = eval( [ '{' eventTypes '}' ] );
            else
                g.eventtype = parsetxt( eventTypes );
            end
        end
        options = { options{:} 'eventtypes' eventTypes };
    end
else
    options = varargin;
end

g = finputcheck(options, { 'addms'       'real'  []    [];
                           'addsamples'  'real'  []    [];
                           'eventtypes'  {'cell' 'string'}  []    {}});
if isstr(g)
    error(g);
end
if ~iscell(g.eventtypes)
    g.eventtypes = { g.eventtypes };
end

% Finding events indices
indx2shift = [];
if isempty(g.eventtypes)
    indx2shift = 1:length({EEG.event.type});
    
else
    for i=1:length(g.eventtypes)
        indxtmp = find(ismember({EEG.event.type}, g.eventtypes{i}));
        
        if isempty(indxtmp)
            warning(['Event type ''' g.eventtypes{i} ''' not found.']);
        end
        indx2shift = cat(2,indx2shift,indxtmp);
    end
end

% Checking eventnums and eventrange
if isempty(indx2shift) || any(isnan(indx2shift))
    error(['Unable to proceed. Event type(s) requested not found']);
end

% Shifting events
if ~isempty(g.addms)
    tmpevents = num2cell([EEG.event(indx2shift).latency]+g.addms/1000*EEG.srate);
    [EEG.event(indx2shift).latency] = tmpevents{:};
elseif ~isempty(g.addsamples)
    tmpevents = num2cell([EEG.event(indx2shift).latency]+g.addsamples);
    [EEG.event(indx2shift).latency] = tmpevents{:};
else
    error('To adjust event latencies, you need to specify a number of samples or ms');
end
    
% Checking new event latencies
EEG = eeg_checkset(EEG, 'eventconsistency');

if nargout > 1
    com = sprintf('[EEG,com] = pop_adjustevent(EEG, %s);', vararg2str(options));
end
