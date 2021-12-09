% pop_rmdat() - Remove continuous data around specific events
%
% Usage:
%   >> OUTEEG = pop_rmdat( EEG); % pop-up a data entry window
%   >> OUTEEG = pop_rmdat( EEG, typerange, timelimits, invertselection);
%
% Graphic interface:
%   "Time-locking event type(s)" - [edit box] Select 'Edit > Event values' 
%                to see a list of event.type values; else use the push button.
%                To use event types containing spaces, enter in single-quotes. 
%   "..."      - [push button] scroll event types.
%   "Time limits" - [edit box] epoch latency range [start, end] in seconds relative
%                to the event type latency.
% 
% Inputs:
%   EEG        - Input dataset. Data may already be epoched; in this case,
%                extract (shorter) subepochs time locked to epoch events.
%   typerange  - Cell array of event types to time lock to.
%                (Note: An event field called 'type' must be defined in 
%                the 'EEG.event' structure).
%   timelimits - Epoch latency limits [start end] in seconds relative to the time-locking event 
%                {default: [-1 2]}%   
%   invertselection - [0|1] Invert selection {default:0 is no}
%   
% Outputs:
%   OUTEEG     - output dataset
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2009-

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = pop_rmdat( EEG, events, timelims, invertsel )

if nargin < 1
   help pop_rmdat;
	return;
end
com = '';
invertsel = 0;

if isempty(EEG(1).event)
    error( [ 'No event. This function removes data' 10 'based on event latencies' ]);
end
if isempty(EEG(1).trials)
    error( [ 'This function only works with continuous data' ]);
end
if ~isfield(EEG(1).event, 'latency')
    error( 'Absent latency field in event array/structure: must name one of the fields ''latency''');
end

if nargin < 3
	% popup window parameters
	% -----------------------
   promptstr    = { strvcat('Time-locking event type(s) ([]=all):', ...
                    'Select ''Edit > Event values'' to see type values.'), ...
                    'Epoch limits [start, end] in seconds:', ... 
                    'Name for the new dataset:', ... 
					'Out-of-bounds EEG rejection limits ([min max], []=none):'  };

   cbevent = ['tmpEEG = get(gcbf, ''userdata'');' ...
           'if ~isfield(tmpEEG.event, ''type'')' ...
				   '   errordlg2(''No type field'');' ...
				   'else' ...
                   '   tmpevent = tmpEEG.event;' ...
                   '   if isnumeric(tmpEEG.event(1).type),' ...
				   '        [tmps,tmpstr] = pop_chansel(unique([ tmpevent.type ]));' ...
				   '   else,' ...
                   '        [tmps,tmpstr] = pop_chansel(unique({ tmpevent.type }));' ...
                   '   end;' ...
				   '   if ~isempty(tmps)' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''events''), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmpEEG tmps tmpv tmpstr tmpevent tmpfieldnames;' ];
   
   geometry = { [2 1 1.2] [2 1 1.2] };
   uilist = { { 'style' 'text'       'string' 'Event type(s) ([]=all)' } ...
              { 'style' 'edit'       'string' '' 'tag' 'events' } ...
              { 'style' 'pushbutton' 'string' '...' 'callback' cbevent } ...
              { 'style' 'text'       'string' 'Time limits [start, end] in sec.' } ...
              { 'style' 'edit'       'string' '-1 1' } ...
              { 'style' 'popupmenu'  'string' 'Keep selected|Remove selected' } };
              
   result = inputgui( 'geometry', geometry, 'uilist', uilist, 'helpcom', 'pophelp(''pop_rmdat'')', 'title', 'Remove data portions around events - pop_rmdat()', 'userdata', EEG);
   if length(result) == 0 return; end
   
   if strcmpi(result{1}, '[]'), result{1} = ''; end
   if ~isempty(result{1})
       if strcmpi(result{1}(1),'''')   % If event type appears to be in single-quotes, use comma
                                       % and single-quote as delimiter between event types. toby 2.24.2006
                                       % fixed Arnaud May 2006
            events = eval( [ '{' result{1} '}' ] );
       else events = parsetxt( result{1});
       end
   else events = {};
   end
   timelims = eval( [ '[' result{2} ']' ] );
   invertsel = result{3}-1;
 
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 3
        [ EEG, com ] = eeg_eval( 'pop_rmdat', EEG, 'warning', 'on', 'params', { events, timelims, invertsel } );
    else
        [ EEG, com ] = eeg_eval( 'pop_rmdat', EEG, 'params', { events, timelims, invertsel } );
    end
    return;
end

tmpevent = EEG.event;

% Checking for numeric values in fieldnames and changing them to string
 for i = 1: length(tmpevent)
     checknum(i) = isnumeric(tmpevent(i).type);
 end
  if sum(checknum)~=0
       tmpindx = find(checknum == 1);
       for i = 1: length(tmpindx)
           tmpevent(tmpindx(i)).type = num2str(tmpevent(tmpindx(i)).type);
       end
  end

% compute event indices
% ---------------------
allinds = [];
for index = 1:length(events)
    inds = strmatch(events{index},{ tmpevent.type }, 'exact');
    allinds = [allinds(:); inds(:) ]';
end
allinds = sort(allinds);
if isempty(allinds)
    disp('No event found');
    return;
end

% compute time limits
% -------------------
array = [];
bnd = strmatch('boundary', lower({tmpevent.type }));
bndlat = [ tmpevent(bnd).latency ];
for bind = 1:length(allinds)
    evtlat = EEG.event(allinds(bind)).latency;
    evtbeg = evtlat+EEG.srate*timelims(1);
    evtend = evtlat+EEG.srate*timelims(2);
    if any(bndlat > evtbeg & bndlat < evtend)
        % find the closer upper and lower boundaries
        bndlattmp = bndlat(bndlat > evtbeg & bndlat < evtend);
        diffbound = bndlattmp-evtlat;
        allneginds = find(diffbound < 0);
        allposinds = find(diffbound > 0);
        if ~isempty(allneginds), evtbeg = bndlattmp(allneginds(1)); end
        if ~isempty(allposinds), evtend = bndlattmp(allposinds(1)); end;       
        fprintf('Boundary found: time limits for event %d reduced from %3.2f to %3.2f\n', allinds(bind), ...
            (evtbeg-evtlat)/EEG.srate, (evtend-evtlat)/EEG.srate);
    end
    if ~isempty(array) && evtbeg < array(end)
        array(end) = evtend;
    else
        array = [ array; evtbeg  evtend];
    end
end

if ~isempty(array) && array(1) < 1, array(1) = 1; end
if ~isempty(array) && array(end) > EEG.pnts, array(end) = EEG.pnts; end
if isempty(array)
    disp('No event found');
    return;
end

if invertsel
    EEG = pop_select(EEG, 'notime', (array-1)/EEG.srate);
else
    EEG = pop_select(EEG, 'time', (array-1)/EEG.srate);
end;    
   
% generate output command
% -----------------------
com = sprintf('EEG = pop_rmdat( EEG, %s);', vararg2str( { events timelims invertsel } ));
