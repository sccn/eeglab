% pop_epoch() - Convert a continuous EEG dataset to epoched data by extracting
%               data epochs time locked to specified event types or event indices. 
%               May also sub-epoch an already epoched dataset (if sub-epochs are 
%               same size or smaller). This pop_function calls epoch().
% Usage:
%   >> OUTEEG = pop_epoch( EEG); % pop-up a data entry window
%   >> OUTEEG = pop_epoch( EEG, events, timelimits);
%   >> [OUTEEG, indices] = pop_epoch( EEG, typerange, timelimits,'key1', value1 ...);
%
% Graphic interface:
%   "Time-locking event type(s)" - [edit box] Select 'Edit > Event values' 
%                to see a list of event.type values; else use the push button.
%                To use event types containing spaces, enter in single-quotes.
%                epoch() function command line equivalent: 'typerange' 
%   "..."      - [push button] scroll event types.
%   "Epoch limits" - [edit box] epoch latency range [start, end] in seconds relative
%                to the time-locking events. epoch() function equivalent: 'timelim' 
%   "Name for the new dataset" - [edit box] 
%                epoch() function equivalent: 'newname'
%   "Out-of-bounds EEG ..." - [edit box] Rejection limits ([min max], []=none).
%                epoch() function equivalent: 'valuelim' 
% Inputs:
%   EEG        - Input dataset. Data may already be epoched; in this case,
%                extract (shorter) subepochs time locked to epoch events.
%   typerange  - Cell array of event types to time lock to. 'eventindices'
%                {default {} --> time lock epochs to any type of event}
%                (Note: An event field called 'type' must be defined in the 
%                'EEG.event' structure. The command line argument is
%                'eventindices' below).
%   timelim    - Epoch latency limits [start end] in seconds relative to 
%                the time-locking event {default: [-1 2]}
%
% Optional inputs:
%   'eventindices'- [integer vector] Extract data epochs time locked to the 
%                indexed event numbers. 
%   'valuelim' - [min max] or [max]. Lower and upper bound latencies for 
%                trial data. Else if one positive value is given, use its 
%                negative as the lower bound. The given values are also 
%                considered outliers (min max) {default: none}
%   'verbose'  - ['yes'|'no'] {default: 'yes'}
%   'newname'  - [string] New dataset name {default: "[old_dataset] epochs"}
%   'epochinfo'- ['yes'|'no'] Propagate event information into the new
%                epoch structure {default: 'yes'}
%   
% Outputs:
%   OUTEEG     - output dataset
%   indices    - indices of accepted events
%
% Authors: Arnaud Delorme and Hilit Serby, SCCN, INC, UCSD, 2001
%
% See also: eeglab, epoch

% deprecated
%   'timeunit' - Time unit ['seconds'|'points'] If 'seconds,' consider events 
%                times to be in seconds. If 'points,' consider events as
%                indices into the data array. {default: 'points'}

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

% 01-25-02 reformated help & license -ad 
% 02-13-02 introduction of 'key', val arguments -ad
% 02-13-02 rereferencing of events -ad
% 03-18-02 interface and debugging -ad
% 03-27-02 interface and debugging -ad & sm

function [EEG, indices, com] = pop_epoch( EEG, events, lim, varargin )

if nargin < 1
   help pop_epoch;
	return;
end
com = '';
indices = [];

if length(EEG) > 1 && isempty(EEG(1).event)
    error('Cannot process empty event structure for the first dataset');
end

if isempty(EEG(1).event)
    if EEG.trials > 1 && EEG.xmin <= 0 && EEG.xmax >=0
        disp('No EEG.event structure found: creating events of type ''TLE'' (Time-Locking Event) at time 0');
        EEG.event(EEG.trials).epoch = EEG.trials; 
        for trial = 1:EEG.trials
            EEG.event(trial).epoch   = trial; 
            EEG.event(trial).type    = 'TLE';
            EEG.event(trial).latency = -EEG.xmin*EEG.srate+1+(trial-1)*EEG.pnts;
        end
    else
        disp('Cannot epoch data with no events'); beep;
        return;
    end
end
if ~isfield(EEG(1).event, 'latency'),
    error( 'Absent latency field in event array/structure: must name one of the fields ''latency''');
end
if size(EEG(1).data,3) > 1
    epochlim = [num2str( round(EEG.xmin)) '  '  num2str(round(EEG.xmax))]; % Units in seconds as in GUI
%   epochlim = [num2str( round(EEG.xmin*1000)) '  '  num2str(round(EEG.xmax*1000))]; % Units in miliseconds
else
    epochlim = '-1 2';
end
OLDEEG = EEG;

if nargin < 3
	% popup window parameters
	% -----------------------
   promptstr    = { strvcat('Time-locking event type(s) ([]=all):', ...
                    'Select ''Edit > Event values'' to see type values.'), ...
                    'Epoch limits [start, end] in seconds:', ... 
                    'Name for the new dataset:', ... 
					'Out-of-bounds EEG rejection limits ([min max], []=none):'  };

   cbevent = ['if ~isfield(EEG(1).event, ''type'')' ...
				   '   errordlg2(''No type field'');' ...
				   'else' ...
                   '   tmpevent = EEG(1).event;' ...
                   '   if isnumeric(EEG(1).event(1).type),' ...
				   '        [tmps,tmpstr] = pop_chansel(unique([ tmpevent.type ]));' ...
				   '   else,' ...
                   '        [tmps,tmpstr] = pop_chansel(unique({ tmpevent.type }));' ...
                   '   end;' ...
				   '   if ~isempty(tmps)' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''events''), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmps tmpevent tmpv tmpstr tmpfieldnames;' ];
   
   geometry = { [2 1 0.5] [2 1 0.5] [2 1.5] [2 1 0.5] };
   uilist = { { 'style' 'text'       'string' 'Time-locking event type(s) ([]=all)' } ...
              { 'style' 'edit'       'string' '' 'tag' 'events' } ...
              { 'style' 'pushbutton' 'string' '...' 'callback' cbevent } ...
              { 'style' 'text'       'string' 'Epoch limits [start, end] in seconds' } ...
              { 'style' 'edit'       'string' epochlim } ...
              { } ...
              { 'style' 'text'       'string' 'Name for the new dataset' } ...
              { 'style' 'edit'       'string'  fastif(length(EEG) > 1 || isempty(EEG(1).setname), '', [ EEG(1).setname ' epochs' ]) 'enable' fastif(length(EEG) == 1, 'on', 'off') } ...
              { 'style' 'text'       'string' 'Out-of-bounds EEG limits if any [min max]' } ...
              { 'style' 'edit'       'string' '' } { } };
              
   result = inputgui( geometry, uilist, 'pophelp(''pop_epoch'')', 'Extract data epochs - pop_epoch()');
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
   lim = eval( [ '[' result{2} ']' ] );

   args = {};
   if ~isempty( result{3} ),  args = { args{:}, 'newname', result{3} }; end
   if ~isempty( result{4} ),  args = { args{:}, 'valuelim', eval( [ '[' result{4} ']' ] ) }; end
   args = { args{:}, 'epochinfo', 'yes' };
 
else % no interactive inputs
    args = varargin;
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_epoch', EEG, 'warning', 'on', 'params', { events lim args{:} } );
    else
        [ EEG, com ] = eeg_eval( 'pop_epoch', EEG, 'params', { events lim args{:} } );
    end
    return;
end

% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('pop_epoch(): wrong syntax in function arguments'); return; end
else
    g = [];
end

% test the presence of variables
% ------------------------------
try, g.epochfield; 	 	  catch, g.epochfield = 'type'; end; % obsolete
try, g.timeunit; 	 	  catch, g.timeunit = 'points'; end
try, g.verbose; 	      catch, g.verbose = 'on'; end
try, g.newname; 	      catch, g.newname = fastif(isempty(EEG.setname), '', [EEG.setname ' epochs' ]); end
try, g.eventindices;      catch, g.eventindices = 1:length(EEG.event); end
try, g.epochinfo;         catch, g.epochinfo = 'yes'; end
try, if isempty(g.valuelim), g.valuelim = [-Inf Inf]; end; catch, g.valuelim = [-Inf Inf]; end

% transform string events into a int array of column indices
% ----------------------------------------------------------
tmpevent = EEG.event;
tmpeventlatency = [ tmpevent(:).latency ];
[tmpeventlatency Itmp] = sort(tmpeventlatency);
EEG.event = EEG.event(Itmp);  % sort by ascending time 
Ievent = g.eventindices;

if ~isempty( events )
    % select the events for epoching
    % ------------------------------
    Ieventtmp = [];
    tmpevent = EEG.event;
    tmpeventtype = { tmpevent.type };
    if ischar(tmpeventtype{1}), tmpeventtype  = deblank(tmpeventtype); end
    if iscell(events)
		if ischar(EEG.event(1).type)
			for index2 = 1:length( events )
				tmpevent = events{index2};
				if ~ischar( tmpevent ), tmpevent = num2str( tmpevent ); end
                tmpEventList = strmatch(deblank(tmpevent), tmpeventtype, 'exact');
				Ieventtmp = [ Ieventtmp ; tmpEventList(:) ];
			end
		else
			for index2 = 1:length( events )
				tmpevent = events{index2};
				if ischar( tmpevent ),tmpevent = str2num( tmpevent ); end
				if isempty( tmpevent ), error('pop_epoch(): string entered in a numeric field'); end
				Ieventtmp = [ Ieventtmp find(tmpevent == [ tmpeventtype{:} ]) ];
			end
		end
    else
        error('pop_epoch(): multiple event types must be entered as {''a'', ''cell'', ''array''}'); return;
    end
    Ievent = sort(intersect(Ievent, Ieventtmp));
end

% select event latencies for epoching
%------------------------------------
Ievent = sort(Ievent);
alllatencies = tmpeventlatency(Ievent);

if isempty(alllatencies)
   error('pop_epoch(): empty epoch range (no epochs were found).'); return;
end
fprintf('pop_epoch():%d epochs selected\n', length(alllatencies));

try
    % ----------------------------------------------------
    % For AMICA probabilities...Temporarily add model probabilities as channels
    %-----------------------------------------------------
    if isfield(EEG.etc, 'amica') && ~isempty(EEG.etc.amica) && isfield(EEG.etc.amica, 'v_smooth') && ~isempty(EEG.etc.amica.v_smooth) && ~isfield(EEG.etc.amica,'prob_added')
        if isfield(EEG.etc.amica, 'num_models') && ~isempty(EEG.etc.amica.num_models)
            if size(EEG.data,2) == size(EEG.etc.amica.v_smooth,2) && size(EEG.data,3) == size(EEG.etc.amica.v_smooth,3) && size(EEG.etc.amica.v_smooth,1) == EEG.etc.amica.num_models

                EEG = eeg_formatamica(EEG);
                %--------------------
                [EEG indices com] = pop_epoch(EEG,events,lim,args{:});
                %---------------------------------

                EEG = eeg_reformatamica(EEG);
                EEG = eeg_checkamica(EEG);
                return;
            else
                disp('AMICA probabilities not compatible with size of data, model probabilities cannot be epoched...')

            end
        end
    end
    % ----------------------------------------------------
catch
    warnmsg = strcat('your dataset contains amica information, but the amica plugin is not installed.  Continuing and ignoring amica information.');
    warning(warnmsg)
end

% change boundaries in rare cases when limits do not include time-locking events
% ------------------------------------------------------------------------------
tmpevents = EEG.event;
if lim(1) > 0 && ischar(EEG.event(1).type)
   % go through all onset latencies
   for Z1 = length(alllatencies):-1:1
      % if there is any event in between trigger and epoch onset which are boundary events
      selEvt = find([tmpevents.latency] > alllatencies(Z1) & [tmpevents.latency] < alllatencies(Z1) + lim(1) * EEG.srate);
      selEvt = selEvt(strcmp({tmpevents(selEvt).type}, 'boundary'));
      if any(selEvt)
          if sum([tmpevents(selEvt).duration]) > lim(1) * EEG.srate
              alllatencies(Z1) = [];
          else
              % correct the latencies by the duration of the data that were cutout
              alllatencies(Z1) = alllatencies(Z1) - sum([tmpevents(selEvt).duration]);
          end
      end
   end
end
if lim(2) < 0 && ischar(EEG.event(1).type)
   % go through all onset latencies
   for Z1 = length(alllatencies):-1:1
      % if there is any event in between trigger and epoch onset which are boundary events
      selEvt = find([tmpevents.latency] < alllatencies(Z1) & [tmpevents.latency] > alllatencies(Z1) + lim(2) * EEG.srate);
      selEvt = selEvt(strcmp({tmpevents(selEvt).type}, 'boundary'));
      if any(selEvt)
          if sum([tmpevents(selEvt).duration]) > -lim(2) * EEG.srate
              alllatencies(Z1) = [];
          else
              % correct the latencies by the duration of the data that were cutout
              alllatencies(Z1) = alllatencies(Z1) + sum([tmpevents(selEvt).duration]);
          end
      end
   end
end

% select event time format and epoch
% ----------------------------------
switch lower( g.timeunit )
 case 'points',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, [lim(1) lim(2)]*EEG.srate, ...
                                                    'valuelim', g.valuelim, 'allevents', tmpeventlatency);
  tmptime = tmptime/EEG.srate;
 case 'seconds',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, lim, 'valuelim', g.valuelim, ...
                                                    'srate', EEG.srate, 'allevents', tmpeventlatency);
	otherwise, disp('pop_epoch(): invalid event time format'); beep; return;
end
alllatencies = alllatencies(indices);
fprintf('pop_epoch():%d epochs generated\n', length(indices));


% update other fields
% -------------------
if lim(1) ~= tmptime(1) && lim(2)-1/EEG.srate ~= tmptime(2)
	fprintf('pop_epoch(): time limits have been adjusted to [%3.3f %3.3f] to fit data points limits\n', ...
		tmptime(1), tmptime(2)+1/EEG.srate);
end
EEG.xmin = tmptime(1);
EEG.xmax = tmptime(2);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.icaact = [];
if ~isempty(EEG.setname)
	if ~isempty(EEG.comments)
		EEG.comments = strvcat(['Parent dataset "' EEG.setname '": ----------'], EEG.comments);
	end
	EEG.comments = strvcat(['Parent dataset: ' EEG.setname ], ' ', EEG.comments);
end
EEG.setname = g.newname;

% count the number of events to duplicate and duplicate them
% ----------------------------------------------------------
totlen = 0;
for index=1:EEG.trials, totlen = totlen + length(epochevent{index}); end
EEG.event(1).epoch = 0;          % create the epoch field (for assignment consistency afterwards)
if totlen ~= 0
    newevent(totlen) = EEG.event(1); % reserve array
else 
    newevent = [];
end

% modify the event structure accordingly (latencies and add epoch field)
% ----------------------------------------------------------------------
allevents = [];
count = 1;
for index=1:EEG.trials
    for indexevent = epochevent{index}
		newevent(count)         = EEG.event(indexevent);
        newevent(count).epoch   = index;
	    newevent(count).latency = newevent(count).latency ...
			- alllatencies(index) - tmptime(1)*EEG.srate + 1 + EEG.pnts*(index-1);
		count = count + 1;
    end
end
EEG.event = newevent;
EEG.epoch = [];
EEG.saved = 'no';
EEG = eeg_checkset(EEG, 'eventconsistency');

% check for boundary events
% -------------------------
disp('pop_epoch(): checking epochs for data discontinuity');
if ~isempty(EEG.event) && ischar(EEG.event(1).type)
    tmpevent = EEG.event;
	boundaryindex = strmatch('boundary', { tmpevent.type });
	if ~isempty(boundaryindex)
		indexepoch = [];
		for tmpindex = boundaryindex
            if isfield(tmpevent, 'epoch')
    			indexepoch = [indexepoch tmpevent(tmpindex).epoch ];
            else 
                indexepoch = 1; % only one epoch
            end
		end
		EEG = pop_select(EEG, 'notrial', indexepoch);
        % update the "indices of accepted events", too
        indices = indices(setdiff(1:length(indices),indexepoch));
	end
end

% generate text command
% ---------------------
com = sprintf('EEG = pop_epoch( EEG, { ');
for j=1:length(events)
    if ischar( events{j} )  com = sprintf('%s ''%s'' ', com, events{j} );
    else                    com = sprintf('%s [%s] ',   com, num2str(events{j}) );
    end
end
com = sprintf('%s }, [%s]', com, num2str(lim)); 
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if ischar( args{i+1} )   com = sprintf('%s, ''%s'', ''%s''', com, args{i}, args{i+1} );
        else                    com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
        end
    end;    
end
com = [com ');'];
return; % text command

