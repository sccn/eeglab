% pop_epoch() - epoching a continuous EEGLAB dataset. 
%
% Usage:
%   >> OUTEEG = pop_epoch( EEG, events, timelim);
%   >> [OUTEEG, indices] = ...
%        pop_epoch( EEG, typerange, timelim, 'key1', val1, 'key2', val1 ...);
%
% Inputs:
%   EEG        - input dataset. Dataset can already be epoched. The
%                epoch function can extract subepochs time locked on
%                an other event.
%   typerange  - cell array of events type to consider for epoching. 
%                {} is all type events. Note that this requires that 
%                field 'type' is defined in 'EEG.event'.
%   timelim    - [init end] in second for trial extraction centered
%                on the events (i.e. [-1 2])
%
% Optional inputs:
%   'timeunit'   - ['seconds'|'points'] if 'seconds', consider events in
%                time events in seconds. If 'points' consider events as
%                data array indexes. Default is 'points'   
%   'valuelim'   - upper and lower limit of values that a trial should not
%                overpass. If one positive value is given, consider the 
%                opposite for lower bound. Given values are also consider
%                outlier (if equal the trial is rejected). Default: none.
%   'verbose'    - ['yes'|'no']. Default is 'yes'.
%   'newname'    - 'string'. Default is (if parent dataset name not empty)
%                "Epoched from "[old dataset name]" dataset" 
%   'eventindices' - [indices], use event indices to epoch data
%   'epochinfo'  - ['yes'|'no']. Propagate all events information into the
%                the epoch structure. Default is 'yes'.
%   
% Outputs:
%   OUTEEG     - output dataset
%   indices    - indices of accepted events
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), epoch() 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.2  2002/04/10 02:42:22  arno
% debuging event selection
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 02-13-02 introduction of 'key', val arguments -ad
% 02-13-02 rereferencing of events -ad
% 03-18-02 interface and debugging -ad
% 03-27-02 interface and debugging -ad & sm

function [EEG, indices, com] = pop_epoch( EEG, events, lim, varargin );

com = '';
if nargin < 1
   help pop_epoch;
	return;
end;	

if isempty(EEG.event)
   disp('Cannot epoch data with no events'); beep;
   return;
end;
if ~isfield(EEG.event, 'latency'),
    disp( 'Absent latency field in event array/strcuture: first rename one of the field to ''latency''');
     beep; return;
end;    
OLDEEG = EEG;

if nargin < 3
	% popup window parameters
	% -----------------------
   promptstr    = { [ 'Enter time-locking event type(s) ([]=all):' 10 ...
                    '(use ''/Edit/Edit event values'' to scan type values)'], ...
                    'Epoch [start, end] in seconds (e.g. [-1 2]):', ... 
                    'Name of the new dataset:', ... 
					'Out-of-bounds EEG rejection limits, [min max] ([]=none):'  };

   inistr       = { '', '[-1 2]', fastif(isempty(EEG.setname), '', ['Epoched from "' EEG.setname '" dataset' ]), '' };
   result       = inputdlg( promptstr, 'Extract epochs -- pop_epochs', 1,  inistr);
   size_result  = size( result );
   if size_result(1) == 0 return; end;
   
   events = parsetxt( result{1} );
   lim = eval( [ '[' result{2} ']' ] );

   args = {};
   if ~isempty( result{3} ),  args = { args{:}, 'newname', result{3} }; end;
   if ~isempty( result{4} ),  args = { args{:}, 'valuelim', eval( [ '[' result{4} ']' ] ) }; end;
   args = { args{:}, 'epochinfo', 'yes' };
 
else % no interactive inputs
    args = varargin;
end;

% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('Setevent: wrong syntax in function arguments'); return; end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try, g.epochfield; 	 	  catch, g.epochfield = 'type'; end; % obsolete
try, g.timeunit; 	 	  catch, g.timeunit = 'points'; end;
try, g.verbose; 	      catch, g.verbose = 'on'; end;
try, g.newname; 	      catch, g.newname = fastif(isempty(EEG.setname), '', ['Epoched from "' EEG.setname '" dataset' ]); end;
try, g.eventindices;      catch, g.eventindices = []; end;
try, g.epochinfo;         catch, g.epochinfo = 'yes'; end;
try, if isempty(g.valuelim), g.valuelim = [-Inf Inf]; end; catch, g.valuelim = [-Inf Inf]; end;

% transform string events into a int array of column indices
% ----------------------------------------------------------
tmpeventlatency = cell2mat( { EEG.event(:).latency } );
Ievent = 1:length(EEG.event);    

if ~isempty( events )
    % select the events for epoching
    % ------------------------------
    Ieventtmp = [];
    tmpeventtype  = { EEG.event.type };
    if iscell(events)
		if isstr(EEG.event(1).type)
			for index2 = 1:length( events )
				tmpevent = events{index2};
				if ~isstr( tmpevent ), tmpevent = num2str( tmpevent ); end;
				Ieventtmp = [ Ieventtmp strmatch(tmpevent, tmpeventtype, 'exact') ];
			end;
		else
			for index2 = 1:length( tmpeventtype )
				tmpevent = events{index2};
				if isstr( tmpevent ),tmpevent = str2num( tmpevent ); end;
				if isempty( tmpevent ), error('Pop_epoch: string type in a numeric field'); end;
				Ieventtmp = [ Ieventtmp strmatch(tmpevent, tmpeventtype, 'exact') ];
			end;
		end;
    else
        error('Pop_epoch: types must be in cell array'); return;
    end;
    Ievent = Ieventtmp;
end;

% select event latencies for epoching
alllatencies = tmpeventlatency(Ievent);
alllatencies =sort(alllatencies(:)');

if isempty(alllatencies)
   error('Pop_epoch error: empty event range'); return;
end;
fprintf('Pop_epoch:%d epochs selected\n', length(alllatencies));

% select event time format and epoch
% ----------------------------------
switch lower( g.timeunit )
	case 'points',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, [lim(1)*EEG.srate lim(2)*EEG.srate], 'valuelim', g.valuelim, 'allevents', tmpeventlatency);
	case 'seconds',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, lim, 'valuelim', g.valuelim, 'srate', EEG.srate, 'allevents', tmpeventlatency);
	otherwise, disp('Pop_epoch error: invalid event time format'); beep; return;
end;
alllatencies = alllatencies(indices);
fprintf('Pop_epoch:%d epochs generated\n', length(indices));

% update other fields
% -------------------
EEG.xmin = lim(1);
EEG.xmax = lim(2)-1/EEG.srate;
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.icaact = [];
if ~isempty(EEG.setname)
	if ~isempty(EEG.comments)
		EEG.comments = strvcat(['Comment of parent dataset ' EEG.setname ': ----------'], EEG.comments);
	end;
	EEG.comments = strvcat(['Parent dataset: ' EEG.setname ], ' ', EEG.comments);
end;
EEG.setname = g.newname;

% modify the event structure accordingly (latencies and add epoch field)
% ----------------------------------------------------------------------
allevents = [];
newevent = EEG.event;
for index=1:EEG.trials
    allevents = union( allevents, epochevent{index});
    for indexevent = epochevent{index}
        newevent(indexevent).epoch   = index;
        switch lower( g.timeunit )
	       case 'points',	newevent(indexevent).latency = EEG.event(indexevent).latency - alllatencies(index) - lim(1)*EEG.srate + EEG.pnts*(index-1);
	       case 'seconds',  newevent(indexevent).latency = EEG.event(indexevent).latency - alllatencies(index) - lim(1) + EEG.pnts*(index-1)/EEG.srate;
        end;
    end;
end;
EEG.event = newevent;
EEG.event = EEG.event(allevents); % remove events that were not in the selection
EEG = eeg_checkset(EEG, 'eventconsistency');

% include the event latencies into EEG.epoch
% ------------------------------------------
maxlen = 0;
EEG.epoch(1).event = [];
for index = 1:length(EEG.event)
    currentepoch = EEG.event(index).epoch;
    if currentepoch <= length(EEG.epoch)
        EEG.epoch(currentepoch).event = [ EEG.epoch(currentepoch).event index ];
    else
        EEG.epoch(currentepoch).event = [ index ];
    end;
    maxlen = max(length(EEG.epoch(currentepoch).event), maxlen);
end;

% copy event information into the epoch array
% -------------------------------------------
switch lower(g.epochinfo)
    case 'yes', eventfields = fieldnames(EEG.event);
                for fieldnum = 1:length(eventfields)
                   for trial = 1:EEG.trials
                        %['EEG.epoch(trial).event' eventfields{fieldnum} '= { EEG.event(EEG.epoch(trial).event).' eventfields{fieldnum} '};' ]
                        if maxlen == 1, eval(['EEG.epoch(trial).event' eventfields{fieldnum} '= EEG.event(EEG.epoch(trial).event).' eventfields{fieldnum} ';' ]);
                        else            eval(['EEG.epoch(trial).event' eventfields{fieldnum} '= { EEG.event(EEG.epoch(trial).event).' eventfields{fieldnum} '};' ]);
                        end;
    			   end;
				end;    
end;

% generate text command
% ---------------------
com = sprintf('%s = pop_epoch( %s, { ', inputname(1), inputname(1));
for j=1:length(events);
    if isstr( events{j} )   com = sprintf('%s ''%s'' ', com, events{j} );
    else                    com = sprintf('%s [%s] ',   com, num2str(events{j}) );
    end;
end;
com = sprintf('%s }, [%s]', com, num2str(lim)); 
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if isstr( args{i+1} )   com = sprintf('%s, ''%s'', ''%s''', com, args{i}, args{i+1} );
        else                    com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
        end;
    end;    
end;
com = [com ');'];
return;

