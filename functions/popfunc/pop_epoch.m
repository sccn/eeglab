% pop_epoch() - Convert a continuous EEG dataset to epoched data by extracting
%               data epochs time locked to the specified events. Calls epoch().
% Usage:
%   >> OUTEEG = pop_epoch( EEG); % pop-up a data entry window
%   >> OUTEEG = pop_epoch( EEG, events, timelimits);
%   >> [OUTEEG, indices] = pop_epoch( EEG, typerange, timelimits, ...
%                               'key1', value1 ...);
%
% Graphic interface:
%   "Time-locking event type(s)" - [edit box] Select 'Edit > Event values' 
%                 to see type values or use the push button. epoch() 
%                 function command line equivalent: 'typerange' 
%   "..." - [push button] scroll event types.
%   "Epoch limits" - [edit box] time range [start, end] in seconds.
%                 epoch() equivalent: 'timelim' 
%   "Name for the new dataset" - [edit box] epoch() equivalent: 'newname'
%   "Out-of-bounds EEG ..." - [edit box] Rejection limits ([min max], []=none).
%                 epoch() equivalent: 'valuelim' 
%
% Inputs:
%   EEG        - Input dataset. Data may already be epoched. In this case,
%                extract (shorter) subepochs time locked to epoch events.
%   typerange  - Cell array of event types to time lock to.
%                {} --> all types of  events. Note: An event field 
%                called 'type' must be defined in the 'EEG.event' structure.
%   timelim    - Epoch limits [start end] in seconds relative to the
%                time-locking event {Default: [-1 2]}
%
% Optional inputs:
%   'valuelim' - [min max] Lower and upper bounds for trial data.  If one 
%                positive value is given, use its negative as the lower bound. 
%                The given values are also considered outliers. {Default: none}
%   'verbose'  - ['yes'|'no']. {Default: 'yes'}
%   'newname'  - [string] New dataset name {Default: "[old_dataset] epochs"}
%   'eventindices'- [indices] Extract data epochs time locked to the indexed events. 
%   'epochinfo'- ['yes'|'no']. Propagate event information into the new
%                epoch structure. {Default: 'yes'}
%   
% Outputs:
%   OUTEEG     - output dataset
%   indices    - indices of accepted events
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), epoch() 

% deprecated
%   'timeunit' - Time unit ['seconds'|'points'] If 'seconds,' consider events 
%                times to be in seconds. If 'points,' consider events as
%                indices into the data array. {Default: 'points'}

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
% Revision 1.40  2004/06/16 22:21:03  arno
% new GUI
%
% Revision 1.39  2004/05/14 23:56:29  arno
% remove debug msg
%
% Revision 1.38  2004/04/07 02:02:07  arno
% same
%
% Revision 1.37  2004/04/07 01:16:29  arno
% resorting events
%
% Revision 1.36  2003/11/18 16:30:18  scott
% same
%
% Revision 1.35  2003/11/18 16:29:31  scott
% Enter ... -> ...
%
% Revision 1.34  2003/08/19 22:07:12  arno
% automatically adding time-locking events if non
%
% Revision 1.33  2003/03/05 01:07:26  arno
% no error if types is "[]"
%
% Revision 1.32  2003/02/23 09:16:09  scott
% header edits -sm
%
% Revision 1.31  2003/02/17 02:40:14  arno
% reformating text for new functionality in help2html
%
% Revision 1.30  2003/02/16 22:53:02  arno
% update header for gui info
%
% Revision 1.29  2002/11/12 00:15:03  arno
% updating header
%
% Revision 1.28  2002/11/12 00:08:26  arno
% debugging eventindices
%
% Revision 1.27  2002/11/10 02:48:01  arno
% debugging g.eventindices
%
% Revision 1.26  2002/09/23 23:36:04  arno
% same
%
% Revision 1.25  2002/09/23 23:35:29  arno
% debug lat
%
% Revision 1.24  2002/09/23 23:22:44  arno
% cle[A
%
% display message
%
% Revision 1.23  2002/09/23 23:06:10  arno
% debugging limits
%
% Revision 1.22  2002/09/04 22:26:03  luca
% removing EEG.epoch mistmatch error message -arno
%
% Revision 1.21  2002/08/19 22:22:55  arno
% default output var
%
% Revision 1.20  2002/08/17 20:07:22  scott
% help message
%
% Revision 1.19  2002/08/17 19:53:02  scott
% Epoch [min,max] -> Epoch limits [min,max]
% ,
%
% Revision 1.18  2002/08/14 18:50:53  arno
% same
%
% Revision 1.17  2002/08/14 18:41:37  arno
% edit text
%
% Revision 1.16  2002/08/14 17:46:18  arno
% debuging
%
% Revision 1.15  2002/08/13 21:27:41  arno
% debug
%
% Revision 1.14  2002/08/12 16:28:13  arno
% inputdlg2
%
% Revision 1.13  2002/08/08 21:56:30  arno
% removing epoch creation
%
% Revision 1.12  2002/08/08 14:46:30  arno
% programming boundary events
%
% Revision 1.11  2002/08/06 21:35:13  arno
% spelling
%
% Revision 1.10  2002/07/21 17:17:15  arno
% debugging extract epochs from epochs
%
% Revision 1.9  2002/06/28 02:15:37  arno
% considering events in different epochs
%
% Revision 1.8  2002/06/25 00:50:17  arno
% debugging several event selection
%
% Revision 1.7  2002/06/25 00:46:39  arno
% adding eventconsistency check
%
% Revision 1.6  2002/06/25 00:38:03  arno
% sort event in ascending time before epoching
%
% Revision 1.5  2002/04/26 20:16:09  arno
% debugging epoch extraction from integer types
%
% Revision 1.4  2002/04/20 00:13:13  arno
% correcting latency computation bug
%
% Revision 1.3  2002/04/11 22:28:54  arno
% adding new dataset name
%
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

if nargin < 1
   help pop_epoch;
	return;
end;	
com = '';
indices = [];

if isempty(EEG.event)
    if EEG.trials > 1 & EEG.xmin <= 0 & EEG.xmax >=0
        disp('No event found: creating events of type ''TLE'' (Time-Locking Event) at time 0');
        EEG.event(EEG.trials).epoch = EEG.trials; 
        for trial = 1:EEG.trials
            EEG.event(trial).epoc    = trial; 
            EEG.event(trial).type    = 'TLE';
            EEG.event(trial).latency = -EEG.xmin*EEG.srate+1+(trial-1)*EEG.pnts;
        end;
    else
        disp('Cannot epoch data with no events'); beep;
        return;
    end;
end;
if ~isfield(EEG.event, 'latency'),
    disp( 'Absent latency field in event array/strcuture: first rename one of the field to ''latency''');
     beep; return;
end;    
OLDEEG = EEG;

if nargin < 3
	% popup window parameters
	% -----------------------
   promptstr    = { strvcat('Time-locking event type(s) ([]=all):', ...
                    'Select ''Edit > Event values'' to see type values.'), ...
                    'Epoch limits [start, end] in seconds:', ... 
                    'Name for the new dataset:', ... 
					'Out-of-bounds EEG rejection limits ([min max], []=none):'  };

   cbevent = ['if ~isfield(EEG.event, ''type'')' ...
				   '   errordlg2(''No type field'');' ...
				   'else' ...
				   '   if isstr(EEG.event(1).type)' ...
				   '       tmpfieldnames = unique({EEG.event.type});' ...
				   '   else' ...
				   '       tmpstr = unique(cell2mat({EEG.event.type}));' ...
				   '       tmpfieldnames = cell(1, length(tmpstr));' ...
                   '       for tmps=1:length(tmpstr), tmpfieldnames{tmps} = num2str(tmpstr(tmps)); end;' ...
				   '   end;' ...
		           '   [tmps,tmpv, tmpstr] = listdlg2(''PromptString'',''Select event types'', ' ...
                   '                                  ''ListString'', tmpfieldnames);' ...
				   '   if tmpv' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''events''), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmps tmpv tmpstr tmpfieldnames;' ];
   
   geometry = { [2 1 0.5] [2 1 0.5] [1.93 1.5] [2 1 0.5] };
   uilist = { { 'style' 'text'       'string' 'Time-locking event type(s) ([]=all)' } ...
              { 'style' 'edit'       'string' '' 'tag' 'events' } ...
              { 'style' 'pushbutton' 'string' '...' 'callback' cbevent } ...
              { 'style' 'text'       'string' 'Epoch limits [start, end] in seconds' } ...
              { 'style' 'edit'       'string' '[-1 2]' } ...
              { } ...
              { 'style' 'text'       'string' 'Name for the new dataset' } ...
              { 'style' 'edit'       'string'  fastif(isempty(EEG.setname), '', [ EEG.setname ' epochs' ]) } ...
              { 'style' 'text'       'string' 'Out-of-bounds EEG limits if any [min max]' } ...
              { 'style' 'edit'       'string' '' } { } };
              
   result = inputgui( geometry, uilist, 'pophelp(''pop_epoch'')', 'Extract data epochs - pop_epoch()');
   if length(result) == 0 return; end;
   
   if strcmpi(result{1}, '[]'), result{1} = ''; end;
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
   catch, disp('pop_epoch(): wrong syntax in function arguments'); return; end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try, g.epochfield; 	 	  catch, g.epochfield = 'type'; end; % obsolete
try, g.timeunit; 	 	  catch, g.timeunit = 'points'; end;
try, g.verbose; 	      catch, g.verbose = 'on'; end;
try, g.newname; 	      catch, g.newname = fastif(isempty(EEG.setname), '', [EEG.setname ' epochs' ]); end;
try, g.eventindices;      catch, g.eventindices = 1:length(EEG.event); end;
try, g.epochinfo;         catch, g.epochinfo = 'yes'; end;
try, if isempty(g.valuelim), g.valuelim = [-Inf Inf]; end; catch, g.valuelim = [-Inf Inf]; end;

% transform string events into a int array of column indices
% ----------------------------------------------------------
tmpeventlatency = cell2mat( { EEG.event(:).latency } );
[tmpeventlatency Itmp] = sort(tmpeventlatency);
EEG.event = EEG.event(Itmp);  % sort in ascending time 
Ievent = g.eventindices;

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
				Ieventtmp = [ Ieventtmp ; strmatch(tmpevent, tmpeventtype, 'exact') ];
			end;
		else
			for index2 = 1:length( events )
				tmpevent = events{index2};
				if isstr( tmpevent ),tmpevent = str2num( tmpevent ); end;
				if isempty( tmpevent ), error('pop_epoch(): string entered in a numeric field'); end;
				Ieventtmp = [ Ieventtmp find(tmpevent == cell2mat(tmpeventtype)) ];
			end;
		end;
    else
        error('pop_epoch(): multiple event types must be entered as {a  cell array}'); return;
    end;
    Ievent = Ieventtmp;
end;

% select event latencies for epoching
Ievent = sort(Ievent);
alllatencies = tmpeventlatency(Ievent);

if isempty(alllatencies)
   error('pop_epoch(): empty event range'); return;
end;
fprintf('pop_epoch():%d epochs selected\n', length(alllatencies));

% select event time format and epoch
% ----------------------------------
switch lower( g.timeunit )
 case 'points',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, [lim(1) lim(2)]*EEG.srate, 'valuelim', g.valuelim, 'allevents', tmpeventlatency);
  tmptime = tmptime/EEG.srate;
 case 'seconds',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, lim, 'valuelim', g.valuelim, 'srate', EEG.srate, 'allevents', tmpeventlatency);
	otherwise, disp('pop_epoch(): invalid event time format'); beep; return;
end;
alllatencies = alllatencies(indices);
fprintf('pop_epoch():%d epochs generated\n', length(indices));

% update other fields
% -------------------
if lim(1) ~= tmptime(1) & lim(2)-1/EEG.srate ~= tmptime(2)
	fprintf('pop_epoch(): time limits have been adjusted to [%3.3f %3.3f] to fit data points limits\n', tmptime(1), tmptime(2)+1/EEG.srate);
end;
EEG.xmin = tmptime(1);
EEG.xmax = tmptime(2);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.icaact = [];
if ~isempty(EEG.setname)
	if ~isempty(EEG.comments)
		EEG.comments = strvcat(['Parent dataset "' EEG.setname '": ----------'], EEG.comments);
	end;
	EEG.comments = strvcat(['Parent dataset: ' EEG.setname ], ' ', EEG.comments);
end;
EEG.setname = g.newname;

% count the number of events to duplicate and duplicate them
% ----------------------------------------------------------
totlen = 0;
for index=1:EEG.trials, totlen = totlen + length(epochevent{index}); end;
EEG.event(1).epoch = 0;          % create the epoch field (for assignment consistency afterwards)
if totlen ~= 0
    newevent(totlen) = EEG.event(1); % reserv array
else 
    newevent = [];
end;

% modify the event structure accordingly (latencies and add epoch field)
% ----------------------------------------------------------------------
allevents = [];
count = 1;
for index=1:EEG.trials
    for indexevent = epochevent{index}
		newevent(count)         = EEG.event(indexevent);
        newevent(count).epoch   = index;
	    newevent(count).latency = newevent(count).latency - alllatencies(index) - tmptime(1)*EEG.srate + 1 + EEG.pnts*(index-1);
		count = count + 1;
    end;
end;
EEG.event = newevent;
EEG.epoch = [];
EEG = eeg_checkset(EEG, 'eventconsistency');

% check for boundary events
% -------------------------
disp('pop_epoch(): checking epochs for data discontinuity');
if ~isempty(EEG.event) & isstr(EEG.event(1).type)
	boundaryindex = strmatch('boundary', { EEG.event.type });
	if ~isempty(boundaryindex)
		indexepoch = [];
		for tmpindex = boundaryindex
			indexepoch = [indexepoch EEG.event(tmpindex).epoch ];
		end;
		EEG = pop_select(EEG, 'notrial', indexepoch);
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

