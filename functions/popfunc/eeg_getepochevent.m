% eeg_getepochevent() - Return event field values for all events of one or more types
%
% Usage:
%       >> epochval = eeg_getepochevent( EEG, types);
%       >> epochval = eeg_getepochevent( EEG, types, timewin, fieldname);
%
% Inputs:
%   EEG       - Input dataset
%
% Optional inputs:
%   types     - Cell array of a subset of event types. 
%               {} is all type events. Note: this requires that 
%               a field named 'type' is defined in 'EEG.event'.
%   timewin   - Event time window [start, end] in milliseconds
%               (default []=whole epoch).
%   fieldname - Name of the field to return the values for. Default 
%               field is the event 'latency' in seconds
%               (though internally this is stored in frames)
% Outputs:
%   epochval  - A value of the selected field for each epoch; this is
%               NaN if no selected event occurred during the epoch.
%               Latencies are measured in msec relative to epoch onset.
%
% Notes: 1) Each epoch structure refers to the events that occurred
%        during its time window. This function allows the user to return 
%        specified field values for a subset of the defined events. 
%
%        2) If several of the selected events occur during a single epoch, 
%        a warning is issued, and value of ONLY THE FIRST event in the epoch 
%        is returned. 
%
%        If NO EVENT is selected in a given epoch, the value returned 
%        is NaN.
%
%        3) If the user elects to return the latency field, eeg_getepochevent()
%        recomputes the latency of each event relative to the epoch time
%        limits.
%
% Example: >> latencies = eeg_getepochevent(EEG, {'target','rare'}, [0 0.3]);
%          % Return the latencies (by default) of 'target' or 'rare' type
%          % events occurring between 0 and 0.3 sec of each epoch.
%          % Returns NaN for epochs with no such events. (See Notes above).
%
% Author: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 15 Feb 2002
%
% See also: eeglab(), epoch() 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 15 Feb 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.7  2002/05/03 01:53:54  arno
% using eeg_point2lat
%
% Revision 1.6  2002/04/22 22:05:53  arno
% debuggig last change
%
% Revision 1.5  2002/04/22 22:01:07  arno
% corrected time limits
%
% Revision 1.4  2002/04/22 21:51:49  arno
% removing error message for latency
%
% Revision 1.3  2002/04/18 18:22:27  arno
% typo can not
%
% Revision 1.2  2002/04/10 03:08:57  arno
% reprogrammed event selection
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 02/15/02 modified function according to new event structure -ad

function epochval = eeg_getepochevent(EEG, type, timewin, fieldname, timeformat);

if nargin <2
    help eeg_getepochevent;
    return;
end;    
if nargin <3
    timewin = [-Inf Inf];
else 
	if isempty(timewin)
        timewin = [-Inf Inf];
	end;
end;
if nargin <4
    fieldname = 'latency';
end;
if nargin <5
    timeformat = 'points';
end;

if isempty(EEG.event)
    disp('Getepochevent: no event structure, abord'); return;
end;
    
% check if EEG.epoch contain references to events
% -----------------------------------------------
if ~isfield( EEG.event, 'epoch' )
    disp('Getepochevent: no epoch indices in events, abord'); return;
end;
    
% check if EEG.epoch and EEG.event contains 'latency' field
% ------------------------------------------
if ~isfield( EEG.event, fieldname)
    disp(['Getepochevent: no ''' fieldname ''' field in events, abord']); return;
end;

% deal with empty types
% ---------------------
if ~isempty(type) & ~iscell(type)
	type = { type };
end;

% convert types
% -------------
for indextype=1:length(type)
     if isstr(type{indextype}) & isnumeric(EEG.event(1).type)
         if ~isempty(str2num(type{indextype}))   
			 type{indextype} = str2num(type{indextype}); 
		 else
			 error('eeg_getepochevent: string type cannot be found in numeric event type array');
		 end;		 
	 elseif isnumeric(type{indextype}) & isstr(EEG.event(1).type)
		  type{indextype} = num2str(type{indextype});
	 end;
end;

% select epochs
% -------------
if ~isempty(type)
	Ieventtmp = [];
	for indextype=1:length(type)
		typeval = type{indextype};
		test = 0;
		if isstr(typeval)
			Ieventtmp = [Ieventtmp strmatch(typeval, { EEG.event.type })' ];
		else
			Ieventtmp = [Ieventtmp find(typeval == cell2mat({ EEG.event.type })) ];
		end;
	end;
else
	Ieventtmp = [1:length(EEG.event)];
end;

% select latencies
% ----------------
if isfield(EEG.event, 'latency') & (timewin(1) ~= -Inf | timewin(2) ~= Inf)
	selected = ones(size(Ieventtmp));
	for index=1:length(Ieventtmp)
		reallat = eeg_point2lat(EEG.event(Ieventtmp(index)).latency, EEG.event(Ieventtmp(index)).epoch, ...
								EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3); 
		if reallat < timewin(1) | reallat > timewin(2)
			selected(index) = 0;
		end;
	end;
	Ieventtmp = Ieventtmp( find(selected == 1) );
end;

% select events
% -------------
epochval = zeros(1,EEG.trials); epochval(:) = nan;
if strcmp(fieldname, 'latency')
	for index = 1:length(Ieventtmp)
		epoch = EEG.event(Ieventtmp(index)).epoch;
		if isnan(epochval(epoch))
			epochval(epoch) = eeg_point2lat(EEG.event(Ieventtmp(index)).latency, epoch, EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
		else
			disp(['Getepochevent warning: more than one event selected in epoch ' int2str(epoch) ' -- only the field value for the first event returned']);
		end;
	end;
else
	for index = 1:length(Ieventtmp)
		eval( [ 'val = EEG.event(Ieventtmp(index)).' fieldname ';']);
		if ~isempty(val)
			epochval(EEG.event(Ieventtmp(index)).epoch) = val;
		end;
	end;
end;    
