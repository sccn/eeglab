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
%   timewin   - Event time window [start, end] in seconds
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

% 02/15/02 modified function according to new event structure -ad

function epochval = eeg_getepochevent(EEG, type, timewin, fieldname, timeformat);

if nargin <2
    help eeg_getepochevent;
    return;
end;    
if nargin <3
    timewin = [-Inf Inf];
else if isempty(timewin)
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
if ~isfield( EEG.event, 'latency' )
    disp('Getepochevent: no ''latency'' field in events, abord'); return;
end;
if ~isfield( EEG.event, fieldname)
    disp(['Getepochevent: no ''' fieldname ''' field in events, abord']); return;
end;
EEG = eeg_checkset(EEG, 'eventconsistency');

% deal with empty types
% ---------------------
if isempty(type)
    type = unique( { EEG.event.type } );
end;

% convert types
% -------------
for indextype=1:length(type)
     if isstr(type{indextype})
         if ~isempty(str2num(type{indextype}))   type{indextype} = str2num(type{indextype}); end;
     end;
 end;

% select epochs
% -------------
epochval = zeros(1,EEG.trials); epochval(:) = nan;
for index = 1:length(EEG.event)
    epoch = EEG.event(index).epoch;

    for indextype=1:length(type)
        typeval = type{indextype};
        test = 0;
        if isstr(typeval) & isstr(EEG.event(index).type)
            if ~isempty(strmatch(typeval, EEG.event(index).type)), test = 1; end;
        elseif isnumeric(typeval) & isnumeric(EEG.event(index).type)
	        if typeval == EEG.event(index).type,  test = 1; end;
	    end;
	    if test           
            if isnan(epochval(epoch))
                switch fieldname,
                     case 'latency', 
                        epochval(epoch) = (mod(EEG.event(index).latency-1, EEG.pnts)/EEG.srate+EEG.xmin)*1000;
                     otherwise, eval( ['epochval(epoch) = EEG.event(index).' fieldname ';']);
                end;    
             else
                disp(['Getepochevent warning: more than one event selected in epoch ' int2str(epoch) ' -- only the field value for the first event returned']);
            end;
        end;    
    end;
end;
           
