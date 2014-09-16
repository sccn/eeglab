% eeg_getepochevent() - Return dataset event field values for all events 
%                                of one or more specified types
% Usage:
%       >> epochval = eeg_getepochevent( EEG );
%       >> epochval = eeg_getepochevent( EEG, 'key', 'val');
%
% Inputs:
%   EEG       - Input dataset
%
% Optional inputs:
%   'type'    - String containing an event type. Cell array of string
%               may be used to select several event types; 
%               {} is all types of events. Note: Requires that 
%               a field named 'type' is defined in 'EEG.event'.
%   'timewin' - [start end] Event time window in milliseconds
%               (default []=whole epoch).
%   'fieldname' - Name of the field to return the values for. 
%               Default field is 'EEG.event.latency' in milliseconds
%               (though internally this information is stored in 
%               real frames).
%   'trials'    - [integer array] return values only for selected trials.
%
% Outputs:
%   epochval    - A value of the selected field for each epoch. This is
%                 NaN if no selected event occurred during the epoch. If
%                 several values are available for each epoch, only the
%                 first one is taken into consideration.
%                 Latencies are measured in msec relative to epoch onset.
%                 Forced to be numerical, where a string is converted by 
%                 double to its ascii number which is normalized to be 
%                 between 0 and 1, and the string is summed together. See 
%                 the subfunction ascii2num for more details.
%   allepochval - cell array with same length as the number of epoch 
%                 containing all values for all epochs. This output is
%                 usefull when several value are found within each epoch.
%                 Not forced to be numerical.
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
% Example: 
%  >> latencies = eeg_getepochevent(EEG, 'rt');
%  % Return the latencies (by default) in milliseconds of events having 
%  % type 'rt' (reaction time)
%
%  >> latencies = eeg_getepochevent(EEG, {'target','rare'}, [0 300], 'position');
%  % Return the position (field 'position') of 'target' or 'rare' type
%  % events occurring between 0 and 300 milliseconds of each epoch.
%  % Returns NaN for epochs with no such events. (See Notes above).
%
% Author: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 15 Feb 2002
%
% See also: eeglab(), epoch() 

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

% 02/15/02 modified function according to new event structure -ad

function [epochval, allepochval] = eeg_getepochevent(EEG, varargin);
    
if nargin < 2
    help eeg_getepochevent;
    return;
end;

% process more than one EEG dataset (for STUDY purposes)
% ------------------------------------------------------
if length(EEG) > 1
    % the trial input may be a cell array; it has to be 
    % extracted before calling the function on each dataset
    trials  = cell(1,length(EEG));
    for iArg = length(varargin)-1:-2:1
        if strcmpi(varargin{iArg}, 'trials')
            trials = varargin{iArg+1};
            varargin(iArg:iArg+1) = [];
        end;
    end;
    
    epochval = [];
    for dat = 1:length(EEG)
        tmpepochval = eeg_getepochevent(EEG(dat), 'trials', trials{dat}, varargin{:});
        epochval = [ epochval tmpepochval ];
    end;
    return;
end;

% deal with old input format
% -------------------------
options = {};
oldformat = 0;
if nargin < 3
    oldformat = 1;
elseif isnumeric(varargin{2})  && length(varargin{2}) == 2
    oldformat = 1;
elseif length(varargin) == 3 && isfield(EEG.event,varargin(3))
    oldformat = 1;
end;
if oldformat
    if nargin > 1, options = { options{:} 'type'      varargin{1} }; end;
    if nargin > 2, options = { options{:} 'timewin'   varargin{2} }; end;
    if nargin > 3, options = { options{:} 'fieldname' varargin{3} }; end;
else
    options = varargin; 
end;
opt = finputcheck(options, { 'type'       { 'string';'cell' } { [] [] } '';
                             'timewin'    'real'              []        [-Inf Inf];
                             'fieldname'  'string'            []        'latency';
                             'trials'     { 'real';'cell' }   []        [] }, 'eeg_getepochevent');
if isstr(opt), error(opt); end;
if iscell(opt.trials) && ~isempty(opt.trials), opt.trials = opt.trials{1}; end;

if isempty(opt.timewin)
    opt.timewin = [-Inf Inf];
end;

if isempty(EEG.event)
    disp('Getepochevent: no event structure, aborting.'); return;
end;
    
% check if EEG.epoch contain references to events
% -----------------------------------------------
if ~isfield( EEG.event, 'epoch' )
    disp('Getepochevent: no epoch indices in events, considering continuous values.');
end;
    
% check if EEG.epoch and EEG.event contains 'latency' field
% ------------------------------------------
if ~isfield( EEG.event, opt.fieldname)
    disp(['Getepochevent: no ''' opt.fieldname ''' field in events, aborting.']); return;
end;

% deal with empty types
% ---------------------
if ~isempty(opt.type) & ~iscell(opt.type)
	opt.type = { opt.type };
end;

% convert types
% -------------
for indextype=1:length(opt.type)
     if isstr(opt.type{indextype}) & isnumeric(EEG.event(1).type)
         if ~isempty(str2num(opt.type{indextype}))   
			 opt.type{indextype} = str2num(opt.type{indextype}); 
		 else
			 error('eeg_getepochevent: string type cannot be found in numeric event type array');
		 end;		 
	 elseif isnumeric(opt.type{indextype}) & isstr(EEG.event(1).type)
		  opt.type{indextype} = num2str(opt.type{indextype});
	 end;
end;

% select epochs
% -------------
if ~isempty(opt.type)
	Ieventtmp = [];
    tmpevent  = EEG.event;
	for indextype=1:length(opt.type)
		typeval = opt.type{indextype};
		if isstr(typeval)
			Ieventtmp = [Ieventtmp strmatch(typeval, { tmpevent.type }, 'exact')' ];
		else
			Ieventtmp = [Ieventtmp find(typeval == [ tmpevent.type ] ) ];
		end;
	end;
else
	Ieventtmp = [1:length(EEG.event)];
end;

% select latencies
% ----------------
if isfield(EEG.event, 'latency') & (opt.timewin(1) ~= -Inf | opt.timewin(2) ~= Inf)
	selected = ones(size(Ieventtmp));
	for index=1:length(Ieventtmp)
        if ~isfield(EEG.event, 'epoch'), epoch = 1;
        else                             epoch = EEG.event(Ieventtmp(index)).epoch;
        end;
		reallat = eeg_point2lat(EEG.event(Ieventtmp(index)).latency, epoch, ...
								EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3); 
		if reallat < opt.timewin(1) | reallat > opt.timewin(2)
			selected(index) = 0;
		end;
	end;
	Ieventtmp = Ieventtmp( find(selected == 1) );
end;

% select events
% -------------
epochval       = cell(1,EEG.trials);  epochval(:) = { nan };
allepochval    = cell(1, EEG.trials); allepochval(:) = { {} };
if strcmp(opt.fieldname, 'latency')
	for index = 1:length(Ieventtmp)
        if ~isfield(EEG.event, 'epoch'), epoch = 1;
        else                             epoch = EEG.event(Ieventtmp(index)).epoch;
        end;
        
        allepochval{epoch}{end+1} = eeg_point2lat(EEG.event(Ieventtmp(index)).latency, epoch, ...
                                            EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
		if length(allepochval{epoch}) == 1
			epochval{epoch} = allepochval{epoch}{end};
		else
            if length(allepochval{epoch}) == 2 & nargout < 2
                disp(['Warning: multiple event latencies found in epoch ' int2str(epoch) ]); 
                %, ignoring event ' int2str(Ieventtmp(index)) ' (''' num2str(EEG.event(Ieventtmp(index)).type) ''' type)' ]);
            end;
		end;
	end;
elseif strcmp(opt.fieldname, 'duration')
	for index = 1:length(Ieventtmp)
		eval( [ 'val = EEG.event(Ieventtmp(index)).' opt.fieldname ';']);
		if ~isempty(val)
            if ~isfield(EEG.event, 'epoch'), epoch = 1;
            else                             epoch = EEG.event(Ieventtmp(index)).epoch;
            end;
            epochval{epoch}           = val/EEG.srate*1000;
            allepochval{epoch}{end+1} = val/EEG.srate*1000;
		end;
	end;
else
    for index = 1:length(Ieventtmp)
        eval( [ 'val = EEG.event(Ieventtmp(index)).' opt.fieldname ';']);
        if ~isempty(val)
            if isstr(val)
                val = ascii2num(val);
                %val_tmp = double(val);  % force epochval output to be numerical
                % **Turn string into number that will sort in alphebetical order**
                %val = 0;
                %for val_count = 1:length(val_tmp)
                    % -48 so that '1' is scaled to ascii number 1, not 49
                    % /74 to scale double('z')=122 to 1
                    % 10^((2-... scale to 0 to 100milliseconds
                    
                %    val = val + (val_tmp(val_count)-48)/74*10^(2-(val_count-1));
                %end
                % **End turn string ...**
            end
            if ~isfield(EEG.event, 'epoch'), 
                epoch = 1;
            else    epoch = EEG.event(Ieventtmp(index)).epoch;
            end
            epochval{epoch}           = val(1);
            allepochval{epoch}{end+1} = val(1);
		end;
	end;
end;    

if isnumeric(epochval{1})
    try 
        epochval = [ epochval{:} ];
        for index = 1:length(allepochval)
            allepochval{index} = [ allepochval{index}{:} ];
        end
    catch 
    end
end

% select specific trials
if ~isempty(opt.trials)
    epochval = epochval(opt.trials);
    allepochval = allepochval(opt.trials);
end;

%% SUBFUNCTION ASCII2NUM
% Maps ascii characters ['0','9'] to [1, 10], ['a','z'] to [11, 36]
% This is intended for alphebetically sorting string arrays numerically
%   ascii_in    [string array]
%   output      [double]
function out = ascii2num(ascii_in)

ascii_vector = double(ascii_in);
out = 0;
% go through each character in the string and scale and add it to output
for val_count = 1:length(ascii_vector)
    ascii_char = ascii_vector(val_count);
    if ascii_char>=48 & ascii_char<=57            % ['0','9'] to [1, 10]
        ascii_adj = ascii_char - 47;
    elseif ascii_char>=65 & ascii_char<=90        % ['A','Z'] to [11, 36]
        ascii_adj = ascii_char - 64;
    elseif ascii_char>=97 & ascii_char<=122       % ['a','z'] to [11, 36]
        ascii_adj = ascii_char - 96;
    else ascii_adj = ascii_char;
    end
    out = out + ascii_adj/36^val_count;
end






