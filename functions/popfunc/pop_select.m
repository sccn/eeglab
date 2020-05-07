% pop_select() - given an input EEG dataset structure, output a new EEG data structure 
%                retaining and/or excluding specified time/latency, data point, channel, 
%                and/or epoch range(s).
% Usage:
%   >> OUTEEG = pop_select(INEEG, 'key1', value1, 'key2', value2 ...);
%
% Graphic interface:
%   "Time range" - [edit box] RETAIN only the indicated epoch latency or continuous data 
%                  time range: [low high] in ms, inclusive. For continuous data, several 
%                  time ranges may be specified, separated by semicolons. 
%                  Example: "5 10; 12 EEG.xmax" will retain the indicated
%                  stretches of continuous data, and remove data portions outside
%                  the indicated ranges, e.g. from 0 s to 5 s and from 10 s to 12 s. 
%                  Command line equivalent: 'time' (or 'notime' - see below)
%   "Time range" - [checkbox] EXCLUDE the indicated latency range(s) from the data.
%                  For epoched data, it is not possible to remove a range of latencies 
%                  from the middle of the epoch, so either the low and/or the high values 
%                  in the specified latency range (see above) must be at an epoch boundary 
%                  (EEG.xmin, EEGxmax).  Command line equivalent: [if checked] 'notime' 
%   "Point range" - [edit box] RETAIN the indicated data point range(s). 
%                  Same options as for the "Time range" features (above).
%                  Command line equivalent: 'point' (or 'nopoint' - see below).
%   "Point range" - [checkbox] EXCLUDE the indicated point range(s).
%                  Command line equivalent: [if checked] 'nopoint' 
%   "Epoch range" - [edit box] RETAIN the indicated data epoch indices in the dataset.
%                  This checkbox is only visible for epoched datasets. 
%                  Command line equivalent: 'trial' (or 'notrial' - see below)
%   "Epoch range" - [checkbox] EXCLUDE the specified data epochs. 
%                   Command line equivalent: [if checked] 'notrial' 
%   "Channel range" - [edit box] RETAIN the indicated vector of data channels 
%                  Command line equivalent: 'channel' (or 'nochannel' - see below)
%   "Channel range" - [checkbox] EXCLUDE the indicated channels.
%                  Command line equivalent: [if checked] 'nochannel' 
%   "..." - [button] select channels by name.
%   "Scroll dataset" - [button] call the eegplot() function to scroll the
%                  channel activities in a new window for visual inspection.
%                  Commandline equivalent: eegplot() - see its help for details.
% Inputs:
%   INEEG         - input EEG dataset structure
%
% Optional inputs
%   'time'        - [min max] in seconds. Epoch latency or continuous data time range 
%                   to retain in the new dataset, (Note: not ms, as in the GUI text entry 
%                   above). For continuous data (only), several time ranges can be specified, 
%                   separated by semicolons. Example: "5 10; 12 EEG.xmax" will retain 
%                   the indicated times ranges, removing data  outside the indicated ranges 
%                   e.g. here from 0 to 5 s and from 10 s to 12 s. (See also, 'notime')
%   'notime'      - [min max] in seconds. Epoch latency or continuous dataset time range 
%                   to exclude from the new dataset. For continuous data, may be 
%                   [min1 max1; min2 max2; ...] to exclude several time ranges. For epoched 
%                   data, the latency range must include an epoch boundary, as latency 
%                   ranges in the middle of epochs cannot be removed from epoched data.
%   'point'       - [min max] epoch or continuous data point range to retain in the new 
%                   dataset. For continuous datasets, this may be [min1 max1; min2 max2; ...] 
%                   to retain several point ranges. (Notes: If both 'point'/'nopoint' and 
%                   'time' | 'notime' are specified, the 'point' limit values take precedence. 
%                   The 'point' argument was originally a point vector, now deprecated).
%   'nopoint'     - [min max] epoch or continuous data point range to exclude in the new dataset. 
%                   For epoched data, the point range must include either the first (0) 
%                   or the last point (EEG.pnts), as a central point range cannot be removed. 
%   'trial'       - array of trial indices to retain in the new dataset
%   'notrial'     - array of trial indices to exclude from the new dataset
%   'sorttrial'   - ['on'|'off'] sort trial indices before extracting them (default: 'on').
%   'channel'     - vector of channel indices to retain in the new 
%                   dataset. Can also be a cell array of channel names.
%   'nochannel'   - vector of channel indices to exclude from the new
%                   dataset. Can also be a cell array of channel names.
%
% Outputs:
%   OUTEEG        - new EEG dataset structure
%
% Note: This function performs a conjunction (AND) of all its optional inputs.
%       Using negative counterparts of all options, any logical combination is
%       possible.
% 
% Author: Arnaud Delorme, CNL/Salk Institute, 2001; SCCN/INC/UCSD, 2002-
% 
% see also: eeglab()

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
% 01-26-02 changed the format for events and trial conditions -ad
% 02-04-02 changed display format and allow for negation of inputs -ad 
% 02-17-02 removed the event removal -ad 
% 03-17-02 added channel info subsets selection -ad 
% 03-21-02 added event latency recalculation -ad 

function [EEG, com] = pop_select( EEG, varargin);

com = '';
if nargin < 1
    help pop_select;
    return;
end
if isempty(EEG(1).data)
    disp('Pop_select error: cannot process empty dataset'); return;
end
    
if nargin < 2
   geometry = { [1 1 1] [1 1 0.25 0.23 0.51] [1 1 0.25 0.23 0.51] [1 1 0.25 0.23 0.51] ...
           [1 1 0.25 0.23 0.51] [1] [1 1 1]};
   uilist = { ...
         { 'Style', 'text', 'string', 'Select data in:', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Input desired range', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'on->remove these', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Time range [min max] (s)', 'fontangle', fastif(length(EEG)>1, 'italic', 'normal') }, ...
         { 'Style', 'edit', 'string', '', 'enable', fastif(length(EEG)>1, 'off', 'on') }, ...
         { }, { 'Style', 'checkbox', 'string', '    ', 'enable', fastif(length(EEG)>1, 'off', 'on') },{ }, ...
         ...
         { 'Style', 'text', 'string', 'Point range (ex: [1 10])', 'fontangle', fastif(length(EEG)>1, 'italic', 'normal') }, ...
         { 'Style', 'edit', 'string', '', 'enable', fastif(length(EEG)>1, 'off', 'on') }, ...
         { }, { 'Style', 'checkbox', 'string', '    ', 'enable', fastif(length(EEG)>1, 'off', 'on') },{ }, ...
         ...
         { 'Style', 'text', 'string', 'Epoch range (ex: 3:2:10)', 'fontangle', fastif(length(EEG)>1, 'italic', 'normal') }, ...
         { 'Style', 'edit', 'string', '', 'enable', fastif(length(EEG)>1, 'off', 'on') }, ...
         { }, { 'Style', 'checkbox', 'string', '    ', 'enable', fastif(length(EEG)>1, 'off', 'on') },{ }, ...
         ...
         { 'Style', 'text', 'string', 'Channel range' }, ...
         { 'Style', 'edit', 'string', '', 'tag', 'chans' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ' }, ...
         { 'style' 'pushbutton' 'string'  '...', 'enable' fastif(isempty(EEG(1).chanlocs), 'off', 'on') ...
           'callback' 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''chans''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' }, ...
           { }, { }, { 'Style', 'pushbutton', 'string', 'Scroll dataset', 'enable', fastif(length(EEG)>1, 'off', 'on'), 'callback', ...
                          'eegplot(EEG.data, ''srate'', EEG.srate, ''winlength'', 5, ''limits'', [EEG.xmin EEG.xmax]*1000, ''position'', [100 300 800 500], ''xgrid'', ''off'', ''eloc_file'', EEG.chanlocs);' } {}};
   results = inputgui( geometry, uilist, 'pophelp(''pop_select'');', 'Select data -- pop_select()' );
   if length(results) == 0, return; end

   
   % decode inputs
   % -------------
   args = {};
   if ~isempty( results{1} )
       if ~results{2}, args = { args{:}, 'time', eval( [ '[' results{1} ']' ] ) };
       else            args = { args{:}, 'notime', eval( [ '[' results{1} ']' ] ) }; end
   end

   if ~isempty( results{3} )
       if ~results{4}, args = { args{:}, 'point', eval( [ '[' results{3} ']' ] ) };
       else            args = { args{:}, 'nopoint', eval( [ '[' results{3} ']' ] ) }; end
   end

   if ~isempty( results{5} )
       if ~results{6}, args = { args{:}, 'trial', eval( [ '[' results{5} ']' ] ) };
       else            args = { args{:}, 'notrial', eval( [ '[' results{5} ']' ] ) }; end
   end

   if ~isempty( results{7} )
       [ chaninds chanlist ] = eeg_decodechan(EEG(1).chanlocs, results{7});
       if isempty(chanlist), chanlist = chaninds; end
       if ~results{8}, args = { args{:}, 'channel'  , chanlist };
       else            args = { args{:}, 'nochannel', chanlist }; end
   end

else
    args = varargin;
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_select', EEG, 'warning', 'on', 'params', args);
    else
        [ EEG, com ] = eeg_eval( 'pop_select', EEG, 'warning', 'off', 'params', args);
    end
    return;
end

%----------------------------AMICA---------------------------------
if isfield(EEG.etc,'amica') && isfield(EEG.etc.amica,'prob_added')
    for index = 1:2:length(args)
       if strcmpi(args{index}, 'channel')
           args{index+1} = [ args{index+1} EEG.nbchan-(0:2*EEG.etc.amica.num_models-1)];
           
       end
       
       
    end
end
%--------------------------------------------------------------------
        
if isempty(EEG.chanlocs), chanlist = [1:EEG.nbchan];
else                      chanlocs = EEG.chanlocs; chanlist = { chanlocs.labels };
end
g = finputcheck(args, { 'time'    'real'      []         []; ...
                        'notime'  'real'      []         []; ...
                        'trial'   'integer'   []         [1:EEG.trials]; ...
                        'notrial' 'integer'   []         []; ...
                        'point'   'integer'   []         []; ...
                        'nopoint' 'integer'   []         []; ...
                        'channel'   { 'integer','cell' }  []  chanlist;
                        'nochannel' { 'integer','cell' }   []  [];
                        'trialcond'   'integer'   []         []; ...
                        'notrialcond' 'integer'   []         []; ...
                        'sort'        'integer'   []         []; ...
                        'sorttrial'   'string'    { 'on','off' } 'on' }, 'pop_select');
if ischar(g), error(g); end
if ~isempty(g.sort)
    if g.sort, g.sorttrial = 'on';
    else       g.sorttrial = 'off';
    end
end
if strcmpi(g.sorttrial, 'on')
    g.trial = sort(setdiff( g.trial, g.notrial ));
    if isempty(g.trial), error('Error: dataset is empty'); end
else
    g.trial(ismember(g.trial,g.notrial)) = [];
    % still warn about & remove duplicate trials (may be removed in the future)
    [p,q] = unique_bc(g.trial);
    if length(p) ~= length(g.trial)
        disp('Warning: trial selection contained duplicated elements, which were removed.'); 
    end    
    g.trial = g.trial(sort(q));
end

if isempty(g.channel) && ~iscell(g.nochannel) && ~iscell(chanlist)
    g.channel = [1:EEG.nbchan];
end

if iscell(g.channel) && ~iscell(g.nochannel) && ~isempty(EEG.chanlocs)
     noChannelAsCell = {};
     for nochanId = 1:length(g.nochannel)
         noChannelAsCell{nochanId} = EEG.chanlocs(g.nochannel(nochanId)).labels;
     end
     g.nochannel =   noChannelAsCell; 
end

if strcmpi(g.sorttrial, 'on')
    if iscell(g.channel)
         g.channel = sort(setdiff( lower(g.channel), lower(g.nochannel) ));
    else g.channel = sort(setdiff( g.channel, g.nochannel ));
    end
else
    g.channel(ismember(lower(g.channel),lower(g.nochannel))) = [];
    % still warn about & remove duplicate channels (may be removed in the future)
    [p,q] = unique_bc(g.channel);
    if length(p) ~= length(g.channel)
        disp('Warning: channel selection contained duplicated elements, which were removed.'); 
    end    
    g.channel = g.channel(sort(q));    
end

if ~isempty(EEG.chanlocs)
    if strcmpi(g.sorttrial, 'on')
        g.channel = eeg_decodechan(EEG.chanlocs, g.channel);
    else
        % we have to protect the channel order against changes by eeg_decodechan
        if iscell(g.channel)
            % translate channel names into indices
            [inds,names] = eeg_decodechan(EEG.chanlocs, g.channel);
            % and sort the indices back into the original order of channel names
            [tmp,I] = ismember_bc(lower(g.channel),lower(names)); 
            g.channel = inds(I);
        end
    end
end

if ~isempty(g.time) && (g.time(1) < EEG.xmin*1000) && (g.time(2) > EEG.xmax*1000)
   error('Wrong time range');
end
if min(g.trial) < 1 || max( g.trial ) > EEG.trials  
   error('Wrong trial range');
end
if ~isempty(g.channel)
    if min(double(g.channel)) < 1 || max(double(g.channel)) > EEG.nbchan  
        error('Wrong channel range');
    end
end

if size(g.point,2) > 2, 
    g.point = [g.point(1) g.point(end)];
    disp('Warning: vector format for point range is deprecated');
end
if size(g.nopoint,2) > 2, 
    g.nopoint = [g.nopoint(1) g.nopoint(end)];
    disp('Warning: vector format for point range is deprecated');
end
if ~isempty( g.point )
    g.time = zeros(size(g.point));
    for index = 1:length(g.point(:))
        g.time(index) = eeg_point2lat(g.point(index), 1, EEG.srate, [EEG.xmin EEG.xmax]);
    end
    g.notime = [];
end
if ~isempty( g.nopoint )
    g.notime = zeros(size(g.nopoint));
    for index = 1:length(g.nopoint(:))
        g.notime(index) = eeg_point2lat(g.nopoint(index), 1, EEG.srate, [EEG.xmin EEG.xmax]);
    end
    g.time = [];
end
if ~isempty( g.notime )
    if size(g.notime,2) ~= 2
        error('Time/point range must contain 2 columns exactly');
    end
    if g.notime(2) == EEG.xmax
        g.time = [EEG.xmin g.notime(1)];
    else
        if g.notime(1) == EEG.xmin
            g.time = [g.notime(2) EEG.xmax];
        elseif EEG.trials > 1
            error('Wrong notime range. Remember that it is not possible to remove a slice of time for data epochs.');
        end
    end
    if g.notime(end) > EEG.xmax, g.notime(end) = EEG.xmax; end
    if g.notime(1)   < EEG.xmin, g.notime(1)   = EEG.xmin; end
    if floor(max(g.notime(:))) > EEG.xmax 
        error('Time/point range exceed upper data limits');
    end
    if min(g.notime(:)) < EEG.xmin
        error('Time/point range exceed lower data limits');
    end
end
if ~isempty(g.time)
    if size(g.time,2) ~= 2
        error('Time/point range must contain 2 columns exactly');
    end
    for index = 1:length(g.time)
        if g.time(index) > EEG.xmax
            g.time(index) = EEG.xmax;
            disp('Upper time limits exceed data, corrected');
        elseif g.time(index) < EEG.xmin
            g.time(index) = EEG.xmin;
            disp('Lower time limits exceed data, corrected');
        end
    end
end

% select trial values
%--------------------
if ~isempty(g.trialcond)
   try, tt = struct( g.trialcond{:} ); catch
      error('Trial conditions format error');
   end
   ttfields = fieldnames (tt);
   for index = 1:length(ttfields)
        if ~isfield( EEG.epoch, ttfields{index} )
            error([ ttfields{index} 'is not a field of EEG.epoch' ]);
        end;    
        tmpepoch = EEG.epoch;
	    eval( [ 'Itriallow  = find( [ tmpepoch(:).' ttfields{index} ' ] >= tt.' ttfields{index} '(1) );' ] );
	    eval( [ 'Itrialhigh = find( [ tmpepoch(:).' ttfields{index} ' ] <= tt.' ttfields{index} '(end) );' ] );
	    Itrialtmp = intersect_bc(Itriallow, Itrialhigh);
	    g.trial = intersect_bc( g.trial(:)', Itrialtmp(:)');
   end;	   
end

if isempty(g.trial)
   error('Empty dataset, no trial');
end
if length(g.trial) ~= EEG.trials
	fprintf('Removing %d trial(s)...\n', EEG.trials - length(g.trial));
end
if length(g.channel) ~= EEG.nbchan
	fprintf('Removing %d channel(s)...\n', EEG.nbchan - length(g.channel));
end

try
    % For AMICA probabilities...
    %-----------------------------------------------------
    if isfield(EEG.etc, 'amica') && ~isempty(EEG.etc.amica) && isfield(EEG.etc.amica, 'v_smooth') && ~isempty(EEG.etc.amica.v_smooth) && ~isfield(EEG.etc.amica,'prob_added')
        if isfield(EEG.etc.amica, 'num_models') && ~isempty(EEG.etc.amica.num_models)
            if size(EEG.data,2) == size(EEG.etc.amica.v_smooth,2) && size(EEG.data,3) == size(EEG.etc.amica.v_smooth,3) && size(EEG.etc.amica.v_smooth,1) == EEG.etc.amica.num_models

                EEG = eeg_formatamica(EEG);

                %-------------------------------------------

                [EEG com] = pop_select(EEG,args{:});

                %-------------------------------------------

                EEG = eeg_reformatamica(EEG);
                EEG = eeg_checkamica(EEG);
                return;
            else
                disp('AMICA probabilities not compatible with size of data, probabilities cannot be rejected')

                disp('Resuming rejection...')
            end
        end

    end
    % ------------------------------------------------------
catch
    warnmsg = strcat('your dataset contains amica information, but the amica plugin is not installed.  Continuing and ignoring amica information.');
    warning(warnmsg)
end


% recompute latency and epoch number for events
% ---------------------------------------------
if length(g.trial) ~= EEG.trials && ~isempty(EEG.event)
    if ~isfield(EEG.event, 'epoch')
        disp('Pop_epoch warning: bad event format with epoch dataset, removing events');
        EEG.event = [];
    else
        if isfield(EEG.event, 'epoch')
            keepevent = [];
            for indexevent = 1:length(EEG.event)
                newindex = find( EEG.event(indexevent).epoch == g.trial );% For AMICA probabilities...
                %-----------------------------------------------------
                try
                    if isfield(EEG.etc, 'amica') && ~isempty(EEG.etc.amica) && isfield(EEG.etc.amica, 'v_smooth') && ~isempty(EEG.etc.amica.v_smooth) && ~isfield(EEG.etc.amica,'prob_added')
                        if isfield(EEG.etc.amica, 'num_models') && ~isempty(EEG.etc.amica.num_models)
                            if size(EEG.data,2) == size(EEG.etc.amica.v_smooth,2) && size(EEG.data,3) == size(EEG.etc.amica.v_smooth,3) && size(EEG.etc.amica.v_smooth,1) == EEG.etc.amica.num_models
                                
                                EEG = eeg_formatamica(EEG);
                                
                                %-------------------------------------------
                                
                                [EEG com] = pop_select(EEG,args{:});
                                
                                %-------------------------------------------
                                
                                EEG = eeg_reformatamica(EEG);
                                EEG = eeg_checkamica(EEG);
                                return;
                            else
                                disp('AMICA probabilities not compatible with size of data, probabilities cannot be rejected')
                                
                                disp('Resuming rejection...')
                            end
                        end
                        
                    end
                catch
                    warnmsg = strcat('your dataset contains amica information, but the amica plugin is not installed.  Continuing and ignoring amica information.');
                    warning(warnmsg)
                end;                
                % ------------------------------------------------------
                
                if ~isempty(newindex)
                    keepevent = [keepevent indexevent];
                    if isfield(EEG.event, 'latency')
                        EEG.event(indexevent).latency = EEG.event(indexevent).latency - (EEG.event(indexevent).epoch-1)*EEG.pnts + (newindex-1)*EEG.pnts;
                    end
                    EEG.event(indexevent).epoch = newindex;
                end
            end
            diffevent = setdiff_bc([1:length(EEG.event)], keepevent);
            if ~isempty(diffevent)
                disp(['Pop_select: removing ' int2str(length(diffevent)) ' unreferenced events']);
                EEG.event(diffevent) = [];
            end
        end
    end
end


% performing removal
% ------------------
if ~isempty(g.time) || ~isempty(g.notime)
    if EEG.trials > 1
        % select new time window
        % ----------------------    
        try,   tmpevent = EEG.event;
               tmpeventlatency = [ tmpevent.latency ];
        catch, tmpeventlatency = [];
        end
        alllatencies = 1-(EEG.xmin*EEG.srate); % time 0 point
        alllatencies = linspace( alllatencies, EEG.pnts*(EEG.trials-1)+alllatencies, EEG.trials);
        [EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, ...
                                                     [g.time(1) g.time(2)]*EEG.srate, 'allevents', tmpeventlatency);
        tmptime = tmptime/EEG.srate;
        if g.time(1) ~= tmptime(1) && g.time(2)-1/EEG.srate ~= tmptime(2)
            fprintf('pop_select(): time limits have been adjusted to [%3.3f %3.3f] to fit data points limits\n', tmptime(1), tmptime(2)+1/EEG.srate);
        end
        EEG.xmin = tmptime(1);
        EEG.xmax = tmptime(2);
        EEG.pnts = size(EEG.data,2);
        alllatencies = alllatencies(indices);
        
        % modify the event structure accordingly (latencies and add epoch field)
        % ----------------------------------------------------------------------
        allevents = [];
        newevent = [];
        count = 1;
        if ~isempty(epochevent)
            newevent = EEG.event(1);
            for index=1:EEG.trials
                for indexevent = epochevent{index}
                    newevent(count)         = EEG.event(indexevent);
                    newevent(count).epoch   = index;
                    newevent(count).latency = newevent(count).latency - alllatencies(index) - tmptime(1)*EEG.srate + 1 + EEG.pnts*(index-1);
                    count = count + 1;
                end
            end
        end
        EEG.event = newevent;
        
        % erase event-related fields from the epochs
        % ------------------------------------------
        if ~isempty(EEG.epoch)
            fn = fieldnames(EEG.epoch);
            EEG.epoch = rmfield(EEG.epoch,{fn{strmatch('event',fn)}});
        end
    else
        if isempty(g.notime)
            if length(g.time) == 2 && EEG.xmin < 0
                disp('Warning: negative minimum time; unchanged to ensure correct latency of initial boundary event');
            end
            g.notime = g.time';
            g.notime = g.notime(:);
            if g.notime(1) ~= 0, g.notime = [EEG.xmin g.notime(:)'];
            else                 g.notime = [g.notime(2:end)'];
            end
            if g.time(end) == EEG.xmax, g.notime(end) = [];
            else                        g.notime(end+1) = EEG.xmax;
            end
            
            for index = 1:length(g.notime)
                if g.notime(index) ~= 0  && g.notime(index) ~= EEG.xmax
                    if mod(index,2), g.notime(index) = g.notime(index) + 1/EEG.srate;
                    else             g.notime(index) = g.notime(index) - 1/EEG.srate;
                    end
                end
            end;        
            g.notime = reshape(g.notime, 2, length(g.notime)/2)';
        end;   
        
        nbtimes = length(g.notime(:));
        [points,flag] = eeg_lat2point(g.notime(:)', ones(1,nbtimes), EEG.srate, [EEG.xmin EEG.xmax]);
        points = reshape(points, size(g.notime));
        
        % fixing if last region is the same
        if flag
            if ~isempty(find((points(end,1)-points(end,2))== 0)), points(end,:) = []; end
        end
        
        EEG = eeg_eegrej(EEG, points);
    end
end

% performing removal
% ------------------
if ~isequal(g.channel,1:size(EEG.data,1)) || ~isequal(g.trial,1:size(EEG.data,3))
    %EEG.data  = EEG.data(g.channel, :, g.trial);
    % this code belows is prefered for memory mapped files
    diff1 = setdiff_bc([1:size(EEG.data,1)], g.channel);
    diff2 = setdiff_bc([1:size(EEG.data,3)], g.trial);
    if ~isempty(diff1)
         EEG.data(diff1, :, :) = [];
    end
    if ~isempty(diff2)
         EEG.data(:, :, diff2) = [];
    end
end
if ~isempty(EEG.icaact), EEG.icaact = EEG.icaact(:,:,g.trial); end
EEG.trials    = length(g.trial);
EEG.pnts      = size(EEG.data,2);
EEG.nbchan    = length(g.channel);
if ~isempty(EEG.chanlocs)
    EEG.chanlocs = EEG.chanlocs(g.channel);
    if ~isfield(EEG.chaninfo, 'removedchans')
        EEG.chaninfo.removedchans = [];
    end
    try EEG.chaninfo.removedchans = [ EEG.chaninfo.removedchans EEG.chanlocs(diff1) ]; catch, end;
end
if ~isempty(EEG.epoch)
   EEG.epoch = EEG.epoch( g.trial );
end
if ~isempty(EEG.specdata)
	if length(g.point) == EEG.pnts
   		EEG.specdata = EEG.specdata(g.channel, :, g.trial);
   	else
   		EEG.specdata = [];
   		fprintf('Warning: spectral data were removed because of the change in the numner of points\n');
   	end;		
end

% ica specific
% ------------
if ~isempty(EEG.icachansind)
    
    rmchans = setdiff_bc( EEG.icachansind, g.channel ); % channels to remove
    
    % channel sub-indices
    % -------------------
    icachans = 1:length(EEG.icachansind);
    for index = length(rmchans):-1:1
        chanind           = find(EEG.icachansind == rmchans(index));
        icachans(chanind) = [];
    end
        
    % new channels indices
    % --------------------
    count   = 1;
    newinds = [];
    for index = 1:length(g.channel)
        if any(EEG.icachansind == g.channel(index))
            newinds(count) = index;
            count          = count+1;
        end
    end
    EEG.icachansind = newinds;
    
else
    icachans = 1:size(EEG.icasphere,2);
end

if ~isempty(EEG.icawinv)
    flag_rmchan = (length(icachans) ~= size(EEG.icawinv,1));
    if  isempty(EEG.icaweights) || flag_rmchan
        EEG.icawinv    = EEG.icawinv(icachans,:);
        EEG.icaweights = pinv(EEG.icawinv);
        EEG.icasphere  = eye(size(EEG.icaweights,2));
    end
end
if ~isempty(EEG.specicaact)
    if length(g.point) == EEG.pnts
        EEG.specicaact = EEG.specicaact(icachans, :, g.trial);
    else
        EEG.specicaact = [];
        fprintf('Warning: spectral ICA data were removed because of the change in the numner of points\n');
    end
end

% check if only one epoch
% -----------------------
if EEG.trials == 1
    if isfield(EEG.event, 'epoch')
        EEG.event = rmfield(EEG.event, 'epoch');
    end
    EEG.epoch = [];
end
if isfield(EEG.reject, 'gcompreject') && isequal(g.channel,1:size(EEG.data,1))
    tmpgcompreject = EEG.reject.gcompreject;
    EEG.reject = [];
    EEG.reject.gcompreject = tmpgcompreject;
else
    EEG.reject = [];
end
EEG.stats  = [];
EEG.reject.rejmanual = [];
% for stats, can adapt remove the selected trials and electrodes
% in the future to gain time -----------------------------------  
EEG.stats.jp = [];
EEG = eeg_checkset(EEG, 'eventconsistency');

% generate command
% ----------------
if nargout > 1
    com = sprintf('EEG = pop_select( EEG, %s);', vararg2str(args));
end

return;

% ********* OLD, do not remove any event any more
% ********* in the future maybe do a pack event to remove events not in the time range of any epoch

if ~isempty(EEG.event)
    % go to array format if necessary
    if isstruct(EEG.event), format = 'struct';
    else                     format = 'array';
    end
    switch format, case 'struct', EEG = eventsformat(EEG, 'array'); end
    
    % keep only events related to the selected trials
    Indexes = [];
    Ievent  = [];
    for index = 1:length( g.trial )
        currentevents = find( EEG.event(:,2) == g.trial(index));
        Indexes = [ Indexes ones(1, length(currentevents))*index ];
        Ievent  = union_bc( Ievent, currentevents );
    end
    EEG.event = EEG.event( Ievent,: );
    EEG.event(:,2) = Indexes(:);
    
    switch format, case 'struct', EEG = eventsformat(EEG, 'struct'); end
end

