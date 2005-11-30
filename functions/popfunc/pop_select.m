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
%   "Channel range" - [edit box] RETAIN the indicdated vector of data channels 
%                  Command line equivalent: 'channel' (or 'nochannel' - see below)
%   "Channel range" - [checkbox] EXCLUDE the indicated channel range [min max].
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
%   'channel'     - vector of channel indices to retain in the new dataset
%   'nochannel'   - vector of channel indices to exclude from the new dataset
%   'newname'     - name for the new dataset (OUTEEG)
%
% Outputs:
%   OUTEEG        - new EEG dataset structure
%
% Note: This function performs a conjunction (AND) of all its optional inputs.
%       Using negative counterparts of all options, any logical combination is
%       possible.
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 2001; SCCN/INC/UCSD, 2002-
% 
% see also: eeglab()

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
% Revision 1.46  2005/11/30 17:57:38  arno
% fixing channels
%
% Revision 1.45  2005/11/30 00:44:21  arno
% fix
%
% Revision 1.44  2005/11/30 00:42:05  arno
% fixing icachansind
%
% Revision 1.43  2005/11/30 00:30:07  arno
% fixing icachansind
%
% Revision 1.42  2005/09/27 22:08:08  arno
% fix out of range time limit error
%
% Revision 1.41  2005/08/31 18:22:42  scott
% edited help message  -sm
%
% Revision 1.40  2005/08/01 17:02:50  arno
% now calling eeg_eval
%
% Revision 1.39  2005/07/30 01:09:47  arno
% allowing to process multiple datasets
% .,
%
% Revision 1.38  2005/05/24 17:24:24  arno
% remove cell2mat
%
% Revision 1.37  2004/11/10 16:37:18  arno
% new channel selector
%
% Revision 1.36  2004/03/20 01:30:30  arno
% nothing
%
% Revision 1.35  2003/05/10 17:41:38  arno
% fixing scroll dataset time unit
%
% Revision 1.34  2003/03/19 20:37:24  scott
% printf->fprintf
%
% Revision 1.33  2003/03/07 23:23:02  arno
% debuging out of data limits
%
% Revision 1.32  2003/03/05 18:50:20  arno
% adding doc for new button
%
% Revision 1.31  2003/03/05 18:47:28  arno
% [Aadding channel selection window
%
% Revision 1.30  2003/03/03 22:06:28  arno
% debuging single interval selection
%
% Revision 1.29  2003/02/28 00:57:21  arno
% error if time out of range
%
% Revision 1.28  2003/02/28 00:51:30  arno
% milliseconds -> seconds
%
% Revision 1.27  2003/02/28 00:18:49  arno
% backward compatibility
%
% Revision 1.26  2003/02/28 00:16:14  arno
% reprograming data point selection
%
% Revision 1.25  2003/02/27 20:29:15  arno
% same
%
% Revision 1.24  2003/02/27 20:14:10  arno
% debugging last
%
% Revision 1.23  2003/02/27 03:35:41  arno
% graphic interface help
%
% Revision 1.22  2003/02/27 03:16:02  arno
% debugging for continuous dataset (time selection)
%
% Revision 1.21  2003/01/14 00:28:21  arno
% implementing new time selection (using the epoch function)
%
% Revision 1.20  2002/12/06 03:21:56  arno
% correcting typos
%
% Revision 1.19  2002/10/11 01:22:37  arno
% debugging min time when time select
%
% Revision 1.18  2002/08/21 01:57:57  arno
% undo changes
%
% Revision 1.17  2002/08/21 01:57:16  arno
% nothing
%
% Revision 1.16  2002/08/13 00:08:12  arno
% button problem
%
% Revision 1.15  2002/08/08 02:04:28  arno
% changing message order
%
% Revision 1.14  2002/08/08 01:47:00  arno
% editing text
%
% Revision 1.13  2002/07/29 18:01:41  arno
% changing message
%
% Revision 1.12  2002/07/08 21:42:41  arno
% correcting icawinv channel selection bug
%
% Revision 1.11  2002/06/25 01:55:02  arno
% adding parameter check
%
% Revision 1.10  2002/05/01 01:58:50  arno
% removing extra gui parameters
%
% Revision 1.9  2002/04/30 18:27:11  arno
% adding scroll button
%
% Revision 1.8  2002/04/22 23:43:44  arno
% debugging for no latency event structure
%
% Revision 1.7  2002/04/18 19:16:19  arno
% modifying gui
%
% Revision 1.6  2002/04/11 02:06:24  arno
% adding event consistency check
%
% Revision 1.5  2002/04/09 03:03:26  arno
% removing debug message
%
% Revision 1.4  2002/04/09 01:31:41  arno
% debuging event latency recalculation
%
% Revision 1.3  2002/04/05 23:50:44  arno
% typo
%
% Revision 1.2  2002/04/05 23:28:22  arno
% typo corection
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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
end;
if isempty(EEG(1).data)
    disp('Pop_select error: cannot process empty dataset'); return;
end;    
    
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
         { 'style', 'pushbutton', 'string', '...', ...
           'callback', 'set(findobj(gcbf, ''tag'', ''chans''), ''string'', int2str(pop_chansel({EEG(1).chanlocs.labels}, ''withindex'', ''on'')));' }, ...
           { }, { }, { 'Style', 'pushbutton', 'string', 'Scroll dataset', 'enable', fastif(length(EEG)>1, 'off', 'on'), 'callback', ...
                          'eegplot(EEG.data, ''srate'', EEG.srate, ''winlength'', 5, ''limits'', [EEG.xmin EEG.xmax]*1000, ''position'', [100 300 800 500], ''xgrid'', ''off'', ''eloc_file'', EEG.chanlocs);' } {}};
   results = inputgui( geometry, uilist, 'pophelp(''pop_select'');', 'Select data -- pop_select()' );
   if length(results) == 0, return; end;

          ...
   % decode inputs
   % -------------
   args = {};
   if ~isempty( results{1} )
       if ~results{2}, args = { args{:}, 'time', eval( [ '[' results{1} ']' ] ) };
       else            args = { args{:}, 'notime', eval( [ '[' results{1} ']' ] ) }; end;
   end;

   if ~isempty( results{3} )
       if ~results{4}, args = { args{:}, 'point', eval( [ '[' results{3} ']' ] ) };
       else            args = { args{:}, 'nopoint', eval( [ '[' results{3} ']' ] ) }; end;
   end;

   if ~isempty( results{5} )
       if ~results{6}, args = { args{:}, 'trial', eval( [ '[' results{5} ']' ] ) };
       else            args = { args{:}, 'notrial', eval( [ '[' results{5} ']' ] ) }; end;
   end;

   if ~isempty( results{7} )
       if ~results{8}, args = { args{:}, 'channel', eval( [ '[' results{7} ']' ] ) };
       else            args = { args{:}, 'nochannel', eval( [ '[' results{7} ']' ] ) }; end;
   end;

else
    args = varargin;
end;

% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG com ] = eeg_eval( 'pop_select', EEG, 'warning', 'on', 'params', args);
    return;
end;

% create structure
if ~isempty(args)
   try, g = struct(args{:});
   catch, error('Wrong syntax in function arguments'); end;
else
    g = [];
end;

try, g.time;     catch, g.time = []; end;
try, g.notime;   catch, g.notime = []; end;
try, g.trial;    catch, g.trial = [1:EEG.trials]; end;
try, g.notrial;  catch, g.notrial = []; end;
try, g.point;    catch, g.point = []; end; % will be set up later
try, g.nopoint;  catch, g.nopoint = []; end;
try, g.channel;  catch, g.channel = [1:EEG.nbchan]; end;
try, g.nochannel;        catch, g.nochannel = []; end;
try, g.trialcond;        catch, g.trialcond = []; end;
try, g.notrialcond;      catch, g.notrialcond = []; end;

allfields = fieldnames(g);
for index = 1:length(allfields)
        switch allfields{index}
         case { 'time' 'notime' 'trial' 'notrial' 'point' 'nopoint' 'channel' 'nochannel' 'trialcond' 'notrialcond' };
         otherwise disp('pop_select error: unrecognized option'); beep; return;
        end;
end;

g.trial   = setdiff( g.trial, g.notrial );
g.channel = setdiff( g.channel, g.nochannel );

if ~isempty(g.time) & (g.time(1) < EEG.xmin*1000) & (g.time(2) > EEG.xmax*1000)
   error('Wrong time range');
end;
if min(g.trial) < 1 | max( g.trial ) > EEG.trials  
   error('Wrong trial range');
end;
if min(g.channel) < 1 | max( g.channel ) > EEG.nbchan  
   error('Wrong channel range');
end;

if size(g.point,2) > 2, 
    g.point = [g.point(1) g.point(end)];
    disp('Warning: vector format for point range is deprecated');
end;
if size(g.nopoint,2) > 2, 
    g.nopoint = [g.nopoint(1) g.nopoint(end)];
    disp('Warning: vector format for point range is deprecated');
end;
if ~isempty( g.point )
    g.time = zeros(size(g.point));
    for index = 1:length(g.point(:))
        g.time(index) = eeg_point2lat(g.point(index), 1, EEG.srate, [EEG.xmin EEG.xmax]);
    end;
    g.notime = [];
end;
if ~isempty( g.nopoint )
    g.notime = zeros(size(g.nopoint));
    for index = 1:length(g.nopoint(:))
        g.notime(index) = eeg_point2lat(g.nopoint(index), 1, EEG.srate, [EEG.xmin EEG.xmax]);
    end;
    g.time = [];
end;
if ~isempty( g.notime )
    if size(g.notime,2) ~= 2
        error('Time/point range must contain 2 columns exactly');
    end;
    if g.notime(2) == EEG.xmax
        g.time = [EEG.xmin g.notime(1)];
    else
        if g.notime(1) == EEG.xmin
            g.time = [g.notime(2) EEG.xmax];
        elseif EEG.trials > 1
            error('Wrong notime range. Remember that it is not possible to remove a slice of time for data epochs.');
        end;
    end;
    if max(g.notime(:)) > EEG.xmax | min(g.notime(:)) < EEG.xmin
        error('Time/point range out of data limits');
    end;
end;
if ~isempty(g.time)
    if size(g.time,2) ~= 2
        error('Time/point range must contain 2 columns exactly');
    end;
    for index = 1:length(g.time)
        if g.time(index) > EEG.xmax
            g.time(index) = EEG.xmax;
            disp('Upper time limits exceed data, corrected');
        elseif g.time(index) < EEG.xmin
            g.time(index) = EEG.xmin;
            disp('Lower time limits exceed data, corrected');
        end;
    end;
end;

% select trial values
%--------------------
if ~isempty(g.trialcond)
   try, tt = struct( g.trialcond{:} ); catch
      error('Trial conditions format error');
   end;
   ttfields = fieldnames (tt);
   for index = 1:length(ttfields)
        if ~isfield( EEG.epoch, ttfields{index} )
            error([ ttfields{index} 'is not a field of EEG.epoch' ]);
        end;    
	    eval( [ 'Itriallow  = find( [ EEG.epoch(:).' ttfields{index} ' ] >= tt.' ttfields{index} '(1) );' ] );
	    eval( [ 'Itrialhigh = find( [ EEG.epoch(:).' ttfields{index} ' ] <= tt.' ttfields{index} '(end) );' ] );
	    Itrialtmp = intersect(Itriallow, Itrialhigh);
	    g.trial = intersect( g.trial(:)', Itrialtmp(:)');
   end;	   
end;

if isempty(g.trial)
   error('Empty dataset, no trial');
end;
if length(g.trial) ~= EEG.trials
	fprintf('Removing %d trial(s)...\n', EEG.trials - length(g.trial));
end;
if length(g.channel) ~= EEG.nbchan
	fprintf('Removing %d channel(s)...\n', EEG.nbchan - length(g.channel));
end;

% recompute latency and epoch number for events
% ---------------------------------------------
if length(g.trial) ~= EEG.trials & ~isempty(EEG.event)
    if ~isfield(EEG.event, 'epoch')
        disp('Pop_epoch warning: bad event format with epoch dataset, removing events');
        EEG.event = [];
    else
		if isfield(EEG.event, 'epoch')
			keepevent = [];
			for indexevent = 1:length(EEG.event)
				newindex = find( EEG.event(indexevent).epoch == g.trial );
				if ~isempty(newindex)
					keepevent = [keepevent indexevent];
					if isfield(EEG.event, 'latency')
						EEG.event(indexevent).latency = EEG.event(indexevent).latency - (EEG.event(indexevent).epoch-1)*EEG.pnts + (newindex-1)*EEG.pnts;
					end;
					EEG.event(indexevent).epoch = newindex;
				end;                
			end;
			diffevent = setdiff([1:length(EEG.event)], keepevent);
			if ~isempty(diffevent)
				disp(['Pop_select: removing ' int2str(length(diffevent)) ' unreferenced events']);
				EEG.event(diffevent) = [];
			end;    
		end;
    end;        
end;

% performing removal
% ------------------
if ~isempty(g.time) | ~isempty(g.notime)
    if EEG.trials > 1
        % select new time window
        % ----------------------    
        try,   tmpeventlatency = [ EEG.event.latency ];
        catch, tmpeventlatency = [];
        end;
        alllatencies = 1-(EEG.xmin*EEG.srate); % time 0 point
        alllatencies = linspace( alllatencies, EEG.pnts*(EEG.trials-1)+alllatencies, EEG.trials);
        [EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, ...
                                                     [g.time(1) g.time(2)]*EEG.srate, 'allevents', tmpeventlatency);
        tmptime = tmptime/EEG.srate;
        if g.time(1) ~= tmptime(1) & g.time(2)-1/EEG.srate ~= tmptime(2)
            fprintf('pop_select(): time limits have been adjusted to [%3.3f %3.3f] to fit data points limits\n', tmptime(1), tmptime(2)+1/EEG.srate);
        end;
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
                end;
            end;
        end;
        EEG.event = newevent;
        EEG.epoch = [];
    else
        if isempty(g.notime)
            g.time = g.time';
            if g.time(1) ~= 0, g.notime = [0 g.time(1:end)'];
            else               g.notime = [g.time(2:end)'];
            end;
            if g.time(end) == EEG.xmax, g.notime(end) = [];
            else                        g.notime(end+1) = EEG.xmax;
            end;
            
            for index = 1:length(g.notime)
                if g.notime(index) ~= 0  & g.notime(index) ~= EEG.xmax
                    if mod(index,2), g.notime(index) = g.notime(index) + 1/EEG.srate;
                    else             g.notime(index) = g.notime(index) - 1/EEG.srate;
                    end;
                end;
            end;        
            g.notime = reshape(g.notime, 2, length(g.notime)/2)';
        end;   
        
        nbtimes = length(g.notime(:));
        points = eeg_lat2point(g.notime(:)', ones(1,nbtimes), EEG.srate, [EEG.xmin EEG.xmax]);
        points = reshape(points, size(g.notime));
        EEG = eeg_eegrej(EEG, points);
    end
end;

% performing removal
% ------------------
EEG.data      = EEG.data(g.channel, :, g.trial);
EEG.trials    = length(g.trial);
EEG.pnts      = size(EEG.data,2);
EEG.nbchan    = length(g.channel);
if ~isempty(EEG.chanlocs)
    EEG.chanlocs = EEG.chanlocs(g.channel);
end;    
if ~isempty(EEG.epoch)
   EEG.epoch = EEG.epoch( g.trial );
end;
if ~isempty(EEG.specdata)
	if length(g.point) == EEG.pnts
   		EEG.specdata = EEG.specdata(g.channel, :, g.trial);
   	else
   		EEG.specdata = [];
   		fprintf('Warning: spectral data were removed because of the change in the numner of points\n');
   	end;		
end;

% ica specific
% ------------
if ~isempty(EEG.icachansind)
    
    rmchans = setdiff( EEG.icachansind, g.channel ); % channels to remove
    
    % channel sub-indices
    % -------------------
    icachans = 1:length(EEG.icachansind);
    for index = 1:length(rmchans)
        chanind           = find(EEG.icachansind == rmchans(index));
        icachans(chanind) = [];
    end;
        
    % new channels indices
    % --------------------
    count   = 1;
    for index = 1:length(g.channel)
        if any(EEG.icachansind == g.channel(index))
            newinds(count) = index;
            count          = count+1;
        end;
    end;
    EEG.icachansind = newinds;
    
else
    icachans = 1:size(EEG.icasphere,2);
end;

if ~isempty(EEG.icasphere)
   EEG.icasphere = EEG.icasphere(:,icachans);
end;
if ~isempty(EEG.icawinv)
   EEG.icawinv = EEG.icawinv(icachans,:);
end;
if ~isempty(EEG.icaact)
    EEG.icaact = EEG.icaact(icachans,:,:);
end;
if ~isempty(EEG.specicaact)
    if length(g.point) == EEG.pnts
        EEG.specicaact = EEG.specicaact(icachans, :, g.trial);
    else
        EEG.specicaact = [];
        fprintf('Warning: spectral ICA data were removed because of the change in the numner of points\n');
    end;	
end;

EEG = rmfield( EEG, 'reject');
EEG.reject.rejmanual = [];

% for stats, can adapt remove the selected trials and electrodes
% in the future to gain time -----------------------------------  
EEG = rmfield( EEG, 'stats');
EEG.stats.jp = [];
EEG = eeg_checkset(EEG, 'eventconsistency');

% generate command
% ----------------
com = sprintf('EEG = pop_select( %s,%s);', inputname(1), vararg2str(args));

return;

% ********* OLD, do not remove any event any more
% ********* in the future maybe do a pack event to remove events not in the time range of any epoch

if ~isempty(EEG.event)
    % go to array format if necessary
    if isstruct(EEG.event), format = 'struct';
    else                     format = 'array';
    end;
    switch format, case 'struct', EEG = eventsformat(EEG, 'array'); end;
    
    % keep only events related to the selected trials
    Indexes = [];
    Ievent  = [];
    for index = 1:length( g.trial )
        currentevents = find( EEG.event(:,2) == g.trial(index));
        Indexes = [ Indexes ones(1, length(currentevents))*index ];
        Ievent  = union( Ievent, currentevents );
    end;
    EEG.event = EEG.event( Ievent,: );
    EEG.event(:,2) = Indexes(:);
    
    switch format, case 'struct', EEG = eventsformat(EEG, 'struct'); end;
end;

