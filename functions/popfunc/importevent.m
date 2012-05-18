% importevent() - Import experimental events from data file or Matlab
%                 array into a structure.
%
% Usage: >> eventstruct = importevent( event, oldevent, srate);
%        >> eventstruct = importevent( event, oldevent, srate, 'key1', 'value1', ...);
%
% Input:
%   event     - [ 'filename'|array ] Filename of a text file, or name of
%               Matlab array in the global workspace containing an
%               array of events in the folowing format: The first column of
%               the cell array is the type of the event, the second the latency. 
%               The others are user-defined. The function can read 
%               either numeric or text entries in ascii files.
%   oldevent  - Old event structure. Used for aligning new events with old
%               ones. Enter [] is no such structure exist.
%   srate     - Sampling rate of the underlying data. Event latencies are
%               expressed in terms of latency in sample point with respect
%               to the data.
%
% Optional file or array input:
%  'fields'   - [Cell array] List of the name of each user-defined column, 
%               optionally followed by a description. Ex: { 'type', 'latency' }
%  'skipline' - [Interger] Number of header rows to skip in the text file 
%  'timeunit' - [ latency unit rel. to seconds ]. Default unit is 1 = seconds. 
%               NaN indicates that the latencies are given in time points.
%  'delim'    - [string] String of delimiting characters in the input file. 
%               Default is tab|space.
%
% Optional oldevent input:
%  'append'   - ['yes'|'no'] 'yes' = Append events to the current events in 
%               the EEG dataset {default}: 'no' = Erase the previous events.
%  'indices'  - {integer vector] Vector indicating the indices of the events to
%               modify. 
%  'align'    - [num] Align the first event latency to the latency of existing 
%               event number (num), and check latency consistency. 
%               A value of 0 indicates that the first events of the pre-defined 
%               events and imported events will be aligned. A positive value (num)
%               aligns the first event to the num-th pre-existing event. 
%               A negative value can also be used; then event number (-num) 
%               is aligned to the first pre-existing event.  Default is 0. 
%               (NaN-> no alignment).
%  'optimalign' - ['on'|'off'] Optimize the sampling rate of the new events so 
%               they best align with old events. Default is 'on'.
%
% Outputs:
%   eventstruct - Event structure containing imported events
%
% Example: >> eventstruct = importevent( 'event_values.txt', [], 256, ...
%         'fields', {'type', 'latency','condition' }, ...
%         'append', 'no', 'align', 0, 'timeunit', 1E-3 );
%
%         This loads the ascii file 'event_values.txt' containing 3 columns 
%         (event_type, latency, and condition). Latencies in the file are
%         in ms (1E-3). The first event latency is re-aligned with the 
%         beginning of the dataset ('align', 0). Any previous events 
%         are erased ('append', 'no').
%
% Author: Arnaud Delorme & Scott Makeig, CNL / Salk Institute, 2004
%
% See also: pop_importevent()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 2004, arno@salk.edu
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

%% REVISION HISTORY
%% function definition, input check
function event = importevent(event, oldevent, srate, varargin)

if nargin < 1
   help importevent;
   return;
end;	

I = [];

% remove the event field
% ----------------------
if ~isempty(oldevent), allfields = fieldnames(oldevent);
else                   allfields = {}; end;
    
g = finputcheck( varargin, { 'fields'    'cell'     []         {};
                         'skipline'  'integer'  [0 Inf]    0;
                         'indices'   'integer'  [1 Inf]    [];
                         'append'    'string'   {'yes';'no';'''yes''';'''no''' }         'yes';
                         'timeunit'  'real'     []         1;
                         'event'     { 'cell';'real';'string' }     []    [];
                         'align'     'integer'  []         NaN;
                         'optimalign' 'string'  { 'on';'off' }         'on';
                         'delim'     {'integer';'string'}   []         char([9 32 44])}, 'importevent');
if isstr(g), error(g); end;
if ~isempty(g.indices), g.append = 'yes'; end;
g.delim = char(g.delim);    

% call from pop_importevent
% -------------------------
if ~isempty(g.event)
    event = g.event;
    if ~isstr(event)
        if size(event,2) > size(event,1), event = event'; end;
    end;
end;
g.event = event;

% determine latency for old event alignment
% -----------------------------------------
tmpalign    = g.align;
g.align     = [];
g.align.val = tmpalign;
if ~isnan(g.align.val)
    if isempty(oldevent)
        error('Setevent: no pre-existing event, cannot perform alignment');
    end;    
    if ~isfield(oldevent, 'latency')
        error('Setevent: pre-existing events do not have a latency field for re-alignment');
    end;    
    switch g.append
        % What is the reasoning behind this warning? It seems misleading. TobyF
        % 3/28/2007
        case {'yes' '''yes'''}, disp('Setevent warning: using align, events should not be appended but erased');
    end;
    if g.align.val < 0
        g.align.event = oldevent(1).latency;
    else
        g.align.event = oldevent(g.align.val+1).latency;
    end
    g.align.nbevent = length(oldevent);
    g.oldevents = oldevent;
    g.align.txt = sprintf([ 'Check alignment between pre-existing (old) and loaded event' ...
                          ' latencies:\nOld event latencies (10 first): %s ...\n' ], ...
                          int2str([ oldevent(1:min(10, length(oldevent))).latency ]));
else
    g.oldevents = [];
end;

tmpfields = fieldnames(g);
event = oldevent;
%% scan all the fields of g
% ------------------------
for curfield = tmpfields'
    if ~isempty(event), allfields = fieldnames(event);
    else                    allfields = {}; end;
    switch lower(curfield{1})
        case {'append', 'fields', 'skipline', 'indices', 'timeunit', 'align', 'delim' }, ; % do nothing now
        case 'event', % load an ascii file
            switch g.append 
                case { '''no''' 'no' } % ''no'' for backward compatibility
                      if isstr(g.event) && ~exist(g.event), g.event = evalin('caller', g.event); end;
                      event = load_file_or_array( g.event, g.skipline, g.delim );
                      allfields = g.fields(1:min(length(g.fields), size(event,2)));
                      event = eeg_eventformat(event, 'struct', allfields);
					  % generate ori fields
					  % -------------------
                      if ~isnan(g.timeunit)
                          for index = 1:length(event)
                              event(index).init_index = index;
                              event(index).init_time  = event(index).latency*g.timeunit;
                          end;
                      end;
                      event = recomputelatency( event, 1:length(event), srate, ...
                                                    g.timeunit, g.align, g.oldevents, g.optimalign);
                case { '''yes''' 'yes' }
                      % match existing fields
                      % ---------------------
                      if isstr(g.event) && ~exist(g.event), g.event = evalin('caller', g.event); end;
                      tmparray = load_file_or_array( g.event, g.skipline, g.delim );
                      if isempty(g.indices) g.indices = [1:size(tmparray,1)] + length(event); end;
                      if length(g.indices) ~= size(tmparray,1)
                            error('Set error: number of row in file does not match the number of event given as input'); 
                      end;

                      % add field
                      % ---------
                      g.fields = getnewfields( g.fields, size(tmparray,2)-length(g.fields));
                      
                      % add new values
                      % ---------------------
                      for eventfield = 1:size(tmparray,2)
                          event = setstruct( event, g.fields{eventfield}, g.indices, { tmparray{:,eventfield} });
                      end;      
					  % generate ori fields
					  % -------------------
                      offset = length(event)-size(tmparray,1);
                      for index = 1:size(tmparray,1)
                          event(index+offset).init_index = index;
                          event(index+offset).init_time  = event(index+offset).latency*g.timeunit;
                      end;
                      event = recomputelatency( event, g.indices, srate, g.timeunit, ...
                                                    g.align, g.oldevents, g.optimalign);
            end;
      end;
end;

if isempty(event) % usefull 0xNB empty structure
    event = [];
end;

%% remove the events wit out-of-bound latencies
% --------------------------------------------
if isfield(event, 'latency')
    try 
        res = cellfun('isempty', { event.latency });
        res = find(res);
        if ~isempty(res)
            fprintf( 'importevent warning: %d/%d event(s) have invalid latencies and were removed\n', ...
                     length(res), length(event));
            event( res ) = [];
        end;
    end;
end;

%% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline, delim );
    if isstr(varname) & exist(varname) == 2  % mean that it is a filename
                                             % --------------------------
        array = loadtxt( varname, 'skipline', skipline, 'delim', delim, 'blankcell', 'off' );
        
    else 
         if ~iscell(varname)
             array = mattocell(varname, ones(1, size(array,1)), ones(1, size(array,2)));
         else
             array = varname;
         end;
    end;     
return;

%% update latency values
% ---------------------
function event = recomputelatency( event, indices, srate, timeunit, align, oldevents, optimalign);

    % update time unit 
    % ----------------
    if ~isfield(event, 'latency'), 
        if isfield(event, 'duration')
            error('A duration field cannot be defined if a latency field has not been defined');
        end;
        return; 
    end;
    if ~isnan(timeunit)
        for index = indices
            event(index).latency  = event(index).latency*srate*timeunit;
            if isfield(event, 'duration')
                event(index).duration = event(index).duration*srate*timeunit;
            end;
        end;
    end;

    % alignment with old events
    % -------------------------
    if ~isnan( align.val )
        if align.val >= 0, alignlatency = event(1).latency;
        else               alignlatency = event(-align.val+1).latency;
        end;
        for index = indices
             event(index).latency = event(index).latency-alignlatency+align.event;
        end;
        if length(event) ~= align.nbevent
            disp('Setevent warning: the number of pre-existing events do not correspond to the number of events');
            disp('                  that were read, so their latencies may have been wrongly re-aligned');
        end;           
        fprintf(align.txt);
        fprintf('New event latencies (10 first): %s ...\n', int2str(round([ event(1:min(10, length(event))).latency ])));
    end;
    if strcmpi(optimalign, 'on') & ~isempty(oldevents)
        newlat = [ event.latency     ];
        oldlat = [ oldevents.latency ];
       
        newlat = repmat(newlat, [length(oldlat) 1]);
        oldlat = repmat(oldlat', [1 size(newlat,2)]);
        if align.val >= 0
            newlat = newlat-newlat(1);
            oldlat = oldlat-oldlat(1+align.val);
        else
            newlat = newlat-newlat(1-align.val);
            oldlat = oldlat-oldlat(1);
        end;
        
        try
            newfactor = fminsearch('eventalign',1,[],newlat, oldlat);
        catch 
            newfactor = fminsearch('eventalign',1,[],[],newlat, oldlat); % Octave
        end;
        fprintf('Best sampling rate ratio found is %1.7f. Below latencies after adjustment\n', newfactor);
        if newfactor > 1.01 | newfactor < 0.99
            disp('Difference is more than 1%, something is wrong; ignoring ratio');
            newfactor = 1;
        else
            difference1 = eventalign(1        ,newlat,oldlat);
            difference2 = eventalign(newfactor,newlat,oldlat);
            fprintf('The average difference before correction was %f sample points\n', difference1);
            fprintf('The average difference after correction is %f sample points\n', difference2);
        end;
        
        %diffarray = abs(newfactor*newlat-oldlat)';
        %[allmins poss] = min(diffarray);
        %figure; hist(allmins);
    else
        newfactor = 1;
    end;
    if ~isnan( align.val ) & newfactor ~= 1 
        if align.val >= 0
            latfirstevent = event(1).latency;
        else
            latfirstevent = event(-align.val+1).latency;
        end;
        for index = setdiff(indices, 1)
            event(index).latency = round(event(index).latency-latfirstevent)*newfactor+latfirstevent;
        end;
        if ~isempty(oldevents)
            fprintf('Old event latencies (10 first): %s ...\n', int2str(round([ oldevents(1:min(10, length(oldevents))).latency ])));
            fprintf('New event latencies (10 first): %s ...\n', int2str(round([ event(1:min(10, length(event))).latency ])));
        end;
    else
        % must add one (because first sample point has latency 0
        % ------------------------------------------------------
        if ~isnan(timeunit)
            for index = indices
                event(index).latency = round((event(index).latency+1)*1000*newfactor)/1000;
            end;
        end;
    end;        

         
%% create new field names
% ----------------------
function epochfield = getnewfields( epochfield, nbfields )
   count = 1;
   while nbfields > 0
       if isempty( strmatch([ 'var' int2str(count) ], epochfield ) )
               epochfield =  { epochfield{:} [ 'var' int2str(count) ] };
               nbfields = nbfields-1;
       else    count = count+1;
       end;                    
   end     
return;

%%
% ----------------------
function var = setstruct( var, fieldname, indices, values )
    if exist('indices') ~= 1, indices = 1:length(var); end
    if ~isempty(values)
        for index = 1:length(indices)
            var = setfield(var, {indices(index)}, fieldname, values{index});
        end
    else
        for index = 1:length(indices)
            var = setfield(var, {indices(index)}, fieldname, '');
        end
    end
return;
