% epoch() - Extract epochs time locked to specified events from continuous EEG data.
%
% Usage:
%   >> epocheddata = epoch( data, events, timelim);
%   >> [epocheddata, newtime, indices, rerefevent, rereflatencies ] = ...
%                      epoch( data, events, timelim, 'key1', value1, ... );
%
% Inputs:
%   data       - input data (chan,frames). In the case, data is a 
%                3D array (chan, frames, epochs), epochs are extracted
%                only if their time windows fall within existing 
%                pre-existing epochs.
%   events     - vector events (expressed in samples)
%   timelim    - [init end] in second centered on the events (i.e. [-1 2])
%
% Optional inputs:
%   'srate'      - sampling rate in Hz for events expressed in seconds.
%                Required.
%   'valuelim'   - [min max] data limits. If one positive value is given,
%                the opposite value is used for lower bound. For example, 
%                use [-50 50] to remove artifactual epoch. Default: none.
%   'verbose'    - ['on'|'off']. Default is 'on'.
%   'allevents'  - event vector containing the latencies of all events
%                (not only those used for epoching). The function
%                return an array 'rerefevent' which contain the latency
%                of these events in each trials (assuming the events are
%                included in the trial time window). These events can
%                be in point or second but they must be in the same format
%                as the events used for epoching.
%   'alleventrange' - for event selection, defines a time range [start end] (in 
%                seconds or data points) relative to the time-locking events. 
%                Default is same as 'timelim'.
%
% Outputs:
%   epocheddata - output (chan, frames, epochs)
%   indices     - indices of accepted events
%   newtime     - new time limits. See notes.
%   rerefevent  - re-referenced event cell array (size nbepochs) of array 
%                 indices for each epochs (note that the number of events 
%                 per trial may vary).
%   rereflatencies - re-referenced latencies event cell array (same as above
%                 but indicates event latencies in epochs instead of event 
%                 indices). 
%
% Note: maximum time limit will be reduced by one point with comparison to the
%       input time limits. For instance at 100 Hz, 3 seconds last 300 points, 
%       but if we assign time 0 to the first point, then we must assign time 
%       299 to the last point.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_epoch(), eeglab()

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

function [epochdat, newtime, indexes, alleventout, alllatencyout, reallim] = epoch( data, events, lim, varargin );

if nargin < 2
   help epoch;
	return;
end;	
alleventout = {};

% create structure
% ----------------
if ~isempty(varargin)
   try, g = struct(varargin{:});
   catch, error('Epoch: wrong syntax in function arguments'); end
else
    g = [];
end

try, g.srate; 	 	     catch, g.srate = 1; end
try, g.valuelim; 	     catch, g.valuelim =  [-Inf Inf]; end
try, g.verbose; 	     catch, g.verbose = 'on'; end
try, g.allevents; 	     catch, g.allevents = []; end
try, g.alleventrange; 	 catch, g.alleventrange = lim; end

% computing point limits
% ----------------------
reallim(1) = round(lim(1)*g.srate);   % compute offset
reallim(2) = round(lim(2)*g.srate-1); % compute offset

% epoching
% --------
fprintf('Epoching...\n');
newdatalength = reallim(2)-reallim(1)+1;

eeglab_options;
if option_memmapdata == 1
     epochdat = mmo([], [size(data,1), newdatalength, length(events)]);
else epochdat = zeros( size(data,1), newdatalength, length(events) );
end
g.allevents =  g.allevents(:)';
datawidth  = size(data,2)*size(data,3);
dataframes = size(data,2);
indexes = zeros(length(events),1);
alleventout = {};
alllatencyout = {};
for index = 1:length(events)
	pos0 = floor(events(index)*g.srate); % offset of time locking event
	posinit = pos0+reallim(1); % compute offset
	posend  = pos0+reallim(2); % compute offset
   
   if floor((posinit-1)/dataframes) == floor((posend-1)/dataframes) && posinit >= 1 && posend <= datawidth % test if within boundaries
      tmpdata = data(:,posinit:posend);
      epochdat(:,:,index) = tmpdata;
      if ~isinf(g.valuelim(1)) || ~isinf(g.valuelim(2))
          tmpmin = min(reshape(tmpdata, prod(size(tmpdata)),1));
          tmpmax = max(reshape(tmpdata, prod(size(tmpdata)),1));
          if (tmpmin > g.valuelim(1)) && (tmpmax < g.valuelim(2))
              indexes(index) = 1;
          else
              switch g.verbose, case 'on', fprintf('Warning: event %d out of value limits\n', index); end
          end;   
      else 
          indexes(index) = 1;
      end
   else
      switch g.verbose, case 'on', fprintf('Warning: event %d out of data boundary\n', index); end
   end

   % rereference events
   % ------------------
   if ~isempty(g.allevents)
        posinit = pos0 + g.alleventrange(1)*g.srate; % compute offset
        posend  = pos0 + g.alleventrange(2)*g.srate; % compute offset
        eventtrial = intersect_bc( find(g.allevents*g.srate >= posinit),  find(g.allevents*g.srate <= posend) );
        alleventout{index} = eventtrial;
        alllatencyout{index} = g.allevents(eventtrial)*g.srate-pos0; 
   end
end;   
newtime(1) = reallim(1)/g.srate;
newtime(2) = reallim(2)/g.srate;

epochdat(:,:,find(indexes == 0)) = [];
indexes = find(indexes == 1);
%epochdat = epochdat(:,:,indexes);
if ~isempty(alleventout)
    alleventout = alleventout(indexes);
    alllatencyout= alllatencyout(indexes);
end
reallim = reallim*g.srate;
return;

%% GENERATION OF NAN IN ARRAYS (old implementation)
%% ------------------------------------------------
%alleventout(index,1:length(eventtrial) ) = eventtrial+1000000;
%%then replace all zeros by Nan and subtract 1000000
%if ~isempty(alleventout)
%    alleventout( find( alleventout == 0) ) = nan;
%    alleventout = alleventout(indexes,:) - 1000000;
%   alllatencyout( find( alllatencyout == 0) ) = nan;
%   alllatencyout = alllatencyout(indexes,:) - 1000000;
%end

function res = lat2point( lat, srate, pnts);

res = lat*srate+1 + (epoch_array-1)*pnts;
