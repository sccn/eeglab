% eeg_lat2point() - convert latencies in time units relative to the
%                   time locking event of an eeglab() data epoch to 
%                   latencies in data points (assuming concatenated epochs).
% Usage:
%       >> [newlat] = eeg_lat2point( lat_array, epoch_array,...
%                                 srate, timelimits, timeunit);
%       >> [newlat] = eeg_lat2point( lat_array, epoch_array,...
%                                 srate, timelimits, '','outrange',1);
% Inputs:
%   lat_array   - latency array in 'timeunit' units (see below)
%   epoch_array - epoch number for each latency
%   srate       - data sampling rate in Hz
%   timelimits  - [min max] epoch timelimits in 'timeunit' units (see below)
%   timeunit    - time unit relative to seconds. Default is 1 = seconds.
%
% Optional inputs:
%   outrange    - [1/0] Replace the points out of the range with the value of
%                 the maximun point in the valid range or raise an error.
%                 Default [1] : Replace point.
%
% Outputs:
%   newlat      - converted latency values in points assuming concatenated
%                 data epochs (see eeglab() event structure)
%   flag        - 1 if any point out of range was replaced.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2 Mai 2002
%
% See also: eeg_point2lat(), eeglab()

% Copyright (C) 2 Mai 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [newlat,flag] = eeg_lat2point( lat_array, epoch_array, srate, timewin, timeunit, varargin);
% -------------------------------------------------------------------------
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end
catch
    error('std_checkdatasession() error: calling convention {''key'', value, ... } error');
end
 try, g.outrange; catch, g.outrange = 1; end; %
 
 flag = 0;
% -------------------------------------------------------------------------

if nargin <4
    help eeg_lat2point;
    return;
end;    
if nargin <5 || isempty(timeunit)
	timeunit = 1;
end

if length(lat_array) ~= length(epoch_array)
	if length(epoch_array)~= 1
		disp('eeg_lat2point: latency and epochs must have the same length'); return;
	else
		epoch_array = ones(1,length(lat_array))*epoch_array;
	end
end
if length(timewin) ~= 2
    disp('eeg_lat2point: timelimits must have length 2'); return;
end
if iscell(epoch_array)
	epoch_array = [ epoch_array{:} ];
end
if iscell(lat_array)
	lat_array = [ lat_array{:} ];
end

timewin = timewin*timeunit;
pnts = (timewin(2)-timewin(1))*srate+1;
newlat  = (lat_array*timeunit-timewin(1))*srate+1 + (epoch_array-1)*pnts;

% Detecting points out of range (RMC)
% Note: This is neccesary since the double precision  multiplication could lead to the
% shifting in one sample out of the valid range

if  and(~isempty(newlat),~isempty(epoch_array)) && max(newlat(:)) > max((epoch_array)*pnts)
    if g.outrange == 1
        IndxOut = find(newlat(:) > max((epoch_array)*pnts));
        newlat(IndxOut) = max((epoch_array)*pnts);
        flag = 1;
        warning('eeg_lat2point(): Points out of range detected. Points replaced with maximum value');
    elseif g.outrange == 0
        error('Error in eeg_lat2point(): Points out of range detected');
    end
end
