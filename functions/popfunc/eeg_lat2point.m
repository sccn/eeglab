% eeg_lat2point() - convert latencies in time units relative to the
%                   time locking event of an eeglab() data epoch to 
%                   latencies in data points (assuming concatenated epochs).
% Usage:
%       >> [newlat] = eeg_lat2point( lat_array, epoch_array,...
%                                 srate, timelimits, timeunit);
% Inputs:
%   lat_array   - latency array in 'timeunit' units (see below)
%   epoch_array - epoch number for each latency
%   srate       - data sampling rate in Hz
%   timelimits  - [min max] epoch timelimits in 'timeunit' units (see below)
%   timeunit    - time unit relative to seconds. Default is 1 = seconds.
%
% Outputs:
%   newlat      - converted latency values in points assuming concatenated
%                 data epochs (see eeglab() event structure)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2 Mai 2002
%
% See also: eeg_point2lat(), eeglab()

% Copyright (C) 2 Mai 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function newlat = eeg_lat2point( lat_array, epoch_array, srate, timewin, timeunit);

if nargin <4
    help eeg_lat2point;
    return;
end;    
if nargin <5
	timeunit = 1;
end;

if length(lat_array) ~= length(epoch_array)
	if length(epoch_array)~= 1
		disp('eeg_lat2point: latency and epochs must have the same length'); return;
	else
		epoch_array = ones(1,length(lat_array))*epoch_array;
	end;
end;
if length(timewin) ~= 2
    disp('eeg_lat2point: timelimits must have length 2'); return;
end;
if iscell(epoch_array)
	epoch_array = [ epoch_array{:} ];
end;
if iscell(lat_array)
	lat_array = [ lat_array{:} ];
end

timewin = timewin*timeunit;
pnts = (timewin(2)-timewin(1))*srate+1;
newlat  = (lat_array*timeunit-timewin(1))*srate+1 + (epoch_array-1)*pnts;
