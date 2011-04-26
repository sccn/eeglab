% eegthresh() - classical trial rejection rejection using a threshold on 
%             the raw signal
%
% Usage:
%   >> [Iin, Iout, newsignal, elec] = eegthresh( signal, frames, ...
%              elecs, negthresh, posthresh, timerange, starttime, entime)
%
% Inputs:
%   signal     - signal 2D, channels x (frames*sweeps) or 3D channels x 
%              frames x sweeps
%   frames     - number of points per epoch
%   elecs      - [e1 e2 ...] electrodes (number) to take into 
%              consideration for rejection
%   negthresh  - negative threshold limit in mV (can be an array if 
%              several electrodes; if less numbe  of values than number 
%              of electrodes the last value is used for the remaining 
%              electrodes)
%   posthresh  - positive threshold limit in mV (same syntax as 
%              negthresh)
%   timerange  - [mintime maxtime] timerange limit in second 
%   starttime  - starting time limit in second (same syntax  as 
%              negthresh)
%   endtime    - ending time limit in second (same syntax  as negthresh)
%
% Outputs:
%   Iin        - Indexes of sweeps not rejected
%   Iout       - Indexes of sweeps rejected
%   newsignal  - signal after 
%              sweeps rejection (same number of dimension as signal)
%   elec       - electrode that served for the rejection (array of 0 and
%              1, same number of column as Iout, rows=number of 
%              electrodes)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_eegthresh(), eeglab()

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

function [Iin, Iout, newsignal, elec] = eegthresh( signal, pnts, electrodes, negthresh, posthresh, timerange, starttime, endtime);

if nargin < 7
	help eegthresh;
	return;
end;

if starttime < timerange(1)
	disp('eegthresh: starting point out of range, adjusted');
    startime = timerange(1);
end;	
if endtime > timerange(2)
	disp('eegthresh: ending point out of range, adjusted');
    endtime = timerange(2);
end;	

% complete thresholds values if necessary
%----------------------------------------
if size(posthresh,2) < size(electrodes,2)
	posthresh = [ posthresh posthresh(end)*ones(1,size(electrodes,2)-size(posthresh,2))];
end;	
if size(negthresh,2) < size(electrodes,2)
	negthresh = [ negthresh negthresh(end)*ones(1,size(electrodes,2)-size(negthresh,2))];
end;
	
% complete timeselect values if necessary
%----------------------------------------
if size(starttime,2) < size(electrodes,2)
	starttime = [ starttime starttime(end)*ones(1,size(electrodes,2)-size(starttime,2))];
end;	
if size(endtime,2) < size(electrodes,2)
	endtime = [ endtime endtime(end)*ones(1,size(electrodes,2)-size(endtime,2))];
end;	

% find the maximum for each trial
%--------------------------------
sweeps = size(signal(1,:),2)/pnts;
signal = reshape(signal(electrodes,:), size(electrodes(:),1), pnts, sweeps);

% reject the selected trials
%---------------------------
elec = zeros(size(electrodes(:),1), sweeps);
allelec = zeros(1, sweeps);
for indexe = 1:size(electrodes(:),1)
	% transform the time range
	% ------------------------
	framelowlimit  = max(1,floor((starttime(indexe)-timerange(1))/(timerange(2)-timerange(1))*(pnts-1))+1);	
	framehighlimit = floor((endtime(indexe)  -timerange(1))/(timerange(2)-timerange(1))*(pnts-1))+1;

	% remove outliers
	% ---------------
	sigtmp = squeeze(signal(indexe,framelowlimit:framehighlimit,:));
	sigmax = max(sigtmp, [], 1);
	sigmin = min(sigtmp, [], 1);
	elec(indexe,:) = ( sigmin < negthresh(indexe) ) | ( sigmax > posthresh(indexe) );
	allelec = allelec | elec(indexe,:);
end;
Iout = find( allelec == 1 );
Iin  = find( allelec == 0 );
elec = elec(:,Iout);

% reject the selected trials
%---------------------------
newsignal = reshape(signal, size(signal,1), pnts, sweeps);
if ~isempty(Iin)
	newsignal = newsignal(:,:,Iin);
	if ndims(signal) == 2
		newsignal = newsignal(:,:);
	end;
else
	newsignal = [];
end;
return;
