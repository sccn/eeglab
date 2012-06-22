% eeg_epoch2continuous() - convert epoched dataset to continuous dataset
%                          with data epochs separated by boundary events.
% Usage:
%           >> EEGOUT = eeg_epoch2continuous(EEGIN);
%
% Inputs:
%   EEGIN  - a loaded epoched EEG dataset structure.
%
% Inputs:
%   EEGOUT - a continuous EEG dataset structure.
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2012, arno@sccn.ucsd.edu
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

function EEG = eeg_epoch2continuous(EEG)

if nargin < 1
    help eeg_epoch2continuous;
    return;
end;

EEG.data = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3));

for index = 1:EEG.trials-1
    EEG.event(end+1).type     = 'boundary';
    EEG.event(end  ).latency  = index*EEG.pnts-0.5;
    EEG.event(end  ).duration = NaN;
end;

EEG.pnts   = size(EEG.data,2);
EEG.trials = 1;

EEG = eeg_checkset(EEG);