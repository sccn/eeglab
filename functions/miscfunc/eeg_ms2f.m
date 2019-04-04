% eeg_ms2f() - convert epoch latency in ms to nearest epoch frame number
%
% Usage:
%         >> outf = eeg_ms2f(EEG,ms);
% Inputs:
%         EEG - EEGLAB data set structure
%         ms  - epoch latency in milliseconds
% Output:
%         outf - nearest epoch frame to the specified epoch latency
%
% Author: Scott Makeig, SCCN/INC, 12/05

% Copyright (C) Scott Makeig, SCCN/INC, 12/05
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

function outf = eeg_ms2f(EEG,ms)
ms = ms/1000;
if ms < EEG.xmin || ms > EEG.xmax
    error('time out of range');
end
outf = 1+round((EEG.pnts-1)*(ms-EEG.xmin)/(EEG.xmax-EEG.xmin));
% else
% [tmp outf] = min(abs(EEG.times-ms));
