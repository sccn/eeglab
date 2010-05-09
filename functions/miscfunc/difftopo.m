% difftopo - compute and plot component decomposition for the difference ERP 
%            between two EEG datasets. Plots into the current axis (gca); 
%            plot into a new empty figure as below.
% Usage:
%          >> figure; difftopo(ALLEEG,eeg1,eeg2,interval);
% Inputs:
%         ALLEEG - array of leaded EEG datasets
%         eeg1   - index of subtrahend (+) dataset
%         eeg2   - index of minuend (-) dataset
%         interval - [minms maxms] latency interval in which to find 
%                  and plot largest contributing components {default: whole epoch}
%         limits - [stms endms] latenc plotting limits (in ms) {default|0: whole epoch}
%         subcomps - array of component indices to remove from data 
%                  (e.g. eye movement components) {default|0: none}
%
% Outputs: none
%
% Author: Scott Makeig, SCCN/INC/UCSD 3/28/05

% Copyright (C) Scott Makeig, SCCN/INC/UCSD 3/28/05
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
function difftopo(ALLEEG,eeg1,eeg2,interval,limits,subcomps);

if nargin < 3
   help difftopo
   return
end
if eeg1 < 1 | eeg1 > length(ALLEEG)
   help difftopo
   return
end
if eeg2 < 1 | eeg2 > length(ALLEEG)
   help difftopo
   return
end
if eeg1 == eeg2
   help difftopo
   return
end
if ndims(ALLEEG(eeg1).data) ~=3  | ndims(ALLEEG(eeg2).data) ~=3
  error('EEG datasets must be epoched data');
end

if nargin < 4
   interval = [ALLEEG(eeg1).xmin*999.9 ALLEEG(eeg1).xmax*999.9];
end
if nargin < 5
   limits = 0; % [ALLEEG(eeg1).xmin*1000 ALLEEG(eeg1).xmax*1000];
end
if nargin < 6
   subcomps = [];
end

if ~isfield(ALLEEG(eeg1),'icaweights')
  error('EEG datasets must have icaweights');
end
if length(interval) ~= 2
   help difftopo
   return
end

if ALLEEG(eeg1).pnts ~= ALLEEG(eeg1).pnts
  error('EEG datasets must have the same number of frames per epoch');
end

set(gcf,'Name','difftopo()');
diff = mean(ALLEEG(eeg1).data,3) - mean(ALLEEG(eeg2).data,3);

plottitle = [ ALLEEG(eeg1).setname ' - ' ALLEEG(eeg2).setname];
envtopo(diff,ALLEEG(eeg1).icaweights*ALLEEG(eeg1).icasphere,...
       'chanlocs',ALLEEG(eeg1).chanlocs, ...
       'timerange',[ALLEEG(eeg1).xmin*1000 ALLEEG(eeg1).xmax*1000],...
       'limits',limits,...
       'limcontrib',interval,...
       'title',plottitle, ...
       'subcomps',subcomps);
