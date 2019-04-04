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
function difftopo(ALLEEG,eeg1,eeg2,interval,limits,subcomps);

if nargin < 3
   help difftopo
   return
end
if eeg1 < 1 || eeg1 > length(ALLEEG)
   help difftopo
   return
end
if eeg2 < 1 || eeg2 > length(ALLEEG)
   help difftopo
   return
end
if eeg1 == eeg2
   help difftopo
   return
end
if ndims(ALLEEG(eeg1).data) ~=3  || ndims(ALLEEG(eeg2).data) ~=3
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
