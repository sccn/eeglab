% eeg_eventtypes()  - find closest channels in a larger EEGLAB chanlocs structure
%                     to channels in a smaller chanlocs structure
% Usage:
%        >> [types,numbers] = eeg_eventtypes(EEG);
% Inputs:
%        EEG        - EEGLAB dataset structure
%                     produce illustrative plots of the BIG and small locations}
% Outputs:
%        types      - cell array of event type strings
%        numbers    - vector giving the numbers of each event type in the data
%
% Author: Scott Makeig, SCCN/INC/UCSD, April 28, 2004

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Scott Makeig, SCCN/INC/UCSD, smakeig@ucsd.edu
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

function [types,numbers] = eeg_eventtypes(EEG)

if ~isstruct(EEG)
   error('EEG argument must be a dataset structure')
end

if ~exist(EEG.event)
   error('EEG.event structure not found');
end

nevents = length(EEG.event);
alltypes = cell(nevents,1)
for k=1:nevents
   alltypes{k} = EEG.event(k).type;
end
[types i j] = unique(alltypes);

ntypes = length(types);
numbers = zeros(ntypes,1);
for k=1:ntypes
 numbers(k) = length(find(j==k));
end
  
if nargout < 1
  for k=1:ntypes
    fprintf('%s\t%d\n',types{k},numbers(k));
  end
end
