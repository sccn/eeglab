% eeg_getica() - get ICA component activation. Recompute if necessary.
%
% >> mergelocs = eeg_getica(EEG, comp);
%
% Inputs: 
%     EEG     - EEGLAB dataset structure
%     comp    - component index
%
% Output: 
%     icaact  - ICA component activity
%
% Author: Arnaud Delorme, 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
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

function icaact = eeg_getica(EEG, comp)

  if nargin < 1
    help eeg_getica;
    return;
  end;
  if nargin < 2
    comp = 1:size(EEG.icaweights,1);
  end;
  
  if ~isempty(EEG.icaact)
    icaact = EEG.icaact(comp,:,:);
  else
    disp('Recomputing ICA activations');
    if isempty(EEG.icachansind)
        EEG.icachansind = 1:EEG.nbchan;
        disp('Channels indices are assumed to be in regular order and arranged accordingly');
    end
    icaact = (EEG.icaweights(comp,:)*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts);
    icaact = reshape( icaact, size(icaact,1), EEG.pnts, EEG.trials);
  end;
