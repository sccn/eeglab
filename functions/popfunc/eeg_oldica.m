% eeg_oldica()  - report, return or add to oldicaweights and oldicasphere 
%                 stored in cell arrays in EEG.etc of an EEGLAB dataset
% Usage:
%        >> eeg_oldica(EEG); % report number of stored oldicaweights 
%        >> [EEG,icaweights, icasphere] = eeg_oldica(EEG,N); % return matrices
%        >> EEG = eeg_oldica(EEG,N,icaweights,icasphere); % add wts and sph 
%                    % matrices to EEG.etc.icaweights and EEG.etc.icasphere
% Inputs:
%        EEG        - EEGLAB dataset structure
%        nreturn    - index of the oldicaweights and sphere to return {default: 1}
%        icaweights - ICA weights matrix to store in EEG.etc.oldicaweights
%        icasphere  - ICA sphere matrix to store in EEG.etc.oldicasphere
% Outputs:
%        icaweights - ICA unmixing matrix (e.g., EEG.icaweights)
%        icasphere  - ICA data sphering matrix (e.g., EEG.icasphere)
%
% See also:   pop_runica()
%
% Author: Scott Makeig, SCCN/INC/UCSD, March 17, 2005

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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

function [EEG,oldicaweights,oldicasphere] = eeg_oldica(EEG, nreturn,icaweights,icasphere)

if nargin< 1
   help eeg_oldica
   return
end
if ~isstruct(EEG)
   error('EEG argument must be a dataset structure')
end

if ~isfield(EEG,'etc')
   error('EEG.etc field not found - no old ICA weights in dataset');
end

if ~isfield(EEG.etc,'oldicaweights') & nargin < 3
   error('EEG.etc.oldicaweights field not found - no old ICA weights in dataset');
elseif nargout < 2 & nargin < 3
   if length(EEG.etc.oldicaweights) > 1
      fprintf('EEG.etc.oldicaweights contains %d weight matrices\n',length(EEG.etc.oldicaweights));
      EEG.etc.oldicaweights
   else
      fprintf('EEG.etc.oldicaweights contains %d weight matrix\n',length(EEG.etc.oldicaweights));
      EEG.etc.oldicaweights
   end
end

if ~isfield(EEG.etc,'oldicasphere') & nargin < 3
   fprintf('EEG.etc.oldicasphere field not found - no old ICA sphere matrix in dataset');
elseif nargout < 2 & nargin < 3
  fprintf('\n');
  if length(EEG.etc.oldicasphere) > 1
   fprintf('EEG.etc.oldicasphere  contains %d weight matrices\n',length(EEG.etc.oldicasphere));
   EEG.etc.oldicasphere
  else
   fprintf('EEG.etc.oldicasphere  contains %d weight matrix\n',length(EEG.etc.oldicasphere));
   EEG.etc.oldicasphere
  end
  fprintf('\n');
end

if nargin < 2
  nreturn = 1;
elseif nargin < 3
  if nreturn< 1 
       error('nreturn must be an oldicaweights index');
  elseif length(EEG.etc.oldicaweights) < nreturn
       fprintf('nreturn (%d) > number of stored oldicaweights (%d) ', nreturn,length(EEG.etc.oldicaweights));
       error(' ');
  end
end

if nargin > 4
   error('too many arguments.\n');
end

if nargin > 2 
  if nargout > 0
    fprintf('New EEG.etc.oldicaweights: ');
    EEG.etc.oldicaweights = [EEG.etc.oldicaweights {icaweights}];
    EEG.etc.oldicaweights 
  else
    error('To update oldicaweights, at least one output required');
  end 
end
if nargin > 3
    fprintf('New EEG.etc.oldicasphere: ');
    EEG.etc.oldicasphere = [EEG.etc.oldicasphere {icasphere}];
    EEG.etc.oldicasphere 
end

if nargout > 1
    % fprintf('\n');
    oldicaweights = EEG.etc.oldicaweights{nreturn};
    % fprintf('Stored oldicaweights matrix (%d) returned (%d,%d).\n',nreturn,size(oldicaweights,1),size(oldicaweights,2));

    if length(EEG.etc.oldicasphere) >= nreturn
      oldicasphere = EEG.etc.oldicasphere{nreturn};
     %  fprintf('Stored oldicasphere  matrix (%d) returned (%d,%d).\n',nreturn,size(oldicasphere,1),size(oldicasphere,2));
    else
      oldicasphere = [];
      % fprintf('No corresponding oldicasphere matrix; [] returned.\n');
   end
   % fprintf('\n');
end
return
