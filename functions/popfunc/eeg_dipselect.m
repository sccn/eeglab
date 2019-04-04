% eeg_dipselect()  - select componet dipoles from an EEG dataset with
%                    reisdual variance (rv) less than a selected threshold
%                    and equivalent dipole location inside the brain volume.
% Usage:
%  >> selctedComponents = eeg_dipselect(EEG, rvThreshold, selectionType, depthThreshold)
%
% Inputs:
%       EEG            - EEGLAB dataset structure
%
% Optional Inputs
%       rvThreshold    - residual variance threshold (%). Dipoles with residual variance
%                        less than this value will be selected. {default = 15}
%       selectionType  - criteria for selecting dipoles:
%                         'rv'  = only by residual variance,
%                         'inbrain' = inside brain volume and residual variance.
%                         {default = 'inbrain'}
%
%       depthThreshold - maximum accepted distance outside brain volume (mm) {default = 1}
%
% Outputs:
%   selctedComponents          - vector of selected components
%
% Example:
%  >> selctedComponents = eeg_dipselect(EEG)             % select in-brain dipoles with rv less than 0.15 (default value)
%  >> selctedComponents = eeg_dipselect(EEG, 20,'rv')   % select dipoles with rv less than 0.2
%
% Author: Nima Bigdely Shamlo, Copyright (C) September 2007
% based on an script from Julie Onton: calls ft_sourcedepth() from Fieldtrip
%
% See also: ft_sourcedepth()

% Copyright (C) Nima Bigdely Shamlo, Copyright (C) September 2007
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

function brainComponents = eeg_dipselect(EEG, rvThreshold, selectionType, depthThreshold);


if nargin<2
    rvThreshold =  0.15;
    fprintf('Maximum residual variance for selected dipoles set to %1.2f (default).\n',rvThreshold);
else
    rvThreshold = rvThreshold/100; % change from percent to value
    if rvThreshold>1
        error('Error: residual variance threshold should be less than 1.\n');
    end
end


if nargin<4
    depthThreshold = 1;
end

% find components with low residual variance

for ic = 1:length(EEG.dipfit.model)
    residualvariance(1,ic) =EEG.dipfit.model(ic).rv;
end
compLowResidualvariance = find(residualvariance <rvThreshold);

if isempty(compLowResidualvariance) || ( (nargin>=3) && strcmp(selectionType, 'rv')) % if only rv is requested (not in-brain)
    brainComponents = compLowResidualvariance;
    return;
else
    if ~exist('ft_sourcedepth')
        selectionType = 'rv';
        tmpWarning = warning('backtrace');
        warning backtrace off;
        warning('You need to install the Fieldtrip extension to be able to select "in brain" dipoles');
        warning(tmpWarning);
        brainComponents = compLowResidualvariance;
        return;
    end
    load(EEG.dipfit.hdmfile);
    
    posxyz = [];
    for c = compLowResidualvariance
        posxyz = cat(1,posxyz,EEG.dipfit.model(c).posxyz(1,:));
    end
    
    depth = ft_sourcedepth(posxyz, vol);
    
    brainComponents = compLowResidualvariance(find(depth<=depthThreshold));
end
