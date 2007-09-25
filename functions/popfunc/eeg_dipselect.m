function brainComponents = eeg_dipselect(EEG, rvThreshold, selectionType, depthThreshold);
% eeg_dipselect()  - select componet dipoles of EEG dataset with
%                    reisdual variance (rv) less than a threshold
%                    and location inside brain volume.
% Usage:  
%  >> selctedComponents = eeg_dipselect(EEG, rvThreshold, selectionType, depthThreshold)
%
% Inputs:
%       EEG            - EEGLAB dataset structure
%
% Optional Inputs
%       rvThreshold    - residual variance threshold. Dipoles will residual variance 
%                        less than this value will be selected. {default = 0.15}
%       selectionType  - criteria for selecting dipoles: 
%                         'rv'  = only by residual variance, 
%                         'inbrain' = inside brain volume and residual variance. 
%                         {default = 'inbrain'}
%
%       depthThreshold - maximum accepted distance outside brain volume {default = 1}
%
% Outputs:
%   selctedComponents          - vector of selected components
%
% Example:
%  >> selctedComponents = eeg_dipselect(EEG)             % select in-brain dipoles with rv less than 0.15 (default value)
%  >> selctedComponents = eeg_dipselect(EEG, 0.2,'rv')   % select dipoles with rv less than 0.2
%
% Author: Nima Bigdely Shamlo, Copyright (C) September 2007 
% based on an script from Julie Onton and sourcedepth() function
% provided by Robert Oostenveld.
%  
% See also: sourcedepth()

if nargin<2
    rvThreshold =  0.15;
    fprintf('Maximum residual variance for selected dipoles set to %1.2f (default).\n',rvThreshold);
end;

if nargin<4
    depthThreshold = 1;
end;

% find components with low residual variance

for ic = 1:length(EEG.dipfit.model)
   residualvariance(1,ic) =EEG.dipfit.model(ic).rv;
end;
compLowResidualvariance = find(residualvariance <rvThreshold);

if (nargin>=3) && strcmp(selectionType, 'rv') % if only rv is requested (not in-brain)
    brainComponents = compLowResidualvariance;
    return;
else
    load(EEG.dipfit.hdmfile);

    posxyz = cell2mat({EEG.dipfit.model(compLowResidualvariance).posxyz}');% select positions for components with low residual variance
    depth = sourcedepth(posxyz, vol);

    brainComponents = compLowResidualvariance(find(depth<=depthThreshold));
end;