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
% based on an script from Julie Onton and sourcedepth() function
% provided by Robert Oostenveld.
%
% See also: sourcedepth()

function brainComponents = eeg_dipselect(EEG, rvThreshold, selectionType, depthThreshold);


if nargin<2
    rvThreshold =  0.15;
    fprintf('Maximum residual variance for selected dipoles set to %1.2f (default).\n',rvThreshold);
else
    rvThreshold = rvThreshold/100; % change from percent to value
    if rvThreshold>1
        error('Error: residual variance threshold should be less than 1.\n');
    end;
end


if nargin<4
    depthThreshold = 1;
end;

% find components with low residual variance

for ic = 1:length(EEG.dipfit.model)
    residualvariance(1,ic) =EEG.dipfit.model(ic).rv;
end;
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
    end;
    load(EEG.dipfit.hdmfile);
    
    posxyz = [];
    for c = compLowResidualvariance
        posxyz = cat(1,posxyz,EEG.dipfit.model(c).posxyz(1,:));
    end;
    
    depth = ft_sourcedepth(posxyz, vol);
    
    brainComponents = compLowResidualvariance(find(depth<=depthThreshold));
end;
