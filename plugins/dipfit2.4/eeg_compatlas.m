% eeg_compatlas() - look up the brain area for each component. The pairwise
%                   distance between the component and each of the vertices
%                   of the surface atlas is being computed. The area
%                   selected is the one of the point closest to the
%                   component equivalent dipole. 
%
% Usage:
%   >> OUTEEG = eeg_compatlas(INEEG, 'key1', value1);
%
% Inputs:
%   INEEG         - input EEG dataset structure
%
% Optional inputs
%   'atlas'       - {'dk'} Surface atlas to use. Default is 
%
% Outputs:
%   OUTEEG        - new EEG dataset structure with updated dipfit structure
% 
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2018-
% 
% see also: eeglab()

% Copyright (C) 2018 Arnaud Delorme, arno@ucsd.edu
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

function EEG = eeg_compatlas(EEG, varargin)

if nargin < 1
    help eeg_compatlas;
    return
end

if ~isfield(EEG, 'dipfit') || isempty(EEG.dipfit) || ~isfield(EEG.dipfit, 'model') || isempty(EEG.dipfit.model)
    error('You must run dipole localization first');
end

% decode options
% --------------
g = finputcheck(varargin, ...
  { 'atlas'      'string'    {'dk' }     'dk';
    'components' 'integer'   []          [1:size(EEG.icaweights,1)] });
if isstr(g), error(g); end;

% loading hm file
if isdeployed
    hm = load('-mat', fullfile( ctfroot, 'functions', 'resources', 'head_modelColin27_5003_Standard-10-5-Cap339.mat'));
    if ~exist(meshfile)
        error(sprintf('headplot(): deployed mesh file "%s" not found\n','head_modelColin27_5003_Standard-10-5-Cap339.mat'));
    end
else
    p  = fileparts(which('eeglab.m'));
    hm = load('-mat', fullfile( p, 'functions', 'resources', 'head_modelColin27_5003_Standard-10-5-Cap339.mat'));
end

% coord transform to the HM file space
if strcmpi(EEG.dipfit.coordformat, 'MNI')
    tf = traditionaldipfit([0.0000000000 -26.6046230000 -46.0000000000 0.1234625600 0.0000000000 -1.5707963000 1000.0000000000 1000.0000000000 1000.0000000000]);
elseif strcmpi(EEG.dipfit.coordformat, 'spherical')
    tf = traditionaldipfit([-5.658258      1.039259     -42.80596   -0.00981033    0.03362692   0.004391199      860.8199      926.6112       858.162]);
else
    error('Unknown coordinate format')
end
tfinv = pinv(tf); % the transformation is from HM to MNI (we need to invert it)

% scan dipoles
fprintf('Looking up brain area in the Desikan-Killiany Atlas\n');
for iComp = g.components(:)'
    if size(EEG.dipfit.model(iComp).posxyz,1) == 1
        atlascoord = tfinv * [EEG.dipfit.model(iComp).posxyz 1]';
        
        % find close location in Atlas
        distance = sqrt(sum((hm.cortex.vertices-repmat(atlascoord(1:3)', [size(hm.cortex.vertices,1) 1])).^2,2));
        
        % compute distance to each brain area
        [~,selectedPt] = min( distance );
        area = hm.atlas.colorTable(selectedPt);
        if area > 0
            EEG.dipfit.model(iComp).areadk = hm.atlas.label{area};
        else
            EEG.dipfit.model(iComp).areadk = 'no area';
        end
        
        fprintf('Component %d: area %s\n', iComp, EEG.dipfit.model(iComp).areadk);
    else
        if ~isempty(EEG.dipfit.model(iComp).posxyz)
            fprintf('Component %d: cannot find brain area for bilateral dipoles\n', iComp);
        else
            fprintf('Component %d: no location (RV too high)\n', iComp);
        end
    end
end
        
    