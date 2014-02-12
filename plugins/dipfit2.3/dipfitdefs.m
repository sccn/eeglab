% dipfitdefs() - default settings and filenames for dipolefitting 
%                to source in the ICA/ERP package functions.
%                Insert local dir reference below. 
%
% Note: Edit this file to change local directories under Unix and Windows 
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@miba.auc.dk
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this file is not a function but a script and is included in the dipfit_XXX functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try,
    if ~isfield(EEG, 'chanlocs')
    error('No electrode locations defined');
    end

    if ~isfield(EEG, 'icawinv')
    error('No ICA components present');
    end
    nchan = length(EEG(1).chanlocs);
    ncomp = size(EEG(1).icawinv, 2);
catch, nchan = 0; end;

% create one-sphere model
% defaultvolume.r = meanradius;
% defaultvolume.c = 0.33http://www.google.com/;
% defaultvolume.o = [0 0 0];

% create three-sphere model
% defaultvolume.r = meanradius * [0.0 0.92 0.88];
% defaultvolume.c = [0.33 0.0042 0.33];
% defaultvolume.o = [0 0 0];

% create four-sphere model that is identical to the default of besa
defaultvolume.r = [85-6-7-1 85-6-7 85-6 85];  % in mm
defaultvolume.c = [0.33 1.00 0.0042 0.33];    % brain/csf/skull/skin
defaultvolume.o = [0 0 0];

% default file locations 
% ----------------------
if ~iseeglabdeployed
    folder = which('pop_dipfit_settings');
    folder = folder(1:end-21);
else
    folder = eeglabexefolder;
end;
try,
    delim  = folder(end);
    template_models(1).name     = 'Spherical Four-Shell (BESA)';
    template_models(1).hdmfile  = fullfile(folder, 'standard_BESA', 'standard_BESA.mat');
    template_models(1).mrifile  = fullfile(folder, 'standard_BESA', 'avg152t1.mat');
    template_models(1).chanfile = fullfile(folder, 'standard_BESA', 'standard-10-5-cap385.elp');
    template_models(1).coordformat = 'spherical';
    template_models(1).coord_transform(1).transform = [ ];
    template_models(1).coord_transform(1).keywords  = { 'standard-10-5-cap385' };
    template_models(1).coord_transform(2).transform = [ 13.4299     0.746361    -0.654923  0.000878113   -0.0818352    0.0023747     0.852832     0.941595      0.85887];
    template_models(1).coord_transform(2).keywords  = { 'standard_1005' };
    template_models(1).coord_transform(3).transform = [ -0.254232 0 -8.4081  0 0.00272526  0  8.59463     -10.9643      10.4963 ];
    template_models(1).coord_transform(3).keywords  = { 'gsn' 'sfp' };
    template_models(1).coord_transform(4).transform = [ 0 0 0 0 0.02 0 85 85 85 ];
    template_models(1).coord_transform(4).keywords  = { 'egi' 'elp' };

    template_models(2).name     = 'Boundary Element Model (MNI)';
    template_models(2).hdmfile  = fullfile(folder, 'standard_BEM', 'standard_vol.mat' );
    template_models(2).mrifile  = fullfile(folder, 'standard_BEM', 'standard_mri.mat' );
    template_models(2).chanfile = fullfile(folder, 'standard_BEM', 'elec', 'standard_1005.elc' );
    template_models(2).coordformat = 'MNI';
    template_models(2).coord_transform(1).transform = [ 0 0 0 0 0 -pi/2  1 1 1];
    template_models(2).coord_transform(1).keywords  = { 'standard_1005' };
    template_models(2).coord_transform(2).transform = [ 0.832146  -15.6287 2.41142 0.0812144 0.000937391 -1.5732 1.17419 1.06011 1.14846];
    template_models(2).coord_transform(2).keywords  = { 'standard-10-5-cap385' };
    template_models(2).coord_transform(3).transform = [ 0.0547605 -17.3653 -8.13178 0.0755019 0.00318357 -1.56963 11.7138 12.7933 12.213 ];
    template_models(2).coord_transform(3).keywords  = { 'gsn' 'sfp' };
    template_models(2).coord_transform(4).transform = [ 0 -15 0 0.08 0 -1.571 102 93 100 ];
    template_models(2).coord_transform(4).keywords  = { 'egi' 'elp' };
    
    template_models(3).name     = 'Spherical Four-Shell (custom conductances - see DIPFIT wiki)';
    template_models(3).hdmfile  = fullfile(folder, 'standard_BESA', 'standard_SCCN.mat');
    template_models(3).mrifile  = fullfile(folder, 'standard_BESA', 'avg152t1.mat');
    template_models(3).chanfile = fullfile(folder, 'standard_BESA', 'standard-10-5-cap385.elp');
    template_models(3).coordformat = 'spherical';
    template_models(3).coord_transform(1).transform = [ ];
    template_models(3).coord_transform(1).keywords  = { 'standard-10-5-cap385' };
    template_models(3).coord_transform(2).transform = [ 13.4299     0.746361    -0.654923  0.000878113   -0.0818352    0.0023747     0.852832     0.941595      0.85887];
    template_models(3).coord_transform(2).keywords  = { 'standard_1005' };
    template_models(3).coord_transform(3).transform = [ -0.254232 0 -8.4081  0 0.00272526  0  8.59463     -10.9643      10.4963 ];
    template_models(3).coord_transform(3).keywords  = { 'gsn' 'sfp' };
    template_models(3).coord_transform(4).transform = [ 0 0 0 0 0.02 0 85 85 85 ];
    template_models(3).coord_transform(4).keywords  = { 'egi' 'elp' };
    
    template_models(4).name     = 'Boundary Element Model (custom conductances - see DIPFIT wiki)';
    template_models(4).hdmfile  = fullfile(folder, 'standard_BEM', 'standard_vol_SCCN.mat' );
    template_models(4).mrifile  = fullfile(folder, 'standard_BEM', 'standard_mri.mat' );
    template_models(4).chanfile = fullfile(folder, 'standard_BEM', 'elec', 'standard_1005.elc' );
    template_models(4).coordformat = 'MNI';
    template_models(4).coord_transform(1).transform = [ 0 0 0 0 0 -pi/2  1 1 1];
    template_models(4).coord_transform(1).keywords  = { 'standard_1005' };
    template_models(4).coord_transform(2).transform = [ 0.832146  -15.6287 2.41142 0.0812144 0.000937391 -1.5732 1.17419 1.06011 1.14846];
    template_models(4).coord_transform(2).keywords  = { 'standard-10-5-cap385' };
    template_models(4).coord_transform(3).transform = [ 0.0547605 -17.3653 -8.13178 0.0755019 0.00318357 -1.56963 11.7138 12.7933 12.213 ];
    template_models(4).coord_transform(3).keywords  = { 'gsn' 'sfp' };
    template_models(4).coord_transform(4).transform = [ 0 -15 0 0.08 0 -1.571 102 93 100 ];
    template_models(4).coord_transform(4).keywords  = { 'egi' 'elp' };
    
catch,
    disp('Warning: problem when setting paths for dipole localization');
end;

template_models(5).name        = 'CTF MEG';
template_models(5).coordformat = 'CTF';
template_models(6).name        = 'Custom model files';
template_models(6).coordformat = 'MNI'; % custom model

% constrain electrode to sphere
% -----------------------------
meanradius = defaultvolume.r(4);

% defaults for GUI pop_dipfit_settings dialog
defaultelectrodes = sprintf('1:%d', nchan);

% these settings determine the symmetry constraint that can be toggled on
% for the second dipole
%defaultconstraint = 'y';      % symmetry along x-axis
% PROBLEM: change with respect to the model used. Now just assume perpendicular to nose

% defaults for GUI pop_dipfit_batch dialogs
rejectstr    = '40';	% in percent
xgridstr     = sprintf('linspace(-%d,%d,11)', floor(meanradius), floor(meanradius));
ygridstr     = sprintf('linspace(-%d,%d,11)', floor(meanradius), floor(meanradius));
zgridstr     = sprintf('linspace(0,%d,6)', floor(meanradius));

% Set DipoleDensity path
DIPOLEDENSITY_STDBEM = fullfile(folder, 'standard_BEM', 'standard_vol.mat');
