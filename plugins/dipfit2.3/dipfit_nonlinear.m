% dipfit_nonlinear() - perform nonlinear dipole fit on one of the components
%                   to improve the initial dipole model. Only selected dipoles
%                   will be fitted.
%
% Usage: 
%  >> EEGOUT = dipfit_nonlinear( EEGIN, optarg)
%
% Inputs:
%    ...
%
% Optional inputs are specified in key/value pairs and can be:
%    ...
%
% Output:
%    ...
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%         Thanks to Nicolas Robitaille for his help on the CTF MEG
%         implementation

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@smi.auc.dk
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
function [EEGOUT] = dipfit_nonlinear( EEG, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert the optional arguments into a configuration structure that can be
% understood by FIELDTRIPs dipolefitting function 
if nargin>2
  cfg = struct(varargin{:});
else
  help dipfit_nonlinear
  return
end

% specify the FieldTrip DIPOLEFITTING configuration
cfg.model      = 'moving';
cfg.gridsearch = 'no';
if ~isfield(cfg, 'nonlinear')
  % if this flag is set to 'no', only the dipole moment will be fitted
  cfg.nonlinear  = 'yes';
end
% add some additional settings from EEGLAB to the configuration
tmpchanlocs    = EEG.chanlocs;
cfg.channel    = { tmpchanlocs(EEG.dipfit.chansel).labels };
if isfield(EEG.dipfit, 'vol')
    cfg.vol        = EEG.dipfit.vol;
elseif isfield(EEG.dipfit, 'hdmfile')
    cfg.hdmfile    = EEG.dipfit.hdmfile;
else
    error('no head model in EEG.dipfit')
end

if isfield(EEG.dipfit, 'elecfile') & ~isempty(EEG.dipfit.elecfile)
    cfg.elecfile = EEG.dipfit.elecfile;
end
if isfield(EEG.dipfit, 'gradfile') & ~isempty(EEG.dipfit.gradfile)
    cfg.gradfile = EEG.dipfit.gradfile;
end

% set up the initial dipole model based on the one in the EEG structure
cfg.dip.pos = EEG.dipfit.model(cfg.component).posxyz;
cfg.dip.mom = EEG.dipfit.model(cfg.component).momxyz';
cfg.dip.mom = cfg.dip.mom(:);

% convert the EEGLAB data structure into a structure that looks as if it
% was computed using FIELDTRIPs componentanalysis function
comp = eeglab2fieldtrip(EEG, 'componentanalysis', 'dipfit');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Added code to handle CTF data with multipleSphere head model           %
%  This code is copy-pasted in dipfit_gridSearch, dipfit_nonlinear        %
%  The flag .isMultiSphere is used by dipplot                             %
%  Nicolas Robitaille, January 2007.                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Do some trick to force fieldtrip to use the multiple sphere model
if strcmpi(EEG.dipfit.coordformat, 'CTF')
   cfg = rmfield(cfg, 'channel');
   comp = rmfield(comp, 'elec');
   cfg.gradfile = EEG.dipfit.chanfile;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit the dipoles to the ICA component(s) of interest using FIELDTRIPs
% dipolefitting function

currentPath = pwd;
ptmp = which('ft_prepare_vol_sens');
ptmp = fileparts(ptmp);
if isempty(ptmp), error('Path to "forward" folder of Fieldtrip missing'); end;
cd(fullfile(ptmp, 'private'));
try,
    source = ft_dipolefitting(cfg, comp);
catch,
    cd(currentPath);
    lasterr
    error(lasterr);
end;
cd(currentPath);

% reformat the output dipole sources into EEGLABs data structure
EEG.dipfit.model(cfg.component).posxyz  = source.dip.pos;
EEG.dipfit.model(cfg.component).momxyz  = reshape(source.dip.mom, 3, length(source.dip.mom)/3)';
EEG.dipfit.model(cfg.component).diffmap = source.Vmodel - source.Vdata;
EEG.dipfit.model(cfg.component).sourcepot = source.Vmodel;
EEG.dipfit.model(cfg.component).datapot   = source.Vdata;
EEG.dipfit.model(cfg.component).rv        = source.dip.rv;
%EEG.dipfit.model(cfg.component).rv = sum((source.Vdata - source.Vmodel).^2) / sum( source.Vdata.^2 );

try 
    EEG = eeg_compatlas(EEG, 'components', cfg.component);
catch
    disp('Fail to look up brain areas');
end
EEGOUT = EEG;
