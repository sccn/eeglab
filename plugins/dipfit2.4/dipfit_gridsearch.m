% dipfit_gridsearch() - do initial batch-like dipole scan and fit to all 
%                       data components and return a dipole model with a 
%                       single dipole for each component.
% 
% Usage: 
%  >> EEGOUT = dipfit_gridsearch( EEGIN, varargin)
%
% Inputs:
%    ...
%
% Optional inputs:
%   'component' - vector with integers, ICA components to scan
%   'xgrid'     - vector with floats, grid positions along x-axis
%   'ygrid'     - vector with floats, grid positions along y-axis
%   'zgrid'     - vector with floats, grid positions along z-axis
%
% Output:
%    ...
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003, load/save by
%         Arnaud Delorme
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
function [EEGOUT] = dipfit_gridsearch(EEG, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert the optional arguments into a configuration structure that can be
% understood by FIELDTRIPs dipolefitting function 
if nargin>2
  cfg = struct(varargin{:});
else
  help dipfit_gridsearch
  return
end

% specify the FieldTrip DIPOLEFITTING configuration
cfg.model      = 'moving';
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'no';
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
if ~isfield(cfg, 'component')
  % default is to scan all components
  cfg.component = 1:size(comp.topo,2);
end

% for each component scan the whole brain with dipoles using FIELDTRIPs
% dipolefitting function
source = ft_dipolefitting(cfg, comp);

% reformat the output dipole sources into EEGLABs data structure
for i=1:length(cfg.component)
  EEG.dipfit.model(cfg.component(i)).posxyz = source.dip(i).pos;
  EEG.dipfit.model(cfg.component(i)).momxyz = reshape(source.dip(i).mom, 3, length(source.dip(i).mom)/3)';
  EEG.dipfit.model(cfg.component(i)).rv     = source.dip(i).rv;
end

EEGOUT = EEG;
