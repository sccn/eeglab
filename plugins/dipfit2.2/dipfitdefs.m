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

% $Log: not supported by cvs2svn $
% Revision 1.3  2005/03/11 18:14:09  arno
% case sensitive problem
%
% Revision 1.2  2005/03/10 18:55:49  arno
% add template files
%
% Revision 1.1  2005/03/10 18:10:27  arno
% Initial revision
%
% Revision 1.16  2003/10/29 16:41:57  arno
% default grid
%
% Revision 1.15  2003/10/29 03:42:55  arno
% same
%
% Revision 1.14  2003/10/29 03:41:30  arno
% meanradius
%
% Revision 1.13  2003/10/29 03:35:20  arno
% remove elc computation
%
% Revision 1.12  2003/09/02 13:01:47  roberto
% added default constraint for symmetry
%
% Revision 1.11  2003/08/01 13:49:49  roberto
% removed 1 and 3 sphere defaults, renamed vol4besa to defaultvolume and added origin
%
% Revision 1.9  2003/06/13 16:48:22  arno
% undo chanlocs checks
%
% Revision 1.8  2003/06/13 01:21:19  arno
% still debuging auto conversion
%
% Revision 1.7  2003/06/13 01:01:34  arno
% debug last
%
% Revision 1.6  2003/06/13 01:00:40  arno
% convert polar to carthesian electrode location strcuture
%
% Revision 1.5  2003/03/12 10:32:12  roberto
% added 4-sphere volume model similar to BESA
%
% Revision 1.4  2003/03/06 15:57:56  roberto
% *** empty log message ***
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this file is not a function but a script and is included in the dipfit_XXX functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(EEG, 'chanlocs')
  error('No electrode locations defined');
end

if ~isfield(EEG, 'icawinv')
  error('No ICA components present');
end

nchan = length(EEG.chanlocs);
ncomp = size(EEG.icawinv, 2);

% create one-sphere model
% defaultvolume.r = meanradius;
% defaultvolume.c = 0.33;
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
folder = which('pop_dipfit_settings');
folder = folder(1:end-21);
delim  = folder(end);
template_models = { ...
    { [ folder 'standard_BESA' delim 'standard_BESA.mat' ] ... % model hdmfile for BESA
      'spherical' ...                                          % coordinate 'spherical' or 'MNI'
      [ folder 'standard_BESA' delim 'avg152t1.mat' ] ...      % MRI MNI normalized file
      [ folder 'standard_BESA' delim 'standard-10-5-cap385.sfp' ] } ... % channel location file
                                                                        % associated with model
    { [ folder 'standard_BEM' delim 'standard_vol.mat' ] ...   % same as above for BEM model
      'MNI' ...
      [ folder 'standard_BEM' delim 'standard_mri.mat' ] ...
      [ folder 'standard_BEM' delim 'elec' delim 'standard_1005.elc' ] } 
    { '' 'MNI' '' '' } }; % custom model

% constrain electrode to sphere
% -----------------------------
meanradius = defaultvolume.r(4);

% defaults for GUI pop_dipfit_settings dialog
defaultelectrodes = sprintf('1:%d', nchan);

% these settings determine the symmetry constraint that can be toggled on
% for the second dipole
defaultconstraint = 'y';      % symmetry along x-axis

% defaults for GUI pop_dipfit_batch dialogs
rejectstr    = '40';	% in percent
xgridstr     = sprintf('linspace(-%d,%d,11)', floor(meanradius), floor(meanradius));
ygridstr     = sprintf('linspace(-%d,%d,11)', floor(meanradius), floor(meanradius));
zgridstr     = sprintf('linspace(0,%d,6)', floor(meanradius));

