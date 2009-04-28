function [eeg,avgref] = eeg_avgref(eeg,remove)

% eeg_avgref - calculate and remove the eeg average reference
%
% [eeg,avgref] = eeg_avgref(eeg,remove)
%
% eeg - input electrode potential timeseries, Ntime x Nelec
% remove - 1, calculate and remove the average reference from eeg (default)
%          0, calculate the average reference (do not remove it)
% 
% This function assumes that the electrode locations are evenly distributed
% across the scalp surface and the total surface potential of the scalp is
% zero, so the average reference provides an indication of the global zero
% potential.
%
% References
% Offner, F.F. The EEG as potential mapping: the value of the average
% monopolar reference. Electroencephalography and Clinical Neurophysiology,
% 1950, 2, 215-216.
% Lehmann, D. & Skrandies, W. Spatial analysis of evoked potentials in
% man--a review. Progress in Neurobiology, 1984, 23, 227-250.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $

% Copyright (C) 2005  Darren L. Weber
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  09/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $';
fprintf('EEG_AVGREF [v %s]\n',ver(11:15));

if nargin < 1,
    help eeg_avgref;
    return
end

if ~exist('eeg','var'),
    error('no eeg input');
end
if isempty(eeg),
    error('no eeg input');
end

[Nt,Ne] = size(eeg);

fprintf('...%d time points, %d sensors',Nt,Ne);

if ~exist('remove','var'), remove = 1; end
if isempty(remove), remove = 1; end

% calculate and remove the average reference
fprintf('...calculating average reference\n');
avgref = mean(eeg,2);
if remove,
    fprintf('...removing average reference\n');
    eeg = eeg - repmat(avgref,1,Ne);
end

return
