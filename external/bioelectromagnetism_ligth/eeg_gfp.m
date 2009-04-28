function [gfp,gd] = eeg_gfp(eeg,avg)

% eeg_gfp - calculate Global Field Power
%
% [gfp,gd] = eeg_gfp(eeg,avg)
%
% eeg - input electrode potential timeseries, Ntime x Nelec
% avg - 1, calculate and subtract the average reference (default)
%       0, assume eeg input is already corrected by the average reference
%       This assumes that the electrode locations have not only high
%       spatial frequency but also an even distribution across the scalp
%       surface and also that the total surface potential of the
%       scalp is zero, so the average reference then provides an
%       indication of the global zero potential.
%
% gfp - Global Field Power (reference-free)
% gd  - Global Dissimilarity Index
%
% -------------------------------------------------------------------------
% Notes
%
% Global Field Power (reference-free), Lehmann & Skrandies (1984, p. 235):
% We use numerical procedures which assess the degree of 'hilliness' (or
% 'relief', or electrical strength) of the fields (Lehmann, 1971, 1972;
% Lehmann and Skrandies, 1980). The principal approach is to consider all
% possible potential differences in the field (for n electrodes, n*(n - 1))
% with equal weight, and thus compute the reference-free, mean potential
% difference (global field power) at each moment in time using formula (1)
%
% GFP (reference-free) = sqrt( (1/2n)* sum_i(sum_j(ui - uj)^2) ) (1)
%
% where n is the number of electrodes which measure the potentials 
% ei and ej; i, j = 1 . . . n; the observed voltages are 
% ui = ei - common reference. We note that meaningful computation of
% electric field power is based on data from electrodes spaced about
% equidistantly over a reasonably large scalp recording area.
% -------------------------------------------------------------------------
% Global Dissimilarity Index, Lehmann & Skrandies (1984, p. 240):
% For the comparison of the general configuration of two maps by a single
% number indicator we have proposed the measure of 'shape dissimilarity'
% (Lehmann and Skrandies, 1980a; Skrandies and Lehmann, 1982b). The two
% maps to be compared are first scaled for equal voltage range or global
% field power to omit differences in overall amplitude, since only shape of
% the distribution should be considered. The index of dissimilarity is
% computed as the mean standard deviation (MSD) over all recording points
% of the average voltages between the two maps at each recording point,
% referred to the average reference of the given map, using formula:
%
% dissimilarity index = MSD = (1/2n)* sum_i( sqrt((u(m1_i)-u(m2_i))^2))
%
% -------------------------------------------------------------------------
% References
% Lehmann, D. & Skrandies, W. Reference-free identification of components
% of checkerboard-evoked multichannel potential fields.
% Electroencephalography and Clinical Neurophysiology, 1980, 48, 609-621.
% Lehmann, D. & Skrandies, W. Spatial analysis of evoked potentials in
% man--a review. Progress in Neurobiology, 1984, 23, 227-250.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

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

ver = '$Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $';
fprintf('EEG_GFP [v %s]\n',ver(11:15));

if nargin < 1,
    help eeg_gfp;
    return
end

if ~exist('eeg','var'),
    error('no eeg input');
end
if isempty(eeg),
    error('no eeg input');
end

[Nt,Ne] = size(eeg);

fprintf('...%d time points, %d sensors\n',Nt,Ne);

if ~exist('avg','var'), avg = 1; end
if isempty(avg), avg = 1; end

% calculate average reference data ('u' in Lehmann & Skrandies,1984,p.235)
if avg,
    fprintf('...calculating average reference correction\n');
    avgref = mean(eeg,2);
    avgref = repmat(avgref,1,Ne);
    u = eeg - avgref;
else
    fprintf('...assuming eeg is average reference corrected\n');
    u = eeg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global field power
fprintf('...computing global field power\n')


% Lehmann & Skrandies (1984, p. 235):
% We use numerical procedures which assess the degree of 'hilliness' (or
% 'relief', or electrical strength) of the fields (Lehmann, 1971, 1972;
% Lehmann and Skrandies, 1980). The principal approach is to consider all
% possible potential differences in the field (for n electrodes, n*(n - 1))
% with equal weight, and thus compute the reference-free, mean potential
% difference (global field power) at each moment in time using formula (1)
% GFP (reference-free) = sqrt( (1/2n)* sum_i(sum_j(ui - uj)^2) ) (1)
% where n is the number of electrodes which measure the potentials 
% ei and ej; i, j = 1 . . . n; the observed voltages are 
% ui = ei - common reference. We note that meaningful computation of
% electric field power is based on data from electrodes spaced about
% equidistantly over a reasonably large scalp recording area.


%         % As defined in Lehmann & Skrandies (1980,1984)
%         % This takes forever to run, so the code below is adapted from
%         % this, but the code here is a useful reference.
%         gfp = zeros(Nt,1);
%         for t = 1:Nt,
%             % u = potential - average reference
%             u = e(t,:) - mean(e,2);
%             sumsqdif = 0;
%             for i=1:Ne,
%                 for j=1:Ne,
%                     sumsqdif = sum( [sumsqdif, (u(i) - u(j))^2] );
%                 end
%             end
%             gfp(t) = sqrt( sumsqdif / (2*Ne) );
%         end



fprintf('...calculating electrode potential differences\n');
% allocate sum variable
sumsqdif = zeros(Nt,1);
for i=1:Ne,
    
    progress = sprintf('...electrode %6d of %6d\n',i,Ne);
    if i>1, progress = [repmat('\b',1,length(progress)),progress]; end
    fprintf(progress);

    % each iteration, extract an electrode waveform and then remove
    % it from further consideration, so
    % a) extract electrode i
    ui = u(:,i);
    ui = repmat(ui,1,length(i:Ne));    
    % b) only consider this and the remaining electrodes
    uj = u(:,i:Ne);
    % differences of this waveform with all others (the
    % difference with itself is zero, it will not add to the sum)
    sqdif = (ui - uj).^2;
    % take the sum over electrodes
    sqdif = sum( sqdif, 2);
    % cumulative pairwise differences
    sumsqdif = sum([sumsqdif, sqdif], 2);
end
% Global Field Power is the root mean squared differences
% between all possible source locations; note here that we only
% divide by Ne because of the above progressive difference
% calculations, wheras Lehmann & Skrandies use 2*Ne because they
% calculate all pairwise differences.
gfp = sqrt(sumsqdif / Ne);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global dissimilarity index
% This requires gfp for scaling topographic maps

% Lehmann & Skrandies (1984, p. 240):
% For the comparison of the general configuration of two maps by a single
% number indicator we have proposed the measure of 'shape dissimilarity'
% (Lehmann and Skrandies, 1980a; Skrandies and Lehmann, 1982b). The two
% maps to be compared are first scaled for equal voltage range or global
% field power to omit differences in overall amplitude, since only shape of
% the distribution should be considered. The index of dissimilarity is
% computed as the mean standard deviation (MSD) over all recording points
% of the average voltages between the two maps at each recording point,
% referred to the average reference of the given map, using formula (3):
%
% dissimilarity index = MSD = (1/2n)* sum( sqrt(u(m1(i)) - u(m2(i)))^2) (3)
%

fprintf('...computing global dissimilarity index\n')

% first normalize the eeg by the gfp
map = eeg ./ repmat(gfp,1,Ne);

% differences of time (t+1 - t).
dif = diff(map,1,1);
% replicate first row, as if calculated backwards (omit the negative sign
% because this is all recitified anyway).
dif = [dif(1,:); dif];
rds = sqrt(dif.^2);
% average across electrodes (spatial dimension)
gd = mean(rds,2); 

% Lehmann & Skrandies (1984) divide by 2*Ne, but the routine above (and
% below) calculates the time differences forward only (not forward and
% back), so it is only necessary to divide by Ne (ie, take the mean over
% electrodes).


% % As defined in Lehmann & Skrandies (1984)
% % This takes a long time, but the code here is a useful reference.
% gd = zeros(Nt,1);
% for t = 2:Nt,
%     % maps are normalized by the Global Field Power
%     map1 = e(t-1,:) / gfp(t-1);
%     map2 = e(t-0,:) / gfp(t-0);
%     sumsqdif = 0;
%     for i=1:Ne,
%         sumsqdif = sum([sumsqdif, sqrt( (map1(i) - map2(i))^2 )] ));
%     end
%     gd(t) = sumsqdif / (2*Ne); % should this be Ne not 2*Ne?
% end
% gd = [gd(1); gd);

return
