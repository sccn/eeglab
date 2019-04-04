% newtimefitc() - Function to compute inter-trial coherence (phase locking
%                 factor is another name for it) from single trial spectral
%                 estimates.
%
% Usage:
%   >>  itc = newtimefitc(tfvals, itctype);
%
% Inputs:
%   tfvals - [3-D or 4-D array] single-trial spectral estimates
%           [freqs x times x trials] or [channels x freqs x times x trials)
%   itctype  - ['coher'|'phasecoher'|'phasecoher2'] Compute either linear
%              coherence ('coher') or phase coherence ('phasecoher').
%              Originall called 'phase-locking factor' {default: 'phasecoher'}
%
% Outputs:
%   itc - (nfreqs,timesout) matrix of complex inter-trial coherencies.
%         tc is complex -- ITC magnitude is abs(itc); ITC phase in radians
%         is angle(itc), or in deg phase(itc)*180/pi.
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, Sep 2016

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2016, arno@sccn.ucsd.edu
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

function [itcvals] = newtimefitc(tfdecomp, itctype);

nd = max(3,ndims(tfdecomp));
switch itctype
    case 'coher',
        try,
            itcvals = sum(tfdecomp,nd) ./ sqrt(sum(tfdecomp .* conj(tfdecomp),nd) * size(tfdecomp,nd));
        catch, % scan rows if out of memory
            for index =1:size(tfdecomp,1)
                itcvals(index,:,:) = sum(tfdecomp(index,:,:,:),nd) ./ sqrt(sum(tfdecomp(index,:,:,:) .* conj(tfdecomp(index,:,:,:)),nd) * size(tfdecomp,nd));
            end
        end
    case 'phasecoher2',
        try,
            itcvals = sum(tfdecomp,nd) ./ sum(sqrt(tfdecomp .* conj(tfdecomp)),nd);
        catch, % scan rows if out of memory
            for index =1:size(tfdecomp,1)
                itcvals(index,:,:) = sum(tfdecomp(index,:,:,:),nd) ./ sum(sqrt(tfdecomp(index,:,:,:) .* conj(tfdecomp(index,:,:,:))),nd);
            end
        end
    case 'phasecoher',
        try,
            itcvals = sum(tfdecomp ./ sqrt(tfdecomp .* conj(tfdecomp)) ,nd) / size(tfdecomp,nd);
        catch, % scan rows if out of memory
            for index =1:size(tfdecomp,1)
                itcvals(index,:,:) = sum(tfdecomp(index,:,:,:) ./ sqrt(tfdecomp(index,:,:,:) .* conj(tfdecomp(index,:,:,:))) ,nd) / size(tfdecomp,nd);
            end
        end
end
