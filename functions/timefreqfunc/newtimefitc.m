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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
