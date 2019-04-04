% std_reset() - Remove all preloaded measures from STUDY
%
% Usage:
%   >> STUDY = std_reset(STUDY);
%
% Inputs:
%   STUDY - EEGLAB STUDY structure
%
% Outputs:
%   STUDY - EEGLAB STUDY structure
%
% Author: Arnaud Delorme, CERCO/CNRS & SCCN, INC, UCSD, 2009-

% Copyright (C) Arnaud Delorme, arno@sccn.ucsd.edu
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

function STUDY = std_reset(STUDY)

if nargin < 1
    help std_reset;
    return;
end

fields = { 'erpdata' 'erptimes' 'specdata' 'specfreqs' 'erspdata' ...
           'ersptimes' 'erspfreqs' 'itcdata' 'itctimes' 'itcfreqs' ...
           'topo' 'topox' 'topoy' 'topoall' 'topopol' 'dipole' };
for ind = 1:length(fields)
    if isfield(STUDY.cluster, fields{ind})
        STUDY.cluster = rmfield(STUDY.cluster, fields{ind});
    end
    if isfield(STUDY, 'changrp')
        if isfield(STUDY.changrp, fields{ind})
            STUDY.changrp = rmfield(STUDY.changrp, fields{ind});
        end
    end
end
