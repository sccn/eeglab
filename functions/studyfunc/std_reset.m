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
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function STUDY = std_reset(STUDY)

if nargin < 1
    help std_reset;
    return;
end;

fields = { 'erpdata' 'erptimes' 'specdata' 'specfreqs' 'erspdata' ...
           'ersptimes' 'erspfreqs' 'itcdata' 'itctimes' 'itcfreqs' ...
           'topo' 'topox' 'topoy' 'topoall' 'topopol' 'dipole' };
for ind = 1:length(fields)
    if isfield(STUDY.cluster, fields{ind})
        STUDY.cluster = rmfield(STUDY.cluster, fields{ind});
    end;
    if isfield(STUDY, 'changrp')
        if isfield(STUDY.changrp, fields{ind})
            STUDY.changrp = rmfield(STUDY.changrp, fields{ind});
        end;
    end;
end;
