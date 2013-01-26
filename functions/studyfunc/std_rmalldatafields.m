% std_rmalldatafields - remove all data fields from STUDY (before saving
%                       it for instance.
%            
% Usage:
%   STUDY = std_rmalldatafields(STUDY, type);
%
% Input:
%   STUDY   - EEGLAB study structure
%
% Optional input:
%   type    - ['chan'|'clust'|'both'] remove from changrp channel location
%             structure, cluster structure or both. Default is 'both'.
%
% Ouput:
%   STUDY   - updated EEGLAB study structure
%
% Author: Arnaud Delorme, CERCO/CNRS, UCSD, 2010-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function STUDY = std_rmalldatafields(STUDY, chanorcomp);

if nargin < 1
    help std_rmalldatafields;
    return;
end;
if nargin < 2
    chanorcomp = 'both';
end;

fields = { 'erpdata'  'erptimes'  'erpdatatrials'  'erptrialinfo' ...
               'specdata' 'specfreqs' 'specdatatrials' 'spectrialinfo' ...
               'erspdata' 'erspbase'  'erspfreqs'      'ersptimes' 'erspdatatrials' 'ersptrialinfo'  ...
               'topo' 'topox' 'topoy' 'topoall' 'topopol' ...
               'itcdata'  'itcfreqs'  'itctimes'       'itcdatatrials' 'itctrialinfo'  ...
               'erpimdata' 'erpimtrials' 'erpimevents' 'erpimtimes' 'erspsubjinds' 'itcsubjinds' 'dipoles' ...
               'data' 'datatimes' 'datasortvals' 'datacontinds' 'centroid' };
for ff = 1:length(fields)
    if strcmpi(chanorcomp, 'data') || strcmpi(chanorcomp, 'both')
        if isfield(STUDY.changrp, fields{ff})
            STUDY.changrp = rmfield(STUDY.changrp, fields{ff} );
        end;
    end;
    if strcmpi(chanorcomp, 'clust') || strcmpi(chanorcomp, 'both')
        if isfield(STUDY.cluster, fields{ff})
            STUDY.cluster = rmfield(STUDY.cluster, fields{ff} );
        end;
    end;
end;
