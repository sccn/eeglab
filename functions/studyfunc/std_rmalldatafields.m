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

function STUDY = std_rmalldatafields(STUDY, chanorcomp);

if nargin < 1
    help std_rmalldatafields;
    return;
end
if nargin < 2
    chanorcomp = 'both';
end

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
        end
    end
    if strcmpi(chanorcomp, 'clust') || strcmpi(chanorcomp, 'both')
        if isfield(STUDY.cluster, fields{ff})
            STUDY.cluster = rmfield(STUDY.cluster, fields{ff} );
        end
    end
end
