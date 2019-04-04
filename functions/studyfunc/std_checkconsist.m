% std_checkconsist() - Create channel groups for plotting.
%
% Usage:    
%                >> boolval = std_checkconsist(STUDY, 'uniform', 'condition');   
%                >> boolval = std_checkconsist(STUDY, 'uniform', 'group');   
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   'uniform'  - ['condition'|'group'] check if there is one group
%                condition per subject 
% Outputs:
%   boolval    - [0|1] 1 if uniform
%
% Authors: Arnaud Delorme, CERCO, 2009

% Copyright (C) Arnaud Delorme, CERCO, arno@salk.edu
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

function [boolval npersubj] = std_checkconsist(STUDY, varargin);

if nargin < 3
    help std_checkconsist;
    return;
end

opt = struct(varargin{:});

if strcmpi(opt.uniform, 'condition')
    allvals = { STUDY.datasetinfo.condition };
    vallist = STUDY.condition;
elseif strcmpi(opt.uniform, 'group')
    allvals = { STUDY.datasetinfo.group };
    vallist = STUDY.group;
elseif strcmpi(opt.uniform, 'session')
    allvals = { STUDY.datasetinfo.session };
    allvals = cellfun(@num2str, allvals, 'uniformoutput', false);
    vallist = STUDY.session;
    if isempty(vallist), boolval = 1; return; end
    vallist = cellfun(@num2str, mattocell(vallist), 'uniformoutput', false);
else
    error('unknown option');
end

if isempty(vallist), boolval = 1; return; end

for index = 1:length(vallist)
    tmplist = strmatch( vallist{index}, allvals, 'exact');
    vallen(index) = length(unique( { STUDY.datasetinfo(tmplist).subject } ));
end
if length(unique(vallen)) == 1
     boolval = 1;
else boolval = 0;
end

