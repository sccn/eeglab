% std_addvarlevel() - add design variable level
%
% Usage:
%   >> [STUDY] = std_makedesign(STUDY);   
%   >> [STUDY] = std_makedesign(STUDY, designind);   
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   designind   - [integer > 0] index (number) of the new STUDY design {default: all}
%
% Author: Arnaud Delorme, Institute for Neural Computation UCSD, 2020-

% Copyright (C) Arnaud Delorme, arno@ucsd.edu
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

function STUDY = std_addvarlevel(STUDY, designind)

if nargin < 1
    help std_addvarlevel;
    return;
end
if nargin < 2
    designind = [1:length(STUDY.design)];
end

%     if ismember(ff, 'condition')
%         disp('Warning: event field named "condition" is ignored because it is used 

varSubject = { STUDY.datasetinfo.subject };
subjects   = unique(varSubject);

for iDes = 1:length(STUDY.design)
    for iVar = 1:length(STUDY.design(iDes).variable)
        if ~isfield(STUDY.datasetinfo(1), STUDY.design(iDes).variable(iVar).label)
            % level is one if not a field of STUDY.datasetinfo
            STUDY.design(iDes).variable(iVar).level = 'one';
        else
            % otherwise level is two
            STUDY.design(iDes).variable(iVar).level = 'two';
            
            % unless the same subject has different values for that variable
            varValue = { STUDY.datasetinfo.(STUDY.design(iDes).variable(iVar).label) };
            for iSubj = 1:length(subjects)
                inds = strmatch(subjects{iSubj}, varSubject, 'exact');
                tmpVals = varValue(inds);
                if ~ischar(tmpVals{1}), tmpVals = [ tmpVals{:} ]; end
                if   length(unique(tmpVals)) > 1
                    STUDY.design(iDes).variable(iVar).level = 'one';    
                end
            end
        end
    end
end

