% std_findgroupvars() - find group vars in STUY structure
%
% Usage:
%   >>  
%
% Inputs:
%    
% Outputs:
%
% See also: 
%   pop_studydesign,
%
% Author: Arnaud Delorme 
%         Ramon Martinez-Cancino
%
% Copyright (C) 2014  Arnaud Delorme, 
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

function  groupvar = std_findgroupvars(STUDY)

fields = setdiff( fieldnames(STUDY.datasetinfo), { 'filepath' 'filename' 'subject' 'comps' 'trialinfo' 'index' });

subjects = unique( { STUDY.datasetinfo.subject } );
groupvar = {};
for iField = 1:length(fields)
    
    subjVals = cell(1, length(subjects));
    
    for iDat = 1:length(STUDY.datasetinfo)
        iSubj = strmatch( STUDY.datasetinfo(iDat).subject, subjects);
        if ~isempty(STUDY.datasetinfo(iDat).(fields{iField}))
            subjVals{iSubj}{end+1} = STUDY.datasetinfo(iDat).(fields{iField});
        end
    end
    
    for iSubj = 1:length(subjVals)
        if ~isempty(subjVals{iSubj})
            if ischar(subjVals{iSubj}{1}) subjVals{iSubj} = unique(subjVals{iSubj});
            else                         subjVals{iSubj} = mattocell(unique([subjVals{iSubj}{:}]));
            end
        end
    end
    if all( cellfun(@length, subjVals) == 1 )
        groupvar{end+1} = fields{iField};
    end
end
