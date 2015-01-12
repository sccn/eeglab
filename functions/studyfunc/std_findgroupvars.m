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
        end;
    end;
    
    for iSubj = 1:length(subjVals)
        if ~isempty(subjVals{iSubj})
            if isstr(subjVals{iSubj}{1}) subjVals{iSubj} = unique(subjVals{iSubj});
            else                         subjVals{iSubj} = mattocell(unique([subjVals{iSubj}{:}]));
            end;
        end;
    end;
    if all( cellfun(@length, subjVals) == 1 )
        groupvar{end+1} = fields{iField};
    end;
end;