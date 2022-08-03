% std_rmdat() - remvoe datasets from STUDY
%
% Usage: 
%   >> STUDY = std_rmdat(STUDY, ALLEEG, 'key', val);  
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
%   'pntsrange'   - [min max] minimum and maximum of samples. Default is
%                   [0 Inf] (no constraint)
%   'chanrange'   - [min max] minimum and maximum of channels. Default is
%                   [0 Inf] (no constraint)
%   'sraterange'  - [min max] minimum and maximum for sampling rate. Default is
%                   [0 Inf] (no constraint)
%   'trialrange'  - [min max] minimum and maximum of trials. Default is
%                   [1 Inf]
%   'rmvarvalues' - {'string' range} remove datasets having variable value
%                   in the selected range. May also be {'string' 'value'}
%                   for non-numerical variables.
%   'checkeventtype' - [cell|array|string] check event type are present.
%   'numeventrange' - [min max] range for number of event of type above.
%                    Default is [1 Inf].
%
% Inputs:
%   STUDY      - EEGLAB STUDY set updated. The fields which is created or
%                updated is STUDY.datasetinfo.trialinfo
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, July 2022

% Copyright (C) Arnaud Delorme arno@ucsd.edu
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

function [STUDY, ALLEEG, rmDats] = std_rmdat(STUDY, ALLEEG, varargin)

if nargin < 3
    help std_rmdat;
    return;
end

g = finputcheck( varargin, { ...
    'chanrange'       'integer'  {} [0 Inf]; ...
    'pntsrange'       'integer'  {} [0 Inf]; ...
    'sraterange'      'float'    {} [0 Inf]; ...
    'trialrange'      'float'    {} [1 Inf]; ...
    'checkeventtype'  ''         {} []; ...
    'numeventrange'    'integer'  {} 1; ...
    'rmvarvalues'     'cell' {} {} }, 'std_rmdat');
if isstr(g)
    error(g);
end

allPnts   = [ALLEEG.pnts];
allSrate  = [ALLEEG.pnts];
allChans  = [ALLEEG.nbchan];
allTrials = [ALLEEG.trials];

% check pnts range
rmDats = g.pntsrange(1)  > allPnts   | allPnts   > g.pntsrange(2);
rmDats = g.sraterange(1) > allSrate  | allSrate  > g.sraterange(2) | rmDats;
rmDats = g.chanrange(1)  > allChans  | allChans  > g.chanrange(2)  | rmDats;
rmDats = g.trialrange(1) > allTrials | allTrials > g.trialrange(2) | rmDats;

% check variable name values
if ~isempty(g.rmvarvalues)
    varName   = g.rmvarvalues{1};
    varValues = g.rmvarvalues{2};
    allValues = { STUDY.datasetinfo.(varName) };
    if ischar(varValues)
        allValues = cellfun(@num2str, allValues, 'uniformoutput', false);
        rmDats = rmDats | cellfun(@(x)isequal(x, varValues),  allValues);
    elseif length(varValues) ~= 2
        error('When providing nmumerical input for variable selection, there must be 2 values - min and max');
    else
        allValues = cellfun(@(x)fastif(ischar(x), str2double(x), x), allValues);
        keepVals = varValues(1) <= allValues & allValues <= varValues(2); % to deal with NaNs
        rmDats = rmDats | ~keepVals;
    end
end     

% check event type present
if ~isempty(g.checkeventtype)
    if ischar(g.checkeventtype) g.checkeventtype = { g.checkeventtype }; end
    rmDatEvents = zeros(1, length(ALLEEG));
    for iDat = 1:length(ALLEEG)
        if isfield(ALLEEG(iDat).event, 'type')
            curEventTypes = { ALLEEG(iDat).event.type };
            if ischar(curEventTypes{1})
                for iType = 1:length(g.checkeventtype)
                    inds = strmatch(g.checkeventtype{iType}, curEventTypes, 'exact');
                    if length(inds) < g.numeventrange(1) || length(inds) > g.numeventrange(2)
                        rmDatEvents(iDat) = true;
                    end
                end
            else
                % numerical or string event types
                curEventTypes = [ ALLEEG(iDat).event.type ];
                for iType = 1:length(g.checkeventtype)
                    inds = find(g.checkeventtype(iType) == curEventTypes);
                    if length(inds) < g.numeventrange(1) || length(inds) > g.numeventrange(2)
                        rmDatEvents(iDat) = true;
                    end
                end
            end
        end
    end
    rmDats = rmDats | rmDatEvents;
end

if sum(rmDats) > 0
    fprintf('%d dataset meet the criteria for removal and have been removed\n', sum(rmDats));
    ALLEEG(rmDats) = [];
    STUDY.datasetinfo(rmDats) = [];
    for iDat = 1:length(STUDY.datasetinfo)
        STUDY.datasetinfo(iDat).index = iDat;
    end
else
    fprintf('No dataset meet the criteria for removal and have been removed\n');
end
STUDY.subject = {};
STUDY = std_checkset(STUDY, ALLEEG);

for iDesign = 1:length(STUDY.design)
    STUDY.design(iDesign).cases.value = intersect(STUDY.design(iDesign).cases.value, STUDY.subject);
end
