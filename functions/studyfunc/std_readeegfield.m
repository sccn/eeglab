% STD_READEEGFIELD - read field from ALLEEG structure
% Usage:
%   >> [dataOut, ~, xvals] = std_readeegfield(datasetinfo, ALLEEG, ...
%                                        designvar, fieldName, 'key', val)
% Inputs:
%  datasetinfo - STUDY.datasetinfo structure array corresponding to the 
%                ALLEEG variable below
%  ALLEEG      - vector of EEG datasets
%  designvar   - [struct] design variable. For example STUDY.design(1)
%  fieldName   - [string or cell] fieldnames to return. For example 'pnts' will
%                return the number of points EEG.pnts. { 'roi' 'MIM' }
%                will return the MIM array.
%
% Optional inputs:
%  'xvalfield'   - [string or cell] fieldnames to return.
%  'xvallimit'   - [min max] select data within the specific range for the 
%                  first dimension of the output
%  'xvalaverage' - ['on'|'off'] average the first dimension of the output 
%                  and squeeze the output. Default is 'off'.
%
% Output:
%  dataOut   - cell array of value
%  xvals     - x values
%  note: the output matches the output of std_readfile because this
%  function can be calle by STD_READDATA
%
% Example:
%  std_precomp(STUDY, ALLEEG, { ALLEEG(1).chanlocs.labels }, 'erp', 'on');
%  [erp times] = std_readdata(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(1).labels });
%
% Author: Arnaud Delorme, CERCO, 2006-

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

function [dataOut, params, xvals, yvals, events] = std_readeegfield(datasetinfo, ALLEEG, designvar, fieldName, varargin)

dataOut = {};
params  = [];
xvals   = [];
yvals   = [];
events  = [];
if nargin < 2
    help std_readroifield;
    return
end

% std_readfile parameters could be useful but most are ignored
opt = finputcheck(varargin, { ...
    'xvalfield'        ''         []    {}; 
    'xvalaverage'      'string'   { 'on','off' }  'off';
    'xvallimits'       'real'     []    [] }, 'std_readeegfield', 'ignore');

if ~iscell(fieldName)
    fieldName = { fieldName };
end
if ~iscell(opt.xvalfield)
    opt.xvalfield = { opt.xvalfield };
end

% check data
for iDes = 1:length(designvar)
    if ~isfield(datasetinfo, designvar(iDes).label)
        error('This function read fields of the EEG structure so it can only handle STUDY designs comparing datasets (not trials within datasets)')
    end
end

if length(designvar) == 1
    for iVar1 = 1:length(designvar(1).value)

        dataOut{iVar1} = getDataCell( designvar(1).value{iVar1},  { datasetinfo.(designvar(1).label) }, ALLEEG, fieldName);

    end
else
    for iVar1 = 1:length(designvar(1).value)
        for iVar2 = 1:length(designvar(2).value)

            dataOut{iVar1, iVar2} = getDataCell2( designvar(1).value{iVar1}, { datasetinfo.(designvar(1).label) }, ...
                                                  designvar(2).value{iVar2}, { datasetinfo.(designvar(2).label) }, ALLEEG, fieldName);

        end
    end
end

% return xvals and yvals
ind = find(~cellfun(@isempty, dataOut(:)));
if ~isempty(ind)
    ind = ind(1);
    if size(dataOut{ind},1) > 1
        xvals = 1:size(dataOut{ind},1);
        if size(dataOut{ind},2) > 2
            yvals = 1:size(dataOut{ind},2);
        end
    elseif size(dataOut{ind},2) > 2
        xvals = 1:size(dataOut{ind},2);
    end
    if ~isempty(opt.xvalfield)
        xvals = getfield( ALLEEG(1), opt.xvalfield{:} );
    end
end

% select xvals
if ~isempty(opt.xvallimits)
    [xvals, indBegin, indEnd] = indicesselect(xvals, opt.xvallimits);
    for iDat = 1:length(dataOut(:))
        if ~isempty(dataOut{iDat})
            dataOut{iDat} = dataOut{iDat}(indBegin:indEnd,:,:,:);
        end
    end
end
if strcmpi(opt.xvalaverage, 'on')
    for iDat = 1:length(dataOut(:))
        if ~isempty(dataOut{iDat})
            dataOut{iDat} = squeeze(mean(dataOut{iDat}));
            if size(dataOut{iDat},1) == 1 
                dataOut{iDat} = transpose(dataOut{iDat});
            end
        end
    end
end

% get data from ALLEEG 1
% ----------------------
function dataOut = getData1( varValue, varList, ALLEEG, fieldName)

if ischar(varList{1}) && ~ischar(varValue), varValue = num2str(varValue); end
datInd = find(cellfun(@(x)isequal(x, varValue), varList));
if isempty(datInd)
    dataOut = {[]};
else
    for iDat = 1:length(datInd)
        dataOut{iDat} = getfield(ALLEEG(datInd(iDat)), fieldName{:});
    end
end

% get data from ALLEEG 2
% ----------------------
function dataOut = getData2( varValue1, varList1, varValue2, varList2, ALLEEG, fieldName)

if ischar(varList1{1}) && ~ischar(varValue1), varValue1 = num2str(varValue1); end
if ischar(varList2{1}) && ~ischar(varValue2), varValue2 = num2str(varValue2); end
datInd1 = cellfun(@(x)isequal(x, varValue1), varList1);
datInd2 = cellfun(@(x)isequal(x, varValue2), varList2);
datInd = find( datInd1 & datInd2 );

if isempty(datInd)
    dataOut = {[]};
else
    for iDat = 1:length(datInd)
        dataOut{iDat} = getfield(ALLEEG(datInd(iDat)), fieldName{:});
    end
end

% combine cell arrays
% -------------------
function dataOut = combinecell( dataTmp )
persistent warnflag

if isempty(dataTmp)
    dataOut = {};
elseif all(cellfun(@isempty, dataTmp))
    dataOut = [];
else
    if any(~cellfun(@isempty, dataTmp)) && isempty(warnflag)
        warning('Some condition have no data for at least one subject; check design')
        warnflag = true;
    end

    % average values
    count = 0;
    dataOut = {};
    for iVar = 1:length(dataTmp)
        if ~isempty(dataTmp{iVar})
            count = count + 1;
            if isempty(dataOut)
                dataOut = dataTmp{iVar};
            else
                dataOut =  + dataTmp{iVar};
            end
        end
    end
    dataOut = dataOut/count;
end

% get data from cells 1
% ---------------------
function dataOut = getDataCell( varValue, varList, ALLEEG, fieldName)

if ~iscell(varValue)
    dataTmp = getData1( varValue, varList, ALLEEG, fieldName);
else
    count = 1;
    for iVar = length(varValue):-1:1
        dataTmp2 = getData1( varValue{iVar}, varList, ALLEEG, fieldName);
        dataTmp(count:count+length(dataTmp2)-1) = dataTmp2;
        count = count+1;
    end
end
dataOut = combinecell(dataTmp);

% get data from cells 2
% ---------------------
function dataOut = getDataCell2( varValue1, varList1, varValue2, varList2, ALLEEG, fieldName)

if ~iscell(varValue1) && ~iscell(varValue2)
    dataTmp = getData2( varValue1, varList1, varValue2, varList2, ALLEEG, fieldName);
else
    count = 1;
    if iscell(varValue1) && ~iscell(varValue2)
        for iVar = length(varValue1):-1:1
            dataTmp2 = getData2( varValue1{iVar}, varList1, varValue2, varList2, ALLEEG, fieldName);
            dataTmp(count:count+length(dataTmp2)-1) = dataTmp2;
            count = count+1;
        end
    elseif ~iscell(varValue1) && iscell(varValue2)
        for iVar = length(varValue2):-1:1
            dataTmp2 = getData2( varValue1, varList1, varValue2{iVar}, varList2, ALLEEG, fieldName);
            dataTmp(count:count+length(dataTmp2)-1) = dataTmp2;
            count = count+1;
        end
    elseif iscell(varValue1) && iscell(varValue2)
        for iVar1 = length(varValue1):-1:1
            for iVar2 = length(varValue2):-1:1
                dataTmp2 = getData2( varValue1{iVar1}, varList1, varValue2{iVar2}, varList2, ALLEEG, fieldName);
                dataTmp(count:count+length(dataTmp2)-1) = dataTmp2;
                count = count+1;
            end
        end
    else
        error('Unknown configuration')
    end

end
dataOut = combinecell(dataTmp);

% select range for xvals
function [measureRange, indBegin, indEnd] = indicesselect(measureRange, measureLimits)
indBegin = 1;
indEnd   = length(measureRange);
if ~isempty(measureRange) && ~isempty(measureLimits) && (measureLimits(1) > measureRange(1) || measureLimits(end) < measureRange(end))
    indBegin   = min(find(measureRange >= measureLimits(1)));
    indEnd     = max(find(measureRange <= measureLimits(end)));
    measureRange = measureRange(indBegin:indEnd);
end
