% STD_CELL2TABLE - convert cell output to table
%
% Usage:
%   >> tab = std_cell2table(xvals, yvals, cellvals, cellinfo, 'key', val)
%
% Inputs:
%  xvals    - [vector] values along the first axis
%  yvals    - [vector] values along the second axis (if any)
%  cellvals - [cell array] second output of STD_READATA
%  cellinfo - [cell array] last output of STD_READATA
%
% Optional inputs:
%  'dimensions' - [cell of string] dimensions of data in cellvals
%                 variable. For example { 'time' 'channel' 'participants' }
%                 These are used as column names in the output table. This
%                 function will infer the last dim to be participants or
%                 trials. By default, the columns will be named 
%                 { 'var1' 'var2' ... }
%  'design'     - [struct] STUDY design corresponding to the data
%  'chanlabels  - [cell fo string] labels for each of the channels.
%  'cellinfofields' - [cell of string] Only use some of the fields of the
%                 cellinfo input. Use { '' } to not include any field.
%                 Default, and {} to include all fields.
%  'csvfile'    - [string] CSV file name to save the table. By default, the
%                 table is not saved.
%
% Output:
%  tab          - MATLAB table
%
% Example:
%  % assumming a STUDY with at least one design and that channel Fz exists
%  [~,spec,xvals,~,~,~,~,specinfo] = std_readdata(STUDY, ALLEEG, 'design', ...
%                       1, 'channels', { 'E1' }, 'datatype', 'spec', 'freqrange', [1 20]);
%  tab = std_cell2table(xvals, [], spec, specinfo, 'design', STUDY.design(1), ...
%         'cellinfofields', { 'run' }, 'dimensions', {'frequency', 'power'});
%
% % Reading custom field and creating a table out of it
% % the code below assumes that alpha asymetry has been computed on all the datasets
% [STUDY,aa,~,~,~,~,~,info] = std_readdata(STUDY, ALLEEG, 'design', 1, 'customread', ...
%  'std_readeegfield', 'customparams', {{ 'etc', 'eegstats', 'alpha_asymmetry' }}, 'ndim', 1);
%  tab = std_cell2table([], [], aa, info, 'design', STUDY.design(1), 'dimensions', {'alpha_asymetry'});
%
% Author: Arnaud Delorme, CERCO, 2024-

%
%  [~,spec,xvals,~,~,~,~,specinfo] = std_readdata(STUDY, ALLEEG, 'design', ...
%                       1, 'channels', { 'E1' 'E2' }, 'datatype', 'spec', 'freqrange', [1 20]);
%  tab = std_cell2table(xvals, [], spec, specinfo, 'design', STUDY.design(1), ...
%         'cellinfofields', { 'run' }, 'dimensions', {'frequency', 'chan', 'power'}, 'chanlabels', { 'E1' 'E2' });
%
%  [~,spec,xvals,~,~,~,~,specinfo] = std_readdata(STUDY, ALLEEG, 'design', ...
%                       1, 'channels', { 'E1' 'E2' 'E3' 'E4' 'E5' }, 'datatype', 'spec', 'freqrange', [10]);
%  tab = std_cell2table(xvals, [], spec, specinfo, 'design', STUDY.design(1), ...
%         'cellinfofields', { 'run' }, 'dimensions', {'frequency', 'chan', 'power'}, 'chanlabels', { 'E1' 'E2' });


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


function tab = std_cell2table(xvals, yvals, cellvals, cellinfo, varargin)

if nargin < 2
    help std_statcond2table;
    return
end

g = finputcheck(varargin, { ...
    'dimensions' 'cell'      {}       {};
    'chanlabels' 'cell'      {}       {};
    'design'     'struct'    {}       struct([]);
    'csvfile'    'string'    {}       '';
    'cellinfofields' 'cell'      {}       {''}}, 'std_statcond2table');
if ischar(g), error(g); end

if isempty(cellvals)
    error('Statistical array cannot be empty')
end

% only keep relevant fields
ind = find(~cellfun(@isempty, cellinfo(:))); % non-empty info
allFields = fieldnames(cellinfo{ind(1)});
if ~isempty(g.cellinfofields)
    if isempty(g.cellinfofields{1})
        g.cellinfofields = {};
    end
    if ~contains(g.cellinfofields, 'subject')
        g.cellinfofields{end+1} = 'subject';
    end
    rmFields = setdiff(allFields, g.cellinfofields);
    for iCell1 = 1:size(cellinfo,1)
        for iCell2 = 1:size(cellinfo,2)
            tmpinfo = cellinfo{iCell1,iCell2};
            if ~isempty(tmpinfo)
                tmpinfo = rmfield(tmpinfo, rmFields);
                cellinfo{iCell1,iCell2} = tmpinfo;
                g.cellinfofields = fieldnames(tmpinfo);
            end
            
        end
    end
else
    g.cellinfofields = allFields;
end

% check input size
% ----------------
warning('off', 'MATLAB:table:RowsAddedExistingVars')
if ~isempty(g.design)
    if size(cellvals,1) > 1 && size(cellvals,2) > 1
        if length(g.design.variable) ~= 2
            error('The STUDY design must contains the same number of variables as the number of dim in the statistical cell array')
        else
            if size(cellvals,1) ~= length(g.design.variable(1).value)
                error('The STUDY design first variable must contain the same number of values as the number of items in the first dim in the statistical cell array')
            end
            if size(cellvals,2) ~= length(g.design.variable(2).value)
                error('The STUDY design second variable must contain the same number of values as the number of items in the second dim in the statistical cell array')
            end
        end
    end
else
    if length(g.design.variable) ~= 1
        error('The STUDY design must contains the same number of variables as the number of dim in the statistical cell array')
    end
    if length(cellvals) ~= length(g.design.variable(1).value)
        error('The STUDY design first variable must contain the same number of values as the number of items in the first dim in the statistical cell array')
    end
end

% create table
% ------------
if size(cellvals,1) > 1 && size(cellvals,2) > 1
    tab = table([]);
    count = 1;
    for iCell1 = 1:size(cellvals,1)
        for iCell2 = 1:size(cellvals,2)
            if ~isempty(g.design)
                designVals = { g.design.variable(1).value{iCell1}, g.design.variable(2).value{iCell2} };
            else
                designVals = { iCell1 iCell2 };
            end
            [tab, count] = adddata(tab, count, designVals, xvals, yvals, g.chanlabels, cellvals{iCell1,iCell2}, cellinfo{iCell1,iCell2});
        end
    end

    % update column names
    curcols = tab.Properties.VariableNames;
    newcols = { g.design.variable.label };
    countCols = length(newcols);
    if ~isempty(g.cellinfofields)
        for iCol = 1:length(g.cellinfofields)
            countCols = countCols + 1;
            if ~contains(newcols, g.cellinfofields{iCol})
                newcols{countCols} = g.cellinfofields{iCol};
            else
                newcols{countCols} = [g.cellinfofields{iCol} 'x'];
            end
        end
    end
    if ~isempty(g.dimensions)
        for iCol = 1:length(g.dimensions)
            countCols = countCols + 1;
            newcols{countCols} = g.dimensions{iCol};
        end
    end
    for iCol = countCols+1:length(curcols)
        newcols{iCol} = [ 'Var' int2str(iCol) ];
    end
    tab.Properties.VariableNames = newcols;

end

% add data to the table
% ---------------------
function [tab, count] = adddata(tab, count, designVals, xvals, yvals, chanlabels, dataArray, infoData)

if ndims(dataArray) > 4
    error('Cannot process 4D array')
end
fields = fieldnames(infoData);
if size(dataArray,1) == 1 && size(dataArray,3) == 1, dataArray = dataArray'; end
if size(dataArray,1) > 1, totDim = 1; end
if size(dataArray,2) > 1, totDim = 2; end
if size(dataArray,3) > 1, totDim = 3; end
if size(dataArray,4) > 1, totDim = 4; end

for iData1 = 1:size(dataArray,1)
    for iData2 = 1:size(dataArray,2)
        for iData3 = 1:size(dataArray,3)
            for iData4 = 1:size(dataArray,4)
                countCols = 1;

                % write cell index information
                for iCol = 1:length(designVals)
                    if isempty(tab) && ischar(designVals{iCol}), tab = table({}); end
                    tab(count, countCols) = designVals(iCol);
                    countCols = countCols+1;
                end

                % write field information
                iLastDim = [iData1 iData2 iData3 iData4];
                iLastDim = iLastDim(totDim);
                for iField = 1:length(fields)
                    tab(count, countCols) = { infoData(iLastDim).(fields{iField}) };
                    countCols = countCols+1;
                end

                % write data
                if totDim >= 2, tab(count, countCols) = { xvals(iData1) }; countCols = countCols+1; end % subject is already added above
                if totDim == 3
                    if ~isempty(chanlabels), tab(count, countCols) = chanlabels(iData2); 
                    elseif ~isempty(yvals) , tab(count, countCols) = { yvals(iData2) }; 
                    else                     tab(count, countCols) = { iData2 }; 
                    end
                    countCols = countCols+1; 
                end
                if totDim == 4
                    if ~isempty(yvals),     tab(count, countCols) = { yvals(iData2) };  else tab(count, countCols) = { iData2 }; end
                    if isempty(chanlabels), tab(count, countCols) = chanlabels(iData3); else tab(count, countCols) = { iData3 }; end
                    countCols = countCols+2; 
                end
                tab(count, countCols) = { dataArray(iData1, iData2, iData3 ) };
                count = count+1;
            end
        end
    end
end


