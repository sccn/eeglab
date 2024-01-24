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
%  'chanlabels] - [cell fo string] labels for each of the channels.
%  'cellinfofields' - [cell of string] Only use some of the fields of the
%                 cellinfo input.
%  'csvfile'    - [string] CSV file name to save the table. By default, the
%                 table is not saved.
%
% Output:
%  tab          - MATLAB table
%
% Example:
%  % assumming a STUDY with at least one design and that channel Fz exists
%  [~,spec,xvals,~,~,~,~,specinfo] = std_readdata(STUDY, ALLEEG, 'design', ...
%                       1, 'channels', { 'E1' }, 'datatype', 'spec');
%  tab = std_cell2table(xvals, [], spec, specinfo);
%
% Author: Arnaud Delorme, CERCO, 2024-

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
    'design'     'struct'    {}       {};
    'csvfile'    'string'    {}       '';
    'cellinfofields' 'cell'      {}       {}}, 'std_statcond2table');
if ischar(g), error(g); end

if isempty(cellvals)
    error('Statistical array cannot be empty')
end

if size(cellvals,1) > 1 && size(cellvals,2) > 1 
    if length(STUDY.design(g.design).variable) ~= 2
        error('The STUDY design must contains the same number of variables as the number of dim in the statistical cell array')
    else
        if size(cellvals,1) ~= length(STUDY.design(g.design).variable(1).value)
            error('The STUDY design first variable must contain the same number of values as the number of items in the first dim in the statistical cell array')
        end
        if size(cellvals,2) ~= length(STUDY.design(g.design).variable(2).value)
            error('The STUDY design second variable must contain the same number of values as the number of items in the second dim in the statistical cell array')
        end
    end

    %create table
    tab = table([]);
    count = 1;
    for iCell1 = 1:size(cellvals,1)
        for iCell1 = 1:size(cellvals,2)
            if ~isempty(g.infoField)
                datInds = g.setinds{iCell1, iCell2};
                fieldData = STUDY.datasetinfo(datInds);
            else
                fieldData = [];
            end
                
            [tab, count] = adddata(tab, count, [iCell1 iCell2], fieldData, g.infoField, cellvals{iCell1,iCell2});
        end
    end

    % update column names
    curcols = tab.Properties.VariableNames;
    newcols = { STUDY.design(g.deisgn).variable.label };
    countCols = length(newcols);
    if ~isempty(g.infoField)
        for iCol = 1:length(infoField)
            countCols = countCols + 1;
            newcols{countCols} = infoField{iCol};
        end
    end
    if ~isempty(g.dimension)
        for iCol = 1:length(g.dimension)
            countCols = countCols + 1;
            newcols{countCols} = g.dimension{iCol};
        end
    end
    for iCol = countCols+1:length(curcols)
        newcols{iCol} = [ 'Var' int2str(iCol) ];
    end
    tab.Properties.VariableNames = newcols;

else
    if length(STUDY.design(g.design).variable) ~= 1)
        error('The STUDY design must contains the same number of variables as the number of dim in the statistical cell array')
    end
    if length(cellvals) ~= length(STUDY.design(g.design).variable(1).value)
        error('The STUDY design first variable must contain the same number of values as the number of items in the first dim in the statistical cell array')
    end
end

% add data to the table
% ---------------------
function [tab, count] = adddata(tab, count, infoData, infoFields, firstCols, dataArray)

if ndims(dataArray) > 3
    error('Cannot process 4D array')
end
countCols = 1;

for iData1 = 1:size(dataArray,1)
    for iData2 = 1:size(dataArray,2)
        for iData3 = 1:size(dataArray,3)

            % write cell index information
            for iCol = 1:length(firstCols)
                tab(count, countCols) = firstCols(iCol);
                countCols = countCols+1;
            end

            % write field information
            if size(dataArray,1) > 1, iLastDim = iData1; end
            if size(dataArray,2) > 1, iLastDim = iData2; end
            if size(dataArray,3) > 1, iLastDim = iData3; end
            for iCol = 1:length(infoFields)
                tab(count, countCols) = infoData(iLastDim).(infoFields{iCol});
                countCols = countCols+1;
            end

            % write data
            if size(dataArray,1) > 1, tab(count, countCols) = iData1; countCols = countCols+1; end
            if size(dataArray,2) > 1, tab(count, countCols) = iData2; countCols = countCols+1; end
            if size(dataArray,3) > 1, tab(count, countCols) = iData3; countCols = countCols+1; end
            tab(count, countCols) = dataArray(iData1, iData2, iData3 );
        end
    end
end


