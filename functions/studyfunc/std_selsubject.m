% std_selsubject() - Helper function for std_erpplot(), std_specplot() 
%                    and std_erspplot() to select specific subject when
%                    plotting channel data.
% Usage:
%  >> data = std_selsubject( data, subject, setinds, allsubjects);
%
% Inputs:
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to compute mean ERPs statistics from a  
%           STUDY for epochs of 800 frames in two conditions from three  
%           groups of 12 subjects,
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
% subject - [string] subject name
% setinds - [cell array] set indices for each of the last dimension of the
%           data cell array.
%           >> setinds = { [12] [12] [12];... % 3 groups, cond 1
%                          [12] [12] [12] };  % 3 groups, cond 2
% allsubject - [cell array] all subjects (same order as in
%              STUDY.datasetinfo)
%
% Output:
%  data       - [cell array] data array with the subject or component selected
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
% 
% See also: std_erpplot(), std_specplot() and std_erspplot()

% Copyright (C) 2012 Arnaud Delorme
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

function [data] = std_selsubject(data, subject, setinds, allsubjects, optndims);

if nargin < 2
    help std_selsubject;
    return;
end

optndims = max(optndims, ndims(data{1}));
if isempty(strmatch(lower(subject), lower(allsubjects)))
    error(sprintf('Cannot select subject %s in list %s', subject, vararg2str({ allsubjects })));
end

% plot specific subject
% ---------------------
if size(setinds{1},1) > 1 && size(setinds{1},2) > 1 % single trials
    % possible subject indices
    selectInds = strmatch(lower(subject), lower(allsubjects),'exact');
    for c = 1:size(data,1)
        for g = 1:size(data,2)
            selectCol = [];
            for ind = 1:length(selectInds)
                selectCol = [ selectCol find(setinds{c,g}  == selectInds') ];
            end
            if optndims == 2
                data{c,g} = data{c,g}(:,selectCol); %2-D
            elseif optndims == 3
                data{c,g} = data{c,g}(:,:,selectCol); %3-D
            else
                data{c,g} = data{c,g}(:,:,:,selectCol); %4-D
            end
        end
    end
else
    for c = 1:size(data,1)
        for g = 1:size(data,2)
            subjectind = strmatch(lower(subject), lower(allsubjects),'exact');
            l = zeros(size(setinds{c,g}));
            for iSubj = 1:length(subjectind), l = l || setinds{c,g} == subjectind(iSubj); end
            if optndims == 2
                data{c,g}(:,~l) = []; %2-D
            elseif optndims == 3
                data{c,g}(:,:,~l) = []; %3-D
            else
                data{c,g}(:,:,:,~l) = []; %4-D
            end
        end
    end
end
