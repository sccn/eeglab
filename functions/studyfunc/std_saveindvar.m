% std_saveindvar - save independent variables of a STUDY
%
% Usage:
%   std_saveindvar(STUDY, indvars, filename);
%
% Input:
%   STUDY - EEGLAB STUDY structure
%   indvars - [cell] cell array of independent variables
%   filename - [string] file name to save variables in
%
% Authors: A. Delorme, CERCO/CNRS and SCCN/UCSD, 2022-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2010, arno@ucsd.edu
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

function std_saveindvar(STUDY, indvars, filename)

if nargin < 1
    help std_saveindvar;
    return;
end

fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing');
end
for iSubj = 1:length(STUDY.subject)
    for iVar  = 1:length(indvars)
        indS = strmatch(STUDY.subject{iSubj}, {STUDY.datasetinfo.subject}, 'exact');
        allVals = [ STUDY.datasetinfo(indS).(indvars{iVar})];
        if isempty(allVals), allVals = NaN; end
        if length(unique(allVals)) > 1
            fprintf('Warning: more than 1 value for "%s" for subject %s, taking the average\n', indvars{iVar}, STUDY.subject{iSubj});
            allVals = mean(allVals);
        end
        fprintf(fid, '%f\t', allVals(1));
    end
    fprintf(fid, '\n');
end
fclose(fid);

