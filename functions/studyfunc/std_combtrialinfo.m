% std_combtrialinfo() - Combine and integrate all the info into the
%                       field trialinfo
%
% Usage:
%   >>   trialinfo = std_combtrialinfo(datasetinfo, [1 2]);
%   >>   trialinfo = std_combtrialinfo(datasetinfo, 'S01');
%
% Inputs:
%      datasetinfo         -  Datasetinfo structure from STUDY or
%                             datafiles
%      inds                - Either subject name  or indices of subject
%                            in datasetinfo structure.
%      trials              - All EEG dataset trial numbers
%
% Outputs:
%    trialinfo  - Updated trialinfo structure
%
% See also:
%   std_plotinfocluster
%
% Author: Arnaud Delorme
%         Ramon Martinez-Cancino
%
% Copyright (C) 2015  Arnaud Delorme, Ramon Martinez-Cancino,  INC, SCCN
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

function trialinfo = std_combtrialinfo(datasetinfo, inds, trials)

if nargin < 1
    help std_combtrialinfo;
    return;
end

if nargin < 3
    trials = ones(1, length(datasetinfo));
end

% Inds or subjectname
if isnumeric(inds)
    if ~isempty(find(~ismember(inds(:), [datasetinfo.index]))) %#ok<EFIND>
        error(['std_combtrialinfo error: Indices '' ' num2str(inds) ' '' aout of range']);
    end
elseif ischar(inds)
    if ~ismember(inds, {datasetinfo.subject})
        error(['std_combtrialinfo error: Subject name '' ' inds ' '' is invalid']);
    else
        % Looking for indices of subject 'ind'
        inds = find(cellfun(@(x) strcmp(x,inds), {datasetinfo.subject}));
    end
    
end

if ~isfield(datasetinfo, 'trialinfo')
    trialinfo = struct([]);
    nvals = [ 1 cumsum( [trials( [ datasetinfo(inds).index ]) ])+ 1 ];
else
    % check if duration field is present
    for iDat = inds(:)'
        if ~isfield(datasetinfo(iDat).trialinfo, 'duration') datasetinfo(iDat).trialinfo(1).duration = []; end
    end
    try
        trialinfo = [ datasetinfo(inds(:)').trialinfo ];
    catch
        warning( [ 'Trial information (event) field differ for subject ' datasetinfo(inds(1)).subject '...attempting to create missing event fields' ]);
        % creating missing fields
        allFields = {};
        for iInd = inds(:)'
            allFields = union(allFields, fieldnames(datasetinfo(iInd).trialinfo));
        end
        for iInd = inds(:)'
            for iField = 1:length(allFields)
                if ~isfield(datasetinfo(iInd).trialinfo, allFields{iField})
                    datasetinfo(iInd).trialinfo(1).(allFields{iField}) = [];
                end
            end
        end
        trialinfo = [ datasetinfo(inds(:)').trialinfo ];
    end
    nvals     = [ 1 cumsum(cellfun(@length, { datasetinfo(inds).trialinfo }))+1 ];
end

fields = fieldnames(datasetinfo);
fields = setdiff( fields, { 'filepath'  'filename' 'subject' 'index' 'comps' 'trialinfo' });
for iDat = 1:length(inds)
    
    for iField = 1:length(fields)
        [trialinfo(nvals(iDat):nvals(iDat+1)-1).(fields{iField})] = deal( datasetinfo(inds(iDat)).(fields{iField}) );
    end
    
    % Checking if field do not exist or if exist and is empty
    % This is old legacy code and it is unclear why the test were so
    % complex so we left it for now to see if the new code generates any
    % problem
    %     if ~isfield(trialinfo,'condition') || length(trialinfo) < nvals(iDat+1)-1 || any(cellfun(@isempty, {trialinfo(nvals(iDat):nvals(iDat+1)-1).condition}))
    %         [trialinfo(nvals(iDat):nvals(iDat+1)-1).condition] = deal( datasetinfo(inds(iDat)).condition );
    %     end
    %
    %     if ~isfield(trialinfo,'group') || any(cellfun(@isempty, {trialinfo(nvals(iDat):nvals(iDat+1)-1).group}))
    %         [trialinfo(nvals(iDat):nvals(iDat+1)-1).group    ] = deal( datasetinfo(inds(iDat)).group     );
    %     end
    %
    %     if ~isfield(trialinfo,'session') || any(cellfun(@isempty, {trialinfo(nvals(iDat):nvals(iDat+1)-1).session}))
    %         [trialinfo(nvals(iDat):nvals(iDat+1)-1).session  ] = deal( datasetinfo(inds(iDat)).session   );
    %     end
    
end
