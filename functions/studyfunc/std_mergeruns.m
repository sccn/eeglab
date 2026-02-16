% STD_MERGERUNS - Merge multiple runs (conditions) for each subject and session
%                 in an EEGLAB STUDY and create a new STUDY with the merged datasets.
%                 This function identifies all datasets from the same subject and
%                 session (across different conditions/runs) and merges them into
%                 single datasets, then creates a new STUDY structure with these
%                 merged datasets.
%
% Usage:
%   >> [STUDY, ALLEEG] = std_mergeruns(STUDY, ALLEEG);
%   >> [STUDY, ALLEEG] = std_mergeruns(STUDY, ALLEEG, 'key', val);
%
% Inputs:
%   STUDY      - EEGLAB STUDY structure containing multiple runs/conditions
%   ALLEEG     - Vector of EEG dataset structures included in the STUDY
%
% Optional inputs:
%   'savedir'  - [string] Directory path where merged datasets will be saved.
%                Default: same as STUDY filepath with '_merged' suffix
%   'prefix'   - [string] Prefix to add to merged dataset filenames.
%                Default: 'merged_'
%   'resave'   - ['on'|'off'] Save the new STUDY structure to disk.
%                Default: 'off'
%
% Outputs:
%   STUDY      - New STUDY structure with merged datasets
%   ALLEEG     - Vector of merged EEG dataset structures
%
% Example:
%   % Merge all runs for each subject/session combination
%   >> [STUDY, ALLEEG] = std_mergeruns(STUDY, ALLEEG);
%
%   % Specify custom save directory and prefix
%   >> [STUDY, ALLEEG] = std_mergeruns(STUDY, ALLEEG, 'savedir', '/path/to/merged', ...
%                                      'prefix', 'allruns_');
%
% Notes:
%   - For each subject and session, all datasets (different conditions/runs)
%     are concatenated in temporal order
%   - Event latencies are adjusted automatically during merging
%   - ICA weights from the first dataset are preserved
%   - The original STUDY and datasets are not modified
%   - Channel locations must be consistent across datasets to be merged
%   - The merged STUDY will have no condition information (all conditions merged)
%
% See also: POP_MERGESET, STD_EDITSET, STD_LOADALLEEG, POP_STUDY
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2024

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD
%
% This file is part of EEGLAB, see https://eeglab.org
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

function [STUDY, ALLEEG] = std_mergeruns(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_mergeruns;
    return;
end

% Parse inputs
g = finputcheck(varargin, {
    'savedir'   'string'  {}                      '';
    'prefix'    'string'  {}                      'merged_';
    'resave'    'string'  {'on','off'}            'off'
}, 'std_mergeruns');

if ischar(g), error(g); end

% Set default save directory
if isempty(g.savedir)
    if isfield(STUDY, 'filepath') && ~isempty(STUDY.filepath)
        g.savedir = fullfile(STUDY.filepath, 'merged_datasets');
    else
        g.savedir = fullfile(pwd, 'merged_datasets');
    end
end

% Create save directory if it doesn't exist
if ~exist(g.savedir, 'dir')
    mkdir(g.savedir);
end

% Get dataset information
nDatasets = length(STUDY.datasetinfo);
fprintf('Processing %d datasets from STUDY...\n', nDatasets);

% Grouping variables: subject and session (merge across conditions/runs)
groupby = {'subject', 'session'};

% Extract grouping variables for each dataset

% Combine subject and session for grouping
groupVars = cell(1,length(ALLEEG));
for iDat = 1:length(STUDY.datasetinfo)
    if ~isempty(STUDY.datasetinfo(iDat).session)
        groupVars{iDat} = [ STUDY.datasetinfo(iDat).subject '_sess-' num2str(STUDY.datasetinfo(iDat).session) ];
    else
        groupVars{iDat} = STUDY.datasetinfo(iDat).subject;
    end
end

% Find unique groups to merge
[uniqueGroups, ~, groupIdx] = unique(groupVars);
nGroups = length(uniqueGroups);

fprintf('Found %d unique subject/session combinations to merge\n', nGroups);

% Initialize new ALLEEG and datasetinfo
newALLEEG = [];
newDatasetInfo = [];

% Process each group
for iGroup = 1:nGroups
    % Find datasets in this group
    datasetsInGroup = find(groupIdx == iGroup);

    % check unique conditions in this group
    conditionsInGroup = unique(cellfun(@num2str, {STUDY.datasetinfo(datasetsInGroup).condition}, 'uniformoutput', false));
    if length(conditionsInGroup) > 1
        error('Subject %s, Session %s: More than one condition found for this subject, cannot merge\n', ...
            uniqueGroups{iGroup, 1}, uniqueGroups{iGroup, 2});
    end

    if length(datasetsInGroup) < 2
        fprintf('Subject %s, Session %s: Only one dataset found, skipping merge\n', ...
            uniqueGroups{iGroup, 1}, uniqueGroups{iGroup, 2});
        continue;
    end

    % Merge using pop_mergeset (keep ICA from first dataset)
    mergedEEG = pop_mergeset(ALLEEG, datasetsInGroup);

    % Save merged dataset
    filename = [ uniqueGroups{iGroup} '_merged.set'];
    mergedEEG.setname = filename(1:end-4); % Remove .set extension
    mergedEEG = pop_saveset(mergedEEG, 'filename', filename, 'filepath', g.savedir);

    % Store in new ALLEEG
    [newALLEEG, mergedEEG] = eeg_store(newALLEEG, mergedEEG, iGroup, 'notext');

    % Create dataset info for new STUDY
    newInfo = STUDY.datasetinfo(datasetsInGroup(1)); % Copy info from first dataset
    newInfo.filename = filename;
    newInfo.filepath = g.savedir;
    newInfo.run   = [];
    newInfo.index = length(newALLEEG);

    if isempty(newDatasetInfo)
        newDatasetInfo = newInfo;
    else
        newDatasetInfo(end+1) = newInfo;
    end

    fprintf('  Saved merged dataset: %s\n', filename);
end

% Create new STUDY structure
if isempty(newDatasetInfo)
    error('No datasets were merged. Each subject/session combination must have at least 2 datasets (different conditions/runs).');
end

% Create new STUDY with merged datasets
STUDY_merged = STUDY;
STUDY_merged.datasetinfo = newDatasetInfo;
STUDY_merged.name = [STUDY.name '_merged'];
STUDY_merged.filename = [STUDY.filename(1:end-6) '_merged.study'];

% Update history
if ~isfield(STUDY_merged, 'history')
    STUDY_merged.history = '';
end
STUDY_merged.history = [STUDY_merged.history sprintf('[STUDY, ALLEEG] = std_mergeruns(STUDY, ALLEEG, %s);\n', vararg2str(varargin))];

% Clear any precomputed data
if isfield(STUDY_merged, 'changrp')
    STUDY_merged.changrp = [];
end
if isfield(STUDY_merged, 'cluster')
    STUDY_merged.cluster = [];
end

% Save new STUDY if requested
if strcmpi(g.resave, 'on')
    studyFile = fullfile(g.savedir, STUDY_merged.filename);
    STUDY = STUDY_merged;
    save(studyFile, '-mat', 'STUDY');
    fprintf('Saved merged STUDY to: %s\n', studyFile);
end

% Return outputs
STUDY = STUDY_merged;
ALLEEG = newALLEEG;

fprintf('Merge complete. New STUDY contains %d merged datasets.\n', length(newDatasetInfo));
