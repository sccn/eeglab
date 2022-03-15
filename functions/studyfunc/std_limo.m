
% std_limo() - Export and run in LIMO the EEGLAB STUDY design.
%           call limo_batch to create all 1st level LIMO_EEG analysis + RFX
%
% Usage:
%   [STUDY LIMO_files] = std_limo(STUDY,ALLEEG,'key',val)
%
% Inputs:
%  STUDY        - studyset structure containing some or all files in ALLEEG
%  ALLEEG       - vector of loaded EEG datasets
%
% Optional inputs:
%  'measure'      - ['daterp'|'icaerp'|'datspec'|'icaspec'|'datersp'|'icaersp']
%                   measure to compute. Currently, only 'daterp' and
%                   'datspec' are supported. Default is 'daterp'.
%  'method'       - ['OLS'|'WTS'|'IRLS'] Ordinary Least Squares (OLS) or Weighted
%                   Least Squares (WTS) or Iterative Reweighted Least Squares'IRLS'.
%                   WTS should be used as it is more robust. IRLS is much slower
%                   and better across subjects than across trials.
%  'design'       - [integer] design index to process. Default is the current
%                   design stored in STUDY.currentdesign.
%  'erase'        - ['on'|'off'] erase previous files. Default is 'on'.
%  'neighboropt'  - [cell] cell array of options for the function computing
%                   the channel neighbox matrix std_prepare_chanlocs(). The file
%                   is saved automatically if channel location are present.
%                   This option allows to overwrite the defaults when computing
%                   the channel neighbox matrix.
%   'chanloc'     - Channel location structure. Must be used with 'neighbormat',
%                   or it will be ignored. If this option is used, it will
%                   ignore 'neighboropt' if used.
%   'neighbormat' - Neighborhood matrix of electrodes. Must be used with 'chanloc',
%                   or it will be ignored. If this option is used, it will
%                   ignore 'neighboropt' if used.
%   'freqlim'     - Frequency trimming in Hz
%   'timelim'     - Time trimming in millisecond
%
% Outputs:
%  STUDY     - modified STUDY structure (the STUDY.design now contains a list
%              of the limo files)
%  LIMO_files a structure with the following fields
%     LIMO_files.LIMO the LIMO folder name where the study is analyzed
%     LIMO_files.mat a list of 1st level LIMO.mat (with path)
%     LIMO_files.Beta a list of 1st level Betas.mat (with path)
%     LIMO_files.con a list of 1st level con.mat (with path)
%     LIMO_files.expected_chanlocs expected channel location neighbor file for
%                                  correcting for multiple comparisons
% Example:
%  [STUDY LIMO_files] = std_limo(STUDY,ALLEEG,'measure','daterp')
%
% Author: Arnaud Delorme (SCCN) & Cyril Pernet (LIMO Team)
%         Based on previous version from Ramon Martinez-Cancino and Arnaud Delorme
%
% Copyright (C) 2018 Arnaud Delorme
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

function [STUDY, LIMO_files] = std_limo(STUDY,ALLEEG,varargin)

LIMO_files = [];

if isempty(STUDY.filepath)
    STUDY.filepath = pwd;
end
cd(STUDY.filepath);
if nargin < 2
    help std_limo;
    return;
end

warning('off', 'MATLAB:lang:cannotClearExecutingFunction');
if ischar(varargin{1}) && ( strcmpi(varargin{1}, 'daterp') || ...
        strcmpi(varargin{1}, 'datspec') || ...
        strcmpi(varargin{1}, 'dattimef') || ...
        strcmpi(varargin{1}, 'icaerp')|| ...
        strcmpi(varargin{1}, 'icaspec')|| ...
        strcmpi(varargin{1}, 'icatimef'))
    opt.measure  = varargin{1};
    opt.design   = varargin{2};
    opt.ow_chanlocfile = 'no';  % if chanloc file exist, do not overwrite
    opt.erase          = 'on';  % erase previous folders/file with the same name
    opt.method         = 'WLS'; % weighted least squares by default
    opt.zscore         = 1;     % zscore regressors
else
    opt = finputcheck( varargin, ...
        { 'measure'        'string'  { 'daterp' 'datspec' 'dattimef' 'icaerp' 'icaspec' 'icatimef' } 'daterp'; ...
        'method'         'string'  { 'OLS' 'WLS' 'IRLS' } 'WLS';
        'design'         'integer' [] STUDY.currentdesign;
        'erase'          'string'  { 'on','off' }   'off';
        'splitreg'       'string'  { 'on','off' }   'off';
        'interaction'    'string'  { 'on','off' }   'off';
        'freqlim'        'real'    []               [] ;
        'timelim'        'real'    []               [] ;
        'neighboropt'    'cell'    {}               {} ;
        'chanloc'        'struct'  {}               struct('no', {}); % default empty structure
        'neighbormat'    'real'    []               [] ;
        'zscore'         'real'    [0,1]            1  ;
        'ow_chanlocfile' 'string'  {'yes','no'}     'no'},...
        'std_limo');
    if ischar(opt), error(opt); end
end
opt.measureori = opt.measure;
if strcmpi(opt.measure, 'datersp')
    opt.measure = 'dattimef';
end

Analysis     = opt.measure;
design_index = opt.design;

% Make sure paths are ok for LIMO (Consider to move this to eeglab.m in a future)
% -------------------------------------------------------------------------
root = fileparts(which('limo_eeg'));
addpath([root filesep 'limo_cluster_functions']);
addpath([root filesep 'external' filesep 'psom']);
addpath([root filesep 'external']);
addpath([root filesep 'help']);

% Checking fieldtrip paths to compute gp channel location
skip_chanlocs   = 0;
chanloc_created = 0;
limoChanlocs    = [];
if exist(fullfile([STUDY.filepath filesep 'derivatives'], 'limo_gp_level_chanlocs.mat'),'file')
    limoChanlocs = fullfile([STUDY.filepath filesep 'derivatives'], 'limo_gp_level_chanlocs.mat');
elseif exist(fullfile(STUDY.filepath, 'limo_gp_level_chanlocs.mat'),'file')
    limoChanlocs = fullfile(STUDY.filepath, 'limo_gp_level_chanlocs.mat');
elseif exist(fullfile([STUDY.filepath filesep 'derivatives'], 'limo_chanlocs.mat'),'file')
    limoChanlocs = fullfile([STUDY.filepath filesep 'derivatives'], 'limo_chanlocs.mat');
elseif exist(fullfile(STUDY.filepath, 'limo_chanlocs.mat'),'file')
    limoChanlocs = fullfile(STUDY.filepath, 'limo_chanlocs.mat');
end

if ~isempty(limoChanlocs)
    if ~strcmpi(opt.ow_chanlocfile,'no') % empty or yes
        opt.ow_chanlocfile = questdlg2('channel location file found, do you want to overwrite','overwrite?','yes','no','no');
    end

    if isempty(opt.ow_chanlocfile) || strcmpi(opt.ow_chanlocfile,'no')
        skip_chanlocs = 1;
    end
end

if skip_chanlocs == 0
    if ~exist('ft_prepare_neighbours','file')
        warndlg('std_limo error: Fieldtrip extension should be installed - chanlocs NOT generated', '', 'non-modal');
        skip_chanlocs = 1;
    else
        if ~exist('eeglab2fieldtrip','file')
            root = fileparts(which('ft_prepare_neighbours'));
            addpath([root filesep 'external' filesep 'eeglab']);
        end
    end
end

% Detecting type of analysis
% -------------------------------------------------------------------------
model.defaults.datatype = opt.measureori(4:end);
if ~isempty(strfind(Analysis,'dat')) %#ok<STREMP>
    model.defaults.type = 'Channels';
elseif  ~isempty(strfind(Analysis,'ica')) %#ok<STREMP>
    [STUDY,flags]=std_checkdatasession(STUDY,ALLEEG);
    if sum(flags)>0
        error('some subjects have data from different sessions - can''t do ICA');
    end
    model.defaults.type = 'Components';
end

% Checking if clusters
% -------------------------------------------------------------------------
if strcmp(model.defaults.type,'Components')
    if isempty(STUDY.cluster(1).child)
        warndlg2(sprintf('Components have not been clustered,\nLIMO will not match them across subjects'), '', 'non-modal')
        model.defaults.icaclustering = 0;
    else
        model.defaults.icaclustering = 1;
    end
end

% computing channel neighbour matrix
% ---------------------------------
if skip_chanlocs == 0
    chanloc_created = 1;
    if isempty(opt.chanloc) && isempty(opt.neighbormat)
        if isfield(ALLEEG(1).chanlocs, 'theta') &&  ~strcmp(model.defaults.type,'Components')
            if  ~isfield(STUDY.etc,'statistic')
                STUDY = pop_statparams(STUDY, 'default');
            end

            try
                [~,~,limoChanlocs] = std_prepare_neighbors(STUDY, ALLEEG, 'force', 'on', opt.neighboropt{:});
                chanlocname = 'limo_gp_level_chanlocs.mat';
            catch neighbors_error
                limoChanlocs = []; chanloc_created = 0;
                warndlg2(neighbors_error.message,'limo_gp_level_chanlocs.mat not created', 'non-modal')
            end
        else
            limoChanlocs = []; chanloc_created = 0;
            if ~isempty(STUDY.cluster(1).child)
                disp('Warning: cannot compute expected channel distance for correction for multiple comparisons');
            end
        end
    else
        limoChanlocs.expected_chanlocs   = opt.chanloc;
        limoChanlocs.channeighbstructmat = opt.neighbormat;
        chanlocname = 'limo_chanlocs.mat';
    end
end

if chanloc_created
    % contains will not work in Octave
    if isempty(strfind(STUDY.filepath,'derivatives'))
        if ~exist([STUDY.filepath filesep 'derivatives'],'dir')
            mkdir([STUDY.filepath filesep 'derivatives']);
        end
        limoChanlocsFile = fullfile([STUDY.filepath filesep 'derivatives'], chanlocname);
    else
        limoChanlocsFile = fullfile(STUDY.filepath, chanlocname);
    end
    save('-mat', limoChanlocsFile, '-struct', 'limoChanlocs');
    fprintf('Saving channel neighbors for correction for multiple comparisons in \n%s\n', limoChanlocsFile);
end

% 1st level analysis
% -------------------------------------------------------------------------
model.cat_files  = [];
model.cont_files = [];

% Cleaning old files from the current design (Cleaning ALL)
% -------------------------------------------------------------------------
% NOTE: Clean up the .lock files to (to be implemented)
% [STUDY.filepath filesep 'derivatives' filesep 'limo_batch_report']
if strcmp(opt.erase,'on')
    [~,filename] = fileparts(STUDY.filename);
    std_limoerase(STUDY.filepath, filename, STUDY.subject, num2str(STUDY.currentdesign));
    STUDY.limo = [];
end

% Check if the measures has been computed
% also find out if the channels are interpolated
% -------------------------------------------------------------------------
interpolated = zeros(1,length(STUDY.datasetinfo));
for iDat = 1:length(STUDY.datasetinfo)
    fileName = fullfile(STUDY.datasetinfo(iDat).filepath, [ STUDY.datasetinfo(iDat).subject '*.' opt.measure ]);
    % fileName should already match unless user moves / rename, hence using dir
    fileName = dir(fileName);
    if isempty(fileName)
        error('std_limo subject %s: Measures must be computed first',STUDY.datasetinfo(iDat).subject);
    else
        if strcmp(model.defaults.type,'Channels')
            tmpChans = load('-mat', fullfile(fileName(1).folder,fileName(1).name), 'labels');
            if length(tmpChans.labels) > ALLEEG(iDat).nbchan, interpolated(iDat) = 1; end
        end
    end
end

measureflags = struct('daterp','off',...
    'datspec','off',...
    'datersp','off',...
    'dattimef','off',...
    'datitc' ,'off',...
    'icaerp' ,'off',...
    'icaspec','off',...
    'icatimef','off',...
    'icaersp','off',...
    'icaitc','off');
measureflags.(lower(opt.measureori))= 'on';
STUDY.etc.measureflags = measureflags;
mergedChanlocs = eeg_mergelocs(ALLEEG.chanlocs);
fprintf('generating temporary files, pulling relevant trials ... \n')

% generate temporary merged datasets needed by LIMO
% -------------------------------------------------
allSubjects    = { STUDY.datasetinfo.subject };
allSessions    = { STUDY.datasetinfo.session };
uniqueSubjects = unique(allSubjects);
nb_subjects    = length(uniqueSubjects);
allSessions(cellfun(@isempty, allSessions)) = { 1 };
allSessions    = cellfun(@num2str, allSessions, 'uniformoutput', false);
uniqueSessions = unique(allSessions);

% by default we create a design matrix with all condition
factors = pop_listfactors(STUDY.design(opt.design), 'gui', 'off', 'level', 'one');

for iSubj = 1:nb_subjects
    for iSess = 1:length(uniqueSessions)
        inds1 = strmatch( uniqueSubjects{iSubj}, allSubjects, 'exact');
        inds2 = strmatch( uniqueSessions{iSess}, allSessions, 'exact');
        inds  = intersect(inds1, inds2);
        if ~isempty(inds)
            if length(inds) ~= 1
                error([ 'Cannot calculate contrast because more than 1 dataset per session' 10 ...
                    'per subject. Merge datasets for each subject and try again.' ]);
            end
            
            % make file-up
            [~,subname] = fileparts(STUDY.datasetinfo(inds).filename);
            if isfield(ALLEEG,'filename')
                if ~strcmp(subname,ALLEEG(inds).filename(1:end-4))
                    error('STUDY and ALLEEG mismatch, can''t figure out which file to use')
                end
            else
                warning('No filename in ALLEEG, pulling data blindly from STUDY')
            end

            if strcmp(subname(1:4),'sub-')
                if contains(subname,'ses-')
                    filename = [subname '_design' num2str(design_index) '.set'];
                else
                    filename = [subname '_ses-' num2str(iSess) '_design' num2str(design_index) '.set'];
                end
            else
                if contains(subname,'ses-')
                    filename = ['sub-' subname '_design' num2str(design_index) '.set'];
                else
                    filename = ['sub-' subname '_ses-' num2str(iSess) '_design' num2str(design_index) '.set'];
                end
            end

            % Creating fields for limo
            % ------------------------
            fprintf('pulling trials for %s ... \n',filename)
            EEGTMP                     = std_lm_seteegfields(STUDY,ALLEEG(inds), inds,'datatype',model.defaults.type,'format', 'cell');
            ALLEEG                     = eeg_store(ALLEEG, EEGTMP, inds);
            file_fullpath              = rel2fullpath(STUDY.filepath,ALLEEG(inds).filepath);
            model.set_files{inds}      = fullfile(file_fullpath , filename);

            OUTEEG = [];
            if all([ALLEEG(inds).trials] == 1)
                OUTEEG.trials = 1;
            else
                OUTEEG.trials = sum([ALLEEG(inds).trials]);
            end

            filepath_tmp               = rel2fullpath(STUDY.filepath,ALLEEG(inds).filepath);
            OUTEEG.filepath            = filepath_tmp;
            OUTEEG.filename            = filename;
            OUTEEG.srate               = ALLEEG(inds).srate;
            OUTEEG.icaweights          = ALLEEG(inds).icaweights;
            OUTEEG.icasphere           = ALLEEG(inds).icasphere;
            OUTEEG.icawinv             = ALLEEG(inds).icawinv;
            OUTEEG.icachansind         = ALLEEG(inds).icachansind;
            OUTEEG.etc                 = ALLEEG(inds).etc;
            OUTEEG.times               = ALLEEG(inds).times;
            if any(interpolated)
                OUTEEG.chanlocs        = mergedChanlocs;
                OUTEEG.etc.interpolatedchannels = setdiff(1:length(OUTEEG.chanlocs), std_chaninds(OUTEEG, { ALLEEG(inds).chanlocs.labels }));
            else
                OUTEEG.chanlocs        = ALLEEG(inds).chanlocs;
            end

            % update EEG.etc
            OUTEEG.etc.merged{1}       = ALLEEG(inds).filename;

            % Def fields
            OUTEEG.etc.datafiles.daterp   = [];
            OUTEEG.etc.datafiles.datspec  = [];
            OUTEEG.etc.datafiles.datersp  = [];
            OUTEEG.etc.datafiles.dattimef = [];
            OUTEEG.etc.datafiles.datitc   = [];
            OUTEEG.etc.datafiles.icaerp   = [];
            OUTEEG.etc.datafiles.icaspec  = [];
            OUTEEG.etc.datafiles.icaersp  = [];
            OUTEEG.etc.datafiles.icatimef = [];
            OUTEEG.etc.datafiles.icaitc   = [];

            % Filling fields
            single_trials_filename = EEGTMP.etc.datafiles.(opt.measureori);
            if exist(single_trials_filename,'file')
                if strcmpi(measureflags.daterp,'on')
                    OUTEEG.etc.datafiles.daterp = single_trials_filename;
                elseif strcmpi(measureflags.datspec,'on')
                    OUTEEG.etc.datafiles.datspec = single_trials_filename;
                elseif strcmpi(measureflags.datersp,'on')
                    OUTEEG.etc.datafiles.datersp = single_trials_filename;
                elseif strcmpi(measureflags.datitc,'on')
                    OUTEEG.etc.datafiles.datitc = single_trials_filename;
                elseif strcmpi(measureflags.icaerp,'on')
                    OUTEEG.etc.datafiles.icaerp = single_trials_filename;
                elseif strcmpi(measureflags.icaspec,'on')
                    OUTEEG.etc.datafiles.icaspec = single_trials_filename;
                elseif strcmpi(measureflags.icaersp,'on')
                    OUTEEG.etc.datafiles.icaersp = single_trials_filename;
                elseif strcmpi(measureflags.icaitc,'on')
                    OUTEEG.etc.datafiles.icaitc = single_trials_filename;
                elseif strcmpi(measureflags.dattimef,'on')
                    OUTEEG.etc.datafiles.dattimef = single_trials_filename;
                end
            end

            % Save info
            EEG = OUTEEG;
            save('-mat', fullfile( filepath_tmp, OUTEEG.filename), 'EEG');
            clear OUTEEG filepath_tmp

            % generate data files
            % -------------------
            fprintf('making up statistical model for %s ... \n',filename)
            % save continuous and categorical data files
            trialinfo = std_combtrialinfo(STUDY.datasetinfo, inds);
            % [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, 'splitreg', opt.splitreg, 'interaction', opt.interaction);
            [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, 'splitreg', 'off', 'interaction', opt.interaction);
            if strcmpi(opt.splitreg,'on')
                for c=size(contMat,2):-1:1
                    splitreg{c} = limo_split_continuous(catMat,contMat(:,c)); % std_limodesign does something else when splitting regressors
                end
                contMat    = cell2mat(splitreg);
                opt.zscore = 0; % regressors are now zscored
            end

            % copy results
            model.cat_files{inds}                 = catMat;
            model.cont_files{inds}                = contMat;
            if isfield(limodesign, 'categorical')
                STUDY.limo.categorical = limodesign.categorical;
            else
                STUDY.limo.categorical = {};
            end
            if isfield(limodesign, 'continuous')
                STUDY.limo.continuous = limodesign.continuous;
            else
                STUDY.limo.continuous = {};
            end
            STUDY.limo.subjects(inds).subject     = STUDY.datasetinfo(inds(1)).subject;
            STUDY.limo.subjects(inds).cat_file    = catMat;
            STUDY.limo.subjects(inds).cont_file   = contMat;
        end
    end % exit session
end % exit subject

% then we add contrasts for conditions that were merged during design selection
% i.e. multiple categorical variables (factors) and yet not matching the number
% of variables (contrasts are then a weighted sum of the crossed factors)
if ~isempty(factors) && isfield(factors, 'value') && ...
        sum(arrayfun(@(x) ~strcmpi(x.label,'group'),STUDY.design(opt.design).variable)) == 1 % only one non-continuous variable other than group
    if length(STUDY.design(opt.design).variable(1).value) ~= length(factors) % and this var has more values than the number of factors
        limocontrast = zeros(length(STUDY.design(opt.design).variable(1).value),length(factors)+1); % length(factors)+1 to add the constant
        for n=length(factors):-1:1
            factor_names{n} = factors(n).value;
        end

        index = find(arrayfun(@(x) ~strcmpi(x.label,'group'),STUDY.design(opt.design).variable)); % which one is not group
        for c=1:length(STUDY.design(opt.design).variable(index).value)
            limocontrast(c,1:length(factors)) = single(ismember(factor_names,STUDY.design(opt.design).variable(index).value{c}));
            limocontrast(c,1:length(factors)) = limocontrast(c,1:length(factors)) ./ sum(limocontrast(c,1:length(factors))); % scale by the number of variables
        end
    end
end

% transpose
model.set_files  = model.set_files';
model.cat_files  = model.cat_files';
model.cont_files = model.cont_files';
if all(cellfun(@isempty, model.cat_files )), model.cat_files  = []; end
if all(cellfun(@isempty, model.cont_files)), model.cont_files = []; end


% set model.defaults - all conditions no bootstrap
% -----------------------------------------------------------------
% to update passing the timing/frequency from STUDY - when computing measures
% -----------------------------------------------------------------
if strcmp(Analysis,'daterp') || strcmp(Analysis,'icaerp')
    model.defaults.analysis = 'Time';
    for s=nb_subjects:-1:1
        vs(s) = ALLEEG(s).xmin*1000;
        ve(s) = ALLEEG(s).xmax*1000;
    end
    model.defaults.start    = max(vs);
    model.defaults.end      = min(ve);

    if length(opt.timelim) == 2 && opt.timelim(1) < opt.timelim(end)
        % start value
        if opt.timelim(1) < model.defaults.start
            fprintf('std_limo: Invalid time lower limit, using default value instead');
        else
            model.defaults.start = opt.timelim(1);
        end
        % end value
        if opt.timelim(end) > model.defaults.end
            fprintf('std_limo: Invalid time upper limit, using default value instead');
        else
            model.defaults.end = opt.timelim(end);
        end
    end
    model.defaults.lowf  = [];
    model.defaults.highf = [];

elseif strcmp(Analysis,'datspec') || strcmp(Analysis,'icaspec')

    model.defaults.analysis= 'Frequency';
    if length(opt.freqlim) == 2
        model.defaults.lowf    = opt.freqlim(1);
        model.defaults.highf   = opt.freqlim(2);
    else
        error('std_limo: Frequency limits need to be specified');
    end
    model.defaults.start    = [];
    model.defaults.end      = [];

elseif strcmp(Analysis,'dattimef') || strcmp(Analysis,'icaersp')
    model.defaults.analysis = 'Time-Frequency';
    for s=nb_subjects:-1:1
        vs(s) = ALLEEG(s).xmin*1000;
        ve(s) = ALLEEG(s).xmax*1000;
    end
    model.defaults.start    = max(vs);
    model.defaults.end      = min(ve);
    model.defaults.lowf     = [];
    model.defaults.highf    = [];

    if length(opt.timelim) == 2
        model.defaults.start    = opt.timelim(1);
        model.defaults.end      = opt.timelim(2);
    end
    if length(opt.freqlim) == 2
        model.defaults.lowf     = opt.freqlim(1);
        model.defaults.highf    = opt.freqlim(2);
    else
        error('std_limo: Frequency limits need to be specified');
    end
end

model.defaults.fullfactorial    = 0;                 % all variables
model.defaults.zscore           = opt.zscore;        % done that already
model.defaults.bootstrap        = 0 ;                % only for single subject analyses - not included for studies
model.defaults.tfce             = 0;                 % only for single subject analyses - not included for studies
model.defaults.method           = opt.method;        % default is OLS - to be updated to 'WLS' once validated
model.defaults.Level            = 1;                 % 1st level analysis
model.defaults.type_of_analysis = 'Mass-univariate'; % option can be multivariate (work in progress)


if ~exist('limocontrast','var')
    [LIMO_files, procstatus] = limo_batch('model specification',model,[],STUDY);
else
    contrast.mat = limocontrast;
    [LIMO_files, procstatus] = limo_batch('both',model,contrast,STUDY);
    if exist(fullfile([STUDY.filepath filesep 'derivatives']),'dir')
 %       save(fullfile([STUDY.filepath filesep 'derivatives']),[STUDY.design(opt.design).name '_contrast.mat']),'limocontrast'));
        save(fullfile([STUDY.filepath filesep 'derivatives'],[STUDY.design(opt.design).name '_contrast.mat']),'limocontrast');
    else
        save(fullfile(STUDY.filepath,[STUDY.design(opt.design).name '_contrast.mat']),'limocontrast');
    end
end

STUDY.limo.model         = model;
STUDY.limo.datatype      = Analysis;
STUDY.limo.chanloc       = limoChanlocs;
if exist('limocontrast','var')
    STUDY.limo.contrast      = limocontrast;
end

% generate between session contrasts
% ----------------------------------
index = 1;
for s = 1:nb_subjects
   sess_index = find(cellfun(@(x) strcmpi(x,uniqueSubjects{s}), allSubjects));
   % matches sess_index = find(contains(LIMO_files.mat,uniqueSubjects{s}))
   if length(sess_index) > 1
        fprintf('std_limo, computing additional between sessions contrasts for subject %s\n',uniqueSubjects{s})
        sess_name = allSessions(sess_index);
        pairs = nchoosek(1:length(sess_index),2); % do all session pairs
        parfor p=1:size(pairs,1)
            strpair = [cell2mat(sess_name(pairs(p,1))) cell2mat(sess_name(pairs(p,2)))];
            strpair(isspace(strpair)) = []; % remove spaces
            filesout{p} = limo_contrast_sessions(cell2mat(LIMO_files.mat(sess_index(pairs(p,1)))), ...
                cell2mat(LIMO_files.mat(sess_index(pairs(p,2)))),strpair);
        end
        
        for f=1:length(filesout)
            for ff=1:length(filesout{f})
                allcon{index} = filesout{f}{ff};
                index = index +1;
            end
        end
        clear filesout
    end
end

% use same glm_name as limo_batch
design_name = STUDY.design(STUDY.currentdesign).name;
design_name(isspace(design_name)) = [];
if findstr(design_name,'STUDY.')
    design_name = design_name(7:end);
end
glm_name = [STUDY.filename(1:end-6) '_' design_name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];

% further split that list per regressor and group
if exist('allcon','var')
    maxcon = max(cellfun(@(x) str2double(x(strfind(x,'con_')+4:strfind(x,'sess_')-1)),allcon));
    for con=1:maxcon
        index = find(cellfun(@(x) ~isempty(x),cellfun(@(x) strfind(x,['con_' num2str(con)]),allcon','UniformOutput',false)));
        cell2csv([LIMO_files.LIMO filesep 'Between_sessions_con_' num2str(con) '_' glm_name '.txt'], allcon(index)');
        if length(STUDY.group) > 1
            for g= 1:length(STUDY.group)
                % find subjects of group g
                subset = find(arrayfun(@(x)(strcmpi(x.group,STUDY.group{g})), STUDY.datasetinfo));
                for s=1:length(subset)
                    sub{s} = STUDY.datasetinfo(subset(s)).subject;
                end
                sub = unique(sub);
                % find subjects of group g and contrast con
                subcon = [];
                for s = 1:length(sub)
                    if strfind(sub{s},'sub-')
                        subindex = find(cellfun(@(x) ~isempty(x),(cellfun(@(x) strfind(x,sub{s}),allcon','UniformOutput',false)))); % subject s group g
                    else
                        subindex = find(cellfun(@(x) ~isempty(x),(cellfun(@(x) strfind(x,['sub-' sub{s} ]),allcon','UniformOutput',false)))); % subject s group g
                    end
                    subcon = [subcon;intersect(index,subindex)];
                end
                % save
                if ~isempty(subcon)
                    cell2csv([LIMO_files.LIMO filesep 'Between_sessions_con_' num2str(con) 'Gp' STUDY.group{g} '_' glm_name '.txt'], allcon(subcon)');
                end
            end
        end
    end
end

% Save STUDY - delete tmp files
% ------------------------------
cd(STUDY.filepath);
STUDY      = pop_savestudy( STUDY, [],'filepath', STUDY.filepath,'savemode','resave');
keep_files = 'no';
if all(procstatus)
    disp('All subjects have been successfully processed.')
else
    if sum(procstatus)==0 % not a WLS issue - limo_batch errors for that and tells the user
        errordlg2('all subjects failed to process, check limo batch report')
    else
        warndlg2('some subjects failed to process, check limo batch report','', 'non-modal')
    end
    % cleanup temp files - except for subjects with errors?
    keep_files = questdlg('Do you want to keep temp files of unsuccessulfully processed subjects','option for manual debugging','yes','no','no');
end

% delete
if isempty(keep_files) || strcmpi(keep_files,'no')
    for s = 1:nb_subjects
        delete(model.set_files{s});
    end
else
    for s = find(procstatus)
        delete(model.set_files{s});
    end
end

% -------------------------------------------------------------------------
% Return full path if 'filepath' is a relative path. The output format will
% fit the one of 'filepath'. That means that if 'filepath' is a cell array,
% then the output will a cell array too, and the same if is a string.
function file_fullpath = rel2fullpath(studypath,filepath)

nit = 1; if iscell(filepath), nit = length(filepath);end

for i = 1:nit
    if iscell(filepath)
        pathtmp = filepath{i};
    else
        pathtmp = filepath;
    end

    if strfind(pathtmp(end),filesep) %#ok<STRIFCND>
        pathtmp = pathtmp(1:end-1);
    end % Getting rid of filesep at the end

    if ~isempty(strfind(pathtmp(1:2),['.' filesep])) || ...
            (isunix && pathtmp(1) ~= '/') || (ispc && pathtmp(2) ~= ':') %#ok<STREMP>
        if iscell(filepath)
            file_fullpath{i} = fullfile(studypath,pathtmp(1:end));
        else
            file_fullpath = fullfile(studypath,pathtmp(1:end));
        end
    else
        if iscell(filepath)
            file_fullpath{i} = pathtmp;
        else
            file_fullpath = pathtmp;
        end
    end
end

% get file base name (from std_precomp)
% ------------------
function filebase = getfilename(filepath, subj, sess, fileSuffix, onlyOneSession)
if onlyOneSession
    filebase = fullfile(filepath, [ subj fileSuffix ] );
else
    sesStr   = [ '0' sess ];
    filebase = fullfile(filepath, [ subj '_ses-' sesStr(end-1:end) fileSuffix ] );
end
