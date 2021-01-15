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
% Author: Arnaud Delorme, SCCN, 2018 based on a previous version of
%         Cyril Pernet (LIMO Team), The university of Edinburgh, 2014
%         Ramon Martinez-Cancino and Arnaud Delorme

% Copyright (C) 2018 Arnaud Delorm
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
    opt.method         = 'WSL'; % weighted least squares by default
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
        warndlg('std_limo error: Fieldtrip extension should be installed - chanlocs NOT generated');
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
if ~isempty(strfind(Analysis,'dat'))
    model.defaults.type = 'Channels';
elseif  ~isempty(strfind(Analysis,'ica'))
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
        warndlg2(sprintf('Components have not been clustered,\nLIMO will not match them across subjects'))
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
                warndlg2(neighbors_error.message,'limo_gp_level_chanlocs.mat not created')
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
    if isempty(contains(STUDY.filepath,'derivatives'))
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

% find out if the channels are interpolated
% -----------------------------------------
interpolated = zeros(1,length(STUDY.datasetinfo));
if strcmp(model.defaults.type,'Channels')
    for iDat = 1:length(STUDY.datasetinfo)
        fileName = fullfile(STUDY.datasetinfo(iDat).filepath, [ STUDY.datasetinfo(iDat).subject '.' opt.measure ]);
        tmpChans = load('-mat', fileName, 'labels');
        if length(tmpChans.labels) > ALLEEG(iDat).nbchan, interpolated(iDat) = 1; end
    end
end

% 1st level analysis
% -------------------------------------------------------------------------
model.cat_files  = [];
model.cont_files = [];
unique_subjects  = STUDY.design(STUDY.currentdesign).cases.value'; % unique_subject uses the STUDY name
nb_subjects      = length(unique_subjects);

% also should split per session
% for s = nb_subjects:-1:1
%     nb_sets(s) = numel(find(strcmp(unique_subjects{s},{STUDY.datasetinfo.subject})));
% end

% simply reshape to read columns
% -------------------------------------------------------------------------
order = cell(1,nb_subjects);
for s = 1:nb_subjects
    order{s} = find(strcmp(unique_subjects{s},{STUDY.datasetinfo.subject}));
end

% Cleaning old files from the current design (Cleaning ALL)
% -------------------------------------------------------------------------
% NOTE: Clean up the .lock files to (to be implemented)
% [STUDY.filepath filesep 'derivatives' filesep 'limo_batch_report']
if strcmp(opt.erase,'on')
    [~,filename] = fileparts(STUDY.filename);
    std_limoerase(STUDY.filepath, filename, unique_subjects, num2str(STUDY.currentdesign));
    STUDY.limo = [];
end

% Check if the measures has been computed
% -------------------------------------------------------------------------
for nsubj = 1 : nb_subjects
    inds     = find(strcmp(unique_subjects{nsubj},{STUDY.datasetinfo.subject}));

    % Checking for relative path
    study_fullpath = rel2fullpath(STUDY.filepath,STUDY.datasetinfo(inds(1)).filepath);
    %---
    subjpath = fullfile(study_fullpath, [unique_subjects{nsubj} '.' lower(Analysis)]);  % Check issue when relative path (remove comment)
    if ~exist(subjpath,'file')
        error('std_limo subject %s: Measures must be computed first',unique_subjects{nsubj});
    end
end
clear study_fullpath pathtmp;

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
FN = fieldnames(measureflags);
measureflags.(lower(opt.measureori))= 'on'; 
STUDY.etc.measureflags = measureflags;

% generate temporary merged datasets needed by LIMO
% -------------------------------------------------
fprintf('generating temporary files, pulling relevant trials ... \n')
mergedChanlocs = eeg_mergelocs(ALLEEG.chanlocs);
for s = 1:nb_subjects
    % field which are needed by LIMO
    % EEGLIMO.etc
    % EEGLIMO.times
    % EEGLIMO.chanlocs
    % EEGLIMO.srate
    % EEGLIMO.filepath
    % EEGLIMO.filename
    % EEGLIMO.icawinv
    % EEGLIMO.icaweights

    filename = [STUDY.datasetinfo(order{s}(1)).subject '_limo_file_tmp' num2str(design_index) '.set'];
    index    = [STUDY.datasetinfo(order{s}).index];
    tmp      = {STUDY.datasetinfo(order{s}).subject};
    if length(unique(tmp)) ~= 1
        error('it seems that sets of different subjects are merged')
    else
        names{s} =  cell2mat(unique(tmp));
    end

    % Creating fields for limo
    % ------------------------
    for sets = 1:length(index)
        EEGTMP = std_lm_seteegfields(STUDY,ALLEEG(index(sets)), index(sets),'datatype',model.defaults.type,'format', 'cell');
        ALLEEG = eeg_store(ALLEEG, EEGTMP, index(sets));
    end

    file_fullpath      = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).filepath);
    model.set_files{s} = fullfile(file_fullpath , filename);

    OUTEEG = [];
    if all([ALLEEG(index).trials] == 1)
         OUTEEG.trials = 1;
    else
        OUTEEG.trials = sum([ALLEEG(index).trials]);
    end

    filepath_tmp           = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).filepath);
    OUTEEG.filepath        = filepath_tmp;
    OUTEEG.filename        = filename;
    OUTEEG.srate           = ALLEEG(index(1)).srate;
    OUTEEG.icaweights      = ALLEEG(index(1)).icaweights;
    OUTEEG.icasphere       = ALLEEG(index(1)).icasphere;
    OUTEEG.icawinv         = ALLEEG(index(1)).icawinv;
    OUTEEG.icachansind     = ALLEEG(index(1)).icachansind;
    OUTEEG.etc             = ALLEEG(index(1)).etc;
    OUTEEG.times           = ALLEEG(index(1)).times;
    if any(interpolated)
        OUTEEG.chanlocs    = mergedChanlocs;
        OUTEEG.etc.interpolatedchannels = setdiff(1:length(OUTEEG.chanlocs), std_chaninds(OUTEEG, { ALLEEG(index(1)).chanlocs.labels }));
    else
        OUTEEG.chanlocs    = ALLEEG(index(1)).chanlocs;
    end

    % update EEG.etc
    OUTEEG.etc.merged{1}   = ALLEEG(index(1)).filename;

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
    % contains will not work in Octave
    single_trials_filename = fullfile(STUDY.datasetinfo(index(1)).filepath,  [STUDY.datasetinfo(index(1)).subject '.' FN{find(contains(FN,opt.measureori))}]);
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
end

% generate data files
% -------------------
fprintf('making up statistical models ... \n')
% by default we create a design matrix with all condition
factors = pop_listfactors(STUDY.design(opt.design), 'gui', 'off', 'level', 'one');
for s = 1:nb_subjects
    % save continuous and categorical data files
    trialinfo = std_combtrialinfo(STUDY.datasetinfo, unique_subjects{s});
    % [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, 'splitreg', opt.splitreg, 'interaction', opt.interaction);
    [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, 'splitreg', 'off', 'interaction', opt.interaction);
    if strcmpi(opt.splitreg,'on')
        for c=1:size(contMat,2)
            splitreg{c} = limo_split_continuous(catMat,contMat(:,c));
        end
        contMat    = cell2mat(splitreg);
        opt.zscore = 0; % regressors are now zscored
    end
    
    % copy results
    model.cat_files{s}                 = catMat;
    model.cont_files{s}                = contMat;
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
    STUDY.limo.subjects(s).subject     = unique_subjects{s};
    STUDY.limo.subjects(s).cat_file    = catMat;
    STUDY.limo.subjects(s).cont_file   = contMat;
end
 
% then we add contrasts for conditions that were merged during design selection
% i.e. multiple categorical variables (factors) and yet not matching the number 
% of variables (contrasts are then a weigthed sum of the crossed factors)
if ~isempty(factors) && length(STUDY.design(opt.design).variable) == 1 && isfield(factors, 'value') % only one non-continuous variable 
    if length(STUDY.design(opt.design).variable(1).value) ~= length(factors) % and this var has more values than the number of factors 
        limocontrast = zeros(length(STUDY.design(opt.design).variable(1).value),length(factors)+1); % length(factors)+1 to add the contant
        for n=1:length(factors)
            factor_names{n} = factors(n).value;
        end
        
        for c=1:length(STUDY.design(opt.design).variable.value)
            limocontrast(c,1:length(factors)) = single(ismember(factor_names,STUDY.design(opt.design).variable.value{c}));
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
    model.defaults.start    = ALLEEG(index(1)).xmin*1000;
    model.defaults.end      = ALLEEG(index(1)).xmax*1000;
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

elseif strcmp(Analysis,'dattimef') || strcmp(Analysis,'icaersp')
    model.defaults.analysis = 'Time-Frequency';
    model.defaults.start    = ALLEEG(index(1)).times(1);
    model.defaults.end      = ALLEEG(index(1)).times(end);
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
    [p,f,~]=fileparts(fullfile(STUDY.filepath,STUDY.filename));
    save(fullfile([p filesep 'LIMO_' f],[STUDY.design(opt.design).name '_contrast.mat']),'limocontrast');
end

STUDY.limo.model         = model;
STUDY.limo.datatype      = Analysis;
STUDY.limo.chanloc       = limoChanlocs;
if exist('limocontrast','var')
    STUDY.limo.contrast      = limocontrast;
end

% Save STUDY, and split saved file per groups
% -------------------------------------------
cd(STUDY.filepath);
STUDY      = pop_savestudy( STUDY, [],'filepath', STUDY.filepath,'savemode','resave');
keep_files = 'no';
if sum(procstatus) == nb_subjects
    disp('All subjects have been successfully processed.')
else
    if sum(procstatus)==0
        errordlg2('all subjects failed to process, check batch report')
    else
        warndlg2('some subjects failed to process, check batch report')
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

% split txt files if more than 1 group
if length(STUDY.group) > 1
    glm_name = [STUDY.design(STUDY.currentdesign).name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];
    for g= 1:length(STUDY.group)
        subset = arrayfun(@(x)(strcmpi(x.group,STUDY.group{g})), STUDY.datasetinfo);
        cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_Gp' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.mat(subset));
        cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_Gp' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.Beta(subset));
        if isfield(LIMO_files,'con')
            tmpcell = LIMO_files.con(subset);
            for c=1:length(tmpcell{1})
                cell2csv(fullfile(LIMO_files.LIMO, ['con' num2str(c) '_Gp' STUDY.group{g} '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
            end
        end
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
    
    if strfind(pathtmp(end),filesep)
        pathtmp = pathtmp(1:end-1); 
    end % Getting rid of filesep at the end
    
    if ~isempty(strfind(pathtmp(1:2),['.' filesep])) || ...
            (isunix && pathtmp(1) ~= '/') || (ispc && pathtmp(2) ~= ':')
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
