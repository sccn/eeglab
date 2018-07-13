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
%  'method'       - ['OLS'|'WTS'] Ordinary Least Square (OLS) or Weighted Least
%                   Square (WTS). WTS should be used as it is more robust. It is
%                   slower though.
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
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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
if isstr(varargin{1}) && ( strcmpi(varargin{1}, 'daterp') || strcmpi(varargin{1}, 'datspec') || strcmpi(varargin{1}, 'icaerp')|| strcmpi(varargin{1}, 'icaspec'))
    opt.measure  = varargin{1};
    opt.design   = varargin{2};
    opt.erase    = 'on';
    opt.method   = 'OSL';
else
    opt = finputcheck( varargin, ...
        { 'measure'        'string'  { 'daterp' 'datspec' 'icaerp' 'icaspec'} 'daterp'; ...
          'method'         'string'  { 'OLS' 'WLS' } 'OLS';
          'design'         'integer' [] STUDY.currentdesign;
          'erase'          'string'  { 'on','off' }   'off';
          'splitreg'       'string'  { 'on','off' }   'off';
          'interaction'    'string'  { 'on','off' }   'off';
          'freqlim'        'real'    []               [] ;
          'timelim'        'real'    []               [] ;
          'neighboropt'    'cell'    {}               {} ;
          'chanloc'        'struct'  {}               struct('no', {}); % default empty structure
          'neighbormat'    'real'    []               [] },...
          'std_limo');
    if isstr(opt), error(opt); end
end

Analysis     = opt.measure;
design_index = opt.design;

% Make sure paths are ok for LIMO (Consider to move this to eeglab.m in a future)
% -------------------------------------------------------------------------
local_path = which('limo_eeg');
root = fileparts(local_path);
addpath([root filesep 'limo_cluster_functions']);
addpath([root filesep 'external' filesep 'psom']);
addpath([root filesep 'external']);
addpath([root filesep 'help']);

% Checking fieldtrip paths
if ~exist('ft_prepare_neighbours')
    error('std_limo error: Fieldtrip extension must be installed');
end

% Detecting type of analysis
% -------------------------------------------------------------------------
if strncmp(Analysis,'dat',3)
    model.defaults.type = 'Channels';
elseif strncmp(Analysis,'ica',3)
    [STUDY,flags]=std_checkdatasession(STUDY,ALLEEG);
    if sum(flags)>0
        error('some subjects have data from different sessions - can''t do ICA');
    end
    model.defaults.type = 'Components';
    model.defaults.icaclustering = 1;
end

% Checking if clusters
% -------------------------------------------------------------------------
if strcmp(model.defaults.type,'Components') && isempty(STUDY.cluster(1).child)
    fprintf(2,'std_limo error: Unable to compute LIMO on unclustered component data \n');
    return;
end
    
% computing channel neighbox matrix
% ---------------------------------
flag_ok = 1;
if isempty(opt.chanloc) && isempty(opt.neighbormat)
    if isfield(ALLEEG(1).chanlocs, 'theta') &&  ~strcmp(model.defaults.type,'Components')
        if  ~isfield(STUDY.etc,'statistic')
            STUDY = pop_statparams(STUDY, 'default');
        end
        
        try
            [~,~,limoChanlocs] = std_prepare_neighbors(STUDY, ALLEEG, 'force', 'on', opt.neighboropt{:});
            chanlocname = 'limo_gp_level_chanlocs.mat';
        catch neighbors_error
            errordlg2(neighbors_error.message,'limo_gp_level_chanlocs.mat not created')
        end
    else
        disp('Warning: cannot compute expected channel distance for correction for multiple comparisons');
        limoChanlocs = [];
        flag_ok = 0;
    end
else
    limoChanlocs.expected_chanlocs   = opt.chanloc;
    limoChanlocs.channeighbstructmat = opt.neighbormat;
    chanlocname = 'limo_chanlocs.mat';
end
if flag_ok
    limoChanlocsFile = fullfile(STUDY.filepath, chanlocname);
    save('-mat', limoChanlocsFile, '-struct', 'limoChanlocs');
    fprintf('Saving channel neighbors for correction for multiple comparisons in %s\n', limoChanlocsFile);
end

% 1st level analysis
% -------------------------------------------------------------------------
model.cat_files = [];
model.cont_files = [];
unique_subjects = STUDY.design(1).cases.value'; % all designs have the same cases
nb_subjects     = length(unique_subjects);

for s = 1:nb_subjects
    nb_sets(s) = numel(find(strcmp(unique_subjects{s},{STUDY.datasetinfo.subject})));
end

% find out if the channels are interpolated
% -----------------------------------------
interpolated = zeros(1,length(STUDY.datasetinfo));
for iDat = 1:length(STUDY.datasetinfo)
    fileName = fullfile(STUDY.datasetinfo(iDat).filepath, [ STUDY.datasetinfo(iDat).subject '.' opt.measure ]);
    tmpChans = load('-mat', fileName, 'labels');
    if length(tmpChans.labels) > ALLEEG(iDat).nbchan, interpolated(iDat) = 1; end        
end

% simply reshape to read columns
% -------------------------------------------------------------------------
for s = 1:nb_subjects
    order{s} = find(strcmp(unique_subjects{s},{STUDY.datasetinfo.subject}));
end

% Cleaning old files from the current design (Cleaning ALL)
% -------------------------------------------------------------------------
% NOTE: Clean up the .lock files to (to be implemented)
if strcmp(opt.erase,'on')
    [tmp,filename] = fileparts(STUDY.filename);
    std_limoerase(STUDY.filepath, filename, unique_subjects, num2str(STUDY.currentdesign));
    STUDY.limo = [];
end

% Check if the measures has been computed
% -------------------------------------------------------------------------
for nsubj = 1 : length(unique_subjects)
    inds     = find(strcmp(unique_subjects{nsubj},{STUDY.datasetinfo.subject}));
    
    % Checking for relative path
    study_fullpath = rel2fullpath(STUDY.filepath,STUDY.datasetinfo(inds(1)).filepath);
    %---
    subjpath = fullfile(study_fullpath, [unique_subjects{nsubj} '.' lower(Analysis)]);  % Check issue when relative path (remove comment)
    if exist(subjpath,'file') ~= 2
        error('std_limo: Measures must be computed first');
    end
end
clear study_fullpath pathtmp;
 
measureflags = struct('daterp','off',...
                     'datspec','off',...
                     'datersp','off',...
                     'datitc' ,'off',...
                     'icaerp' ,'off',...
                     'icaspec','off',...
                     'icaersp','off',...
                     'icaitc','off');
                 
measureflags.(lower(Analysis))= 'on';
STUDY.etc.measureflags = measureflags;

% generate temporary merged datasets needed by LIMO
% -------------------------------------------------
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
    
    filename = [ STUDY.datasetinfo(order{s}(1)).subject '_limo_file_tmp' num2str(design_index) '.set'];
    index = [STUDY.datasetinfo(order{s}).index];
    tmp   = {STUDY.datasetinfo(order{s}).subject};
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
    
    file_fullpath = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).filepath);
    model.set_files{s} = fullfile(file_fullpath , filename);
    
    OUTEEG = [];    
    if all([ALLEEG(index).trials] == 1)
         OUTEEG.trials = 1;
    else OUTEEG.trials = sum([ALLEEG(index).trials]);
    end
    
    filepath_tmp = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).filepath);
    OUTEEG.filepath    = filepath_tmp;
    OUTEEG.filename    = filename;
    OUTEEG.srate       = ALLEEG(index(1)).srate;
    OUTEEG.icaweights  = ALLEEG(index(1)).icaweights;
    OUTEEG.icasphere   = ALLEEG(index(1)).icasphere;
    OUTEEG.icachansind = ALLEEG(index(1)).icachansind;
    OUTEEG.etc         = ALLEEG(index(1)).etc;
    OUTEEG.times       = ALLEEG(index(1)).times;
    if any(interpolated)
        OUTEEG.chanlocs    = mergedChanlocs;
        OUTEEG.etc.interpolatedchannels = setdiff([1:length(OUTEEG.chanlocs)], std_chaninds(OUTEEG, { ALLEEG(index(1)).chanlocs.labels }));
    else
        OUTEEG.chanlocs    = ALLEEG(index(1)).chanlocs;
    end
    
    % update EEG.etc
    OUTEEG.etc.merged{1} = ALLEEG(index(1)).filename;
    
    % Def fields
    OUTEEG.etc.datafiles.daterp   = [];
    OUTEEG.etc.datafiles.datspec  = [];
    OUTEEG.etc.datafiles.dattimef = [];
    OUTEEG.etc.datafiles.datitc   = [];
    OUTEEG.etc.datafiles.icaerp   = [];
    OUTEEG.etc.datafiles.icaspec  = [];
    OUTEEG.etc.datafiles.icatimef = [];
    OUTEEG.etc.datafiles.icaitc   = [];
    
    % Filling fields
    if isfield(ALLEEG(index(1)).etc.datafiles,'daterp')
        OUTEEG.etc.datafiles.daterp{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.daterp);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'datspec')
        OUTEEG.etc.datafiles.datspec{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.datspec);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'dattimef')
        OUTEEG.etc.datafiles.datersp{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.dattimef);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'datitc')
        OUTEEG.etc.datafiles.datitc{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.datitc);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'icaerp')
        OUTEEG.etc.datafiles.icaerp{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.icaerp);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'icaspec')
        OUTEEG.etc.datafiles.icaspec{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.icaspec);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'icatimef')
        OUTEEG.etc.datafiles.icaersp{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.icatimef);
    end
    if isfield(ALLEEG(index(1)).etc.datafiles,'icaitc')
        OUTEEG.etc.datafiles.icaitc{1} = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).etc.datafiles.icaitc);
    end
    
    % Save info
    EEG = OUTEEG;
    save('-mat', fullfile( filepath_tmp, OUTEEG.filename), 'EEG');
    clear OUTEEG filepath_tmp
end

% generate data files
% -------------------
factors = pop_listfactors(STUDY.design, 'gui', 'off');
for s = 1:nb_subjects     
    % save continuous and categorical data files
    % ************* PLURAL IMPORTANT IN FILE? continuous_variable.txt vs continuous_variables.txt
    %filepath_tmp = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).filepath);
    trialinfo = std_combtrialinfo(STUDY.datasetinfo, unique_subjects{s});
    [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, 'splitreg', opt.splitreg, 'interaction', opt.interaction); %, 'filepath', filepath_tmp); 

    % copy results
    model.cat_files{s}  = catMat;
    model.cont_files{s} = contMat;
    STUDY.limo.categorical             = limodesign.categorical;
    STUDY.limo.continuous              = limodesign.continuous;
    STUDY.limo.subjects(s).subject     = unique_subjects{s};
    STUDY.limo.subjects(s).cat_file    = catMat;
    STUDY.limo.subjects(s).cont_file   = contMat;
end
    
% transpose
model.set_files = model.set_files';
model.cat_files = model.cat_files';
model.cont_files = model.cont_files';
if all(cellfun(@isempty, model.cat_files )), model.cat_files  = []; end
if all(cellfun(@isempty, model.cont_files)), model.cont_files = []; end

% set model.defaults - all conditions no bootstrap
% -----------------------------------------------------------------
% to update passing the timing/frequency from STUDY - when computing measures
% -----------------------------------------------------------------
if strcmp(Analysis,'daterp') || strcmp(Analysis,'icaerp')
    model.defaults.analysis= 'Time';
    model.defaults.start = ALLEEG(index(1)).xmin*1000;
    model.defaults.end   = ALLEEG(index(1)).xmax*1000;
    if length(opt.timelim) == 2 && opt.timelim(1) < opt.timelim(end)
        % start value
        if opt.timelim(1) > model.defaults.start && opt.timelim(1) < model.defaults.end
            model.defaults.start = opt.timelim(1);
        else
            display('std_limo: Invalid time lower limit, using default value instead');
        end
        % end value
        if opt.timelim(end) < model.defaults.end && opt.timelim(end) > model.defaults.start
            model.defaults.end = opt.timelim(end);
        else
            display('std_limo: Invalid time upper limit, using default value instead');
        end
    end
    model.defaults.lowf  = [];
    model.defaults.highf = [];
    
elseif strcmp(Analysis,'datspec') || strcmp(Analysis,'icaspec')
    
    model.defaults.analysis= 'Frequency';
    model.defaults.start   = -10;
    model.defaults.end     = ALLEEG(index(1)).xmax*1000;
    if length(opt.freqlim) == 2 && opt.freqlim(1) < opt.freqlim(end)
        model.defaults.lowf    = opt.freqlim(1);
        model.defaults.highf   = opt.freqlim(end);
    else
        model.defaults.lowf    = [];
        model.defaults.highf   = [];
    end;
    
elseif strcmp(Analysis,'datersp') || strcmp(Analysis,'icaersp')
    model.defaults.analysis = 'Time-Frequency';
    model.defaults.start    = [];
    model.defaults.end      = [];
    model.defaults.lowf     = [];
    model.defaults.highf    = [];
end

model.defaults.fullfactorial = 0;                    % factorial only for single subject analyses - not included for studies
model.defaults.zscore = 0;                           % done that already
model.defaults.bootstrap = 0 ;                       % only for single subject analyses - not included for studies
model.defaults.tfce = 0;                             % only for single subject analyses - not included for studies
model.defaults.method = opt.method;                  % default is OLS - to be updated to 'WLS' once validated
model.defaults.Level= 1;                             % 1st level analysis
model.defaults.type_of_analysis = 'Mass-univariate'; % future version will allow other techniques

limocontrast.mat = []; %{ [1 1 -1 -1] [ 1 -1] };
[LIMO_files, procstatus] = limo_batch('model specification',model,limocontrast,STUDY);
STUDY.limo.model    = model;
STUDY.limo.datatype = Analysis;
STUDY.limo.chanloc  = limoChanlocs.expected_chanlocs;

%LIMO_files.expected_chanlocs = limoChanlocs;
%procOK = find(procstatus);


% 
% if isempty(find(procstatus==0)) % all test succeded
%     % Computing univariate one sample t-test for each parameters
%     % ----------------------------------------------------------
%     nparams = 0;
%     if exist('categ','var'),   nparams = max(categ); end;
%     if exist('contvar','var'), nparams = nparams + size(contvar,2); end;
%     
%     foldername = [STUDY.filename(1:end-6) '_GLM' num2str(STUDY.currentdesign) '_' model.defaults.type '_' model.defaults.analysis '_'];
%     nbootval   = 0;
%     tfceval    = 0;
%     
%     try
%         filesout = limo_random_select(1,LIMO_files.expected_chanlocs,'nboot'         ,nbootval...
%                                                                     ,'type'          ,model.defaults.type ...
%                                                                     ,'tfce'          ,tfceval...
%                                                                     ,'analysis_type' ,'fullchan'...
%                                                                     ,'parameters'    ,{[1:nparams]}...
%                                                                     ,'limofiles'     ,{LIMO_files.Beta}...
%                                                                     ,'folderprefix'  ,foldername...
%                                                                     ,'folderpath'    ,LIMO_files.LIMO);
%         testype    = cell([length(filesout),1]);
%         testype(:) = {'ttest'};
%         for i = 1: nparams
%             testname{i,1}   = [testype{i} '_parameter' num2str(i)];
%         end
% 
%         % 2- Computing Paired ttest for each combination of parameters
%         ncomb        = combnk(1:nparams,2);
%         limofiles{1} = LIMO_files.Beta;
%         limofiles{2} = LIMO_files.Beta;
% 
%         for i = 1:size(ncomb,1)
%             tmpname          = [foldername 'par_' num2str(ncomb(i,1)) '_' num2str(ncomb(i,2))];
%             mkdir(LIMO_files.LIMO,tmpname);
%             folderpath       = fullfile(LIMO_files.LIMO,tmpname);
%             filesout{end+1}  = limo_random_select(2,LIMO_files.expected_chanlocs,'nboot'         ,nbootval...
%                                                                                 ,'type'          ,model.defaults.type ...
%                                                                                 ,'tfce'          ,tfceval...
%                                                                                 ,'analysis_type' ,'fullchan'...
%                                                                                 ,'parameters'    ,{[ncomb(i,1)] [ncomb(i,2)]}...
%                                                                                 ,'limofiles'     ,limofiles...
%                                                                                 ,'folderpath'    ,folderpath);
%             testype{end+1}      = 'p_ttest';
%             testname{end+1,1}   = [testype{end} '_parameter_' num2str(ncomb(i,1)) '-' num2str(ncomb(i,2))];
%         end
%         % Assigning Level 2 info to STUDY
%         STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).l2files = [testype testname filesout' ];
% 
%         for i = 1:size(filesout,2)
%             STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).groupmodel(i).filename = filesout{1,i};
%             STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).groupmodel(i).guiname  = testname{i,1};
%         end
%     catch,
%         disp('2nd level LIMO function failed - call from the command line');
%     end
%  end
% 
% % Saving STUDY
% STUDY = pop_savestudy( STUDY, [],'filepath', STUDY.filepath,'savemode','resave');
% cd(STUDY.filepath);
% end

% -------------------------------------------------------------------------
% Return full path if 'filepath' is a relative path. The output format will
% fit the one of 'filepath'. That means that if 'filepath' is a cell array,
% then the output will a cell array too, and the same if is a string.
function file_fullpath = rel2fullpath(studypath,filepath)

nit = 1; if iscell(filepath), nit = length(filepath);end

for i = 1:nit
    if iscell(filepath),pathtmp = filepath{i}; else pathtmp = filepath; end
    if strfind(pathtmp(end),filesep), pathtmp = pathtmp(1:end-1); end % Getting rid of filesep at the end
    if ~isempty(strfind(pathtmp(1:2),['.' filesep])) || (isunix && pathtmp(1) ~= '/') || (ispc && pathtmp(2) ~= ':')
        if iscell(filepath),
            file_fullpath{i} = fullfile(studypath,pathtmp(1:end));
        else
            file_fullpath = fullfile(studypath,pathtmp(1:end));
        end
    else
        if iscell(filepath),
            file_fullpath{i} = pathtmp;
        else
            file_fullpath = pathtmp;
        end
    end
end
