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
%   'freqlim'     - Frequency trimming
%   'timelim'     - Time trimming
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
% Author: Cyril Pernet (LIMO Team), The university of Edinburgh, 2014
%         Ramon Martinez-Cancino and Arnaud Delorme, SCCN, 2014
%
% Copyright (C) 2015  Ramon Martinez-Cancino,INC, SCCN
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

cd(STUDY.filepath);
if nargin < 2
    help std_limo;
    return;
end;

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
          'freqlim'        'real'    []               [] ;
          'timelim'        'real'    []               [] ;
          'neighboropt'    'cell'    {}               {} ;
          'chanloc'        'struct'  {}               struct('no', {}); % default empty structure
          'neighbormat'    'real'    []               [] },...
          'std_limo');
    if isstr(opt), error(opt); end;
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
        [tmp1 tmp2 limostruct] = std_prepare_neighbors(STUDY, ALLEEG, 'force', 'on', opt.neighboropt{:});
        chanlocname = 'limo_gp_level_chanlocs.mat';
    else
        disp('Warning: cannot compute expected channel distance for correction for multiple comparisons');
        limoChanlocs = [];
        flag_ok = 0;
    end
else
    limostruct.expected_chanlocs   = opt.chanloc;
    limostruct.channeighbstructmat = opt.neighbormat;
    chanlocname = 'limo_chanlocs.mat';
end

if flag_ok
    limoChanlocs = fullfile(STUDY.filepath, chanlocname);
    save('-mat', limoChanlocs, '-struct', 'limostruct');
    fprintf('Saving channel neighbors for correction for multiple comparisons in %s\n', limoChanlocs);
end
clear flag_ok tmp1 tmp2

% 1st level analysis
% -------------------------------------------------------------------------
model.cat_files = [];
model.cont_files = [];
if isempty(STUDY.design(design_index).filepath)
    STUDY.design(design_index).filepath = STUDY.filepath;
end
unique_subjects = STUDY.design(design_index).cases.value';
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
    if exist([STUDY.filepath filesep 'limo_batch_report'],'dir')
        try
            rmdir([STUDY.filepath filesep 'limo_batch_report'],'s');
        catch
        end
    end
    
    [tmp,filename] = fileparts(STUDY.filename);
    
    % Cleaning level 1 folders
    for i = 1:nb_subjects
        tmpfiles = dir([STUDY.filepath filesep 'LIMO_' filename filesep unique_subjects{i} filesep 'GLM' num2str(STUDY.currentdesign) '*']);
        tmpfiles = {tmpfiles.name};
        if ~isempty(tmpfiles)
            for j = 1:length(tmpfiles)
                try
                rmdir([STUDY.filepath filesep 'LIMO_' filename filesep unique_subjects{i} filesep tmpfiles{j}],'s');
                catch
                end
            end
        end
    end
    
    % Cleaning level 2 folders
    if isfield(STUDY.design(design_index),'limo') && isfield(STUDY.design(design_index).limo,'groupmodel')
        for nlimos = 1:length(STUDY.design(design_index).limo)
            for ngroupmodel = 1:size(STUDY.design(design_index).limo(nlimos).groupmodel,2)
                path2file = fileparts(STUDY.design(design_index).limo(nlimos).groupmodel(ngroupmodel).filename);
                try
                    rmdir(path2file,'s');
                catch
                    fprintf(2,['Fail to remove. Folder ' path2file ' does not exist or has been removed\n']);
                end
            end
            
        end
    end
    
    % Cleaning files related with design (innecesary??)
    tmpfiles = dir([STUDY.filepath filesep 'LIMO_' filename]);
    tmpfiles = {tmpfiles.name};
    
    filecell = strfind(tmpfiles,['GLM' num2str(STUDY.currentdesign)]);
    indxtmp = find(not(cellfun('isempty',filecell)));
    if ~isempty(indxtmp)
        for nfile = 1:length(indxtmp)
            try
                file2delete = [STUDY.filepath filesep 'LIMO_' filename filesep tmpfiles{indxtmp(nfile)}];
                if isdir(file2delete)
                    rmdir(file2delete,'s');
                else
                    delete(file2delete);
                end
            catch
                fprintf(2,['Fail to remove. File ' tmpfiles{indxtmp(nfile)} ' does not exist or has been removed\n']);
            end
        end
    end
    
    % Cleaning limo fields
    STUDY.design(STUDY.currentdesign).limo = [];
    
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

% Checking if continuous variables
% -------------------------------------------------------------------------
if isfield(STUDY.design(design_index).variable,'vartype')
    if any(strcmp(unique({STUDY.design(design_index).variable.vartype}),'continuous'))
        cont_var_flag = 1;
    else
        cont_var_flag = 0;
    end
else
    error('std_limo: Define a valid design');
end
% -------------------------------------------------------------------------
mergedChanlocs = eeg_mergelocs(ALLEEG.chanlocs);
for s = 1:nb_subjects     
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
    
    % field which are needed
    % EEGLIMO.etc
    % EEGLIMO.times
    % EEGLIMO.chanlocs
    % EEGLIMO.srate
    % EEGLIMO.filepath
    % EEGLIMO.filename
    % EEGLIMO.icawinv
    % EEGLIMO.icaweights
    
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
    
    catvar_matrix = std_lm_getvars(STUDY,STUDY.datasetinfo(order{s}(1)).subject,'design',design_index,'vartype','cat');

    if ~isempty(catvar_matrix)
        if isvector(catvar_matrix) % only one factor of n conditions
            categ = catvar_matrix;
            model.cat_files{s} = catvar_matrix;
%         elseif 1
%             
%             trialinfo = std_combtrialinfo(STUDY.datasetinfo, order{s});
%             categ = std_builddesignmat(STUDY.design(design_index), trialinfo);
%             model.cat_files{s} = catvar_matrix;
%             
        else % multiple factors (=multiple columns)
            X = []; nb_row = size(catvar_matrix,1);
            for cond = 1:size(catvar_matrix,2)
                nb_col(cond) = max(unique(catvar_matrix(:,cond)));
                x = NaN(nb_row,nb_col(cond));
                for c=1:nb_col(cond)
                    x(:,c) = [catvar_matrix(:,cond) == c];
                end
                X = [X x];
            end
            
            tmpX = limo_make_interactions(X, nb_col); dim_interactions(s) = size(tmpX,2);
            tmpX = tmpX(:,sum(nb_col)+1:end);
            getnan = find(isnan(catvar_matrix(:,1)));
            if ~isempty(getnan)
                tmpX(getnan,:) = NaN;
            end
            categ = sum(repmat([1:size(tmpX,2)],nb_row,1).*tmpX,2);
            clear x X tmpX nb_col; 
            model.cat_files{s} = categ;
        end 
        filepath_tmp = rel2fullpath(STUDY.filepath,ALLEEG(index(1)).filepath);
        save(fullfile(filepath_tmp, 'categorical_variable.txt'),'categ','-ascii');
    end
    
    % model.cont_files: a cell array of continuous variable files
    if cont_var_flag
        
        [contvar,contvar_info] = std_lm_getvars(STUDY,STUDY.datasetinfo(order{s}(1)).subject,'design',design_index,'vartype','cont'); clear tmp; %#ok<ASGLU>
        
        if nb_sets ~= 1;
%             [contvar_matrix,contvar_info] = std_lm_getvars(STUDY,STUDY.datasetinfo(order(1,s)).subject,'design',design_index,'vartype','cont');
            contvar_matrix = contvar; clear contvar;
            % we need to know the nb of trials in each of the sets
            set_positions = [0];
            getnan = find(isnan(contvar_matrix(:,1)));
            for dataset = 1:nb_sets
                cond_size = cellfun(@length,contvar_info.datasetinfo_trialindx(contvar_info.dataset == index(dataset)));
                set_size = sum(cond_size);
                set_positions(dataset+1) = set_size + sum(getnan<set_size);
                % update getnan as the next loop restart at one
                getnan((getnan - set_size)<0) = [];
            end
            set_positions = cumsum(set_positions);
            
            % create contvar updating with set_positions
            contvar = NaN(size(contvar_matrix));
            for cond = 1:size(contvar_info.datasetinfo_trialindx,2)
                which_dataset = contvar_info.dataset(cond);                    % dataset nb
                dataset_nb = find(which_dataset == order{s});                  % position in the concatenated data
                position = cell2mat(contvar_info.datasetinfo_trialindx(cond)); % position in the set
                values = contvar_matrix(position+set_positions(dataset_nb),:); % covariable values
                contvar(position+set_positions(dataset_nb),:) = values;        % position in the set + position in the concatenated data
            end
        end

        % --> split per condition and zscore (if requested)
        if exist('categ','var')
            if strcmp(opt.splitreg,'on')
                model.cont_files{s} = limo_split_continuous(categ,contvar);
            else
                model.cont_files{s} = contvar;
            end
        end
        
        if size(contvar,2) > 1
            save([ALLEEG(index(1)).filepath filesep 'continuous_variables.txt'],'contvar','-ascii')
        else
            save([ALLEEG(index(1)).filepath filesep 'continuous_variable.txt'],'contvar','-ascii')
        end
    end
end

if exist('dim_interactions','var') == 1 && length(unique(dim_interactions)) ~= 1
   error('Number of interactions are not the same. Check design');
end

model.set_files = model.set_files';
model.cat_files = model.cat_files';
if cont_var_flag
    model.cont_files = model.cont_files';
end

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

STUDY.design_info = STUDY.design(design_index).variable;
STUDY.design_index = design_index; % Limo batch have to be fixed to use the field STUDY.currentdesig
STUDY.names = names;

% Contrast
% -------------------------------------------------------------------------
% if cont_var_flag && exist('categ','var')
%     if strcmp(opt.splitreg,'on')
%         limocontrast.mat = [zeros(1,max(categ)) ones(1,max(categ)) 0]; % this already computed ie this the 1.*Beta ! 
%     else
%         limocontrast.mat = [zeros(1,max(categ))  ones(1,size(contvar,2)) 0];  % doesn't make sense to add continuous regressors 
%     end
%     LIMO_files = limo_batch('both',model,limocontrast,STUDY);
% else
%     limocontrast.mat = [];
%     LIMO_files = limo_batch('model specification',model,limocontrast,STUDY);
% end
limocontrast.mat = [];
[LIMO_files, procstatus] = limo_batch('model specification',model,limocontrast,STUDY);
LIMO_files.expected_chanlocs = limoChanlocs;
procOK = find(procstatus);

for i = 1:length(procOK)
    limofolders{i} = fileparts(LIMO_files.mat{procOK(i)});
end

% Getting indices to save LIMO in STUDY structure (save multiples analysis)
% -------------------------------------------------------------------------
if isfield(STUDY.design(STUDY.currentdesign),'limo') && isfield(STUDY.design(STUDY.currentdesign).limo, 'datatype')
    stdlimo_indx = find(strcmp({STUDY.design(STUDY.currentdesign).limo.datatype},Analysis));
     if isempty(stdlimo_indx)
         stdlimo_indx = length(STUDY.design(STUDY.currentdesign).limo) + 1;
     end   
else
    stdlimo_indx = 1;
end

% Cleaning
% -------------------------------------------------------------------------
rmfield(STUDY,'design_index');
rmfield(STUDY,'design_info');
rmfield(STUDY,'names');

% Clenainig level 2 folders associated to specific analaysis if opt.erase =
% 'off' (IF EXIST BEFORE)
if strcmp(opt.erase,'off')
    if isfield(STUDY.design(STUDY.currentdesign),'limo') && length(STUDY.design(STUDY.currentdesign).limo) >= stdlimo_indx
        STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).model = [];
        STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).groupmodel = [];
    end
end

% Assigning info to STUDY
% -------------------------------------------------------------------------
STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).datatype = Analysis;
STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).chanloc  = LIMO_files.expected_chanlocs;

for i =1:length(procOK)
    STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).basefolder_level1{procOK(i)}  = fileparts(LIMO_files.mat{procOK(i)});
    tmp = dir(fileparts(LIMO_files.mat{procOK(i)}));
    if ~isempty({tmp.name})
        outputfiles = {tmp.name}';
        
        % Cleaning out the list
        %----------------------
        cleanoutlist = {'Betas.mat','LIMO.mat','Yr.mat','Yhat.mat','Res.mat','.','..'};
        for j = 1:length(cleanoutlist)
            ind2delete{j} = find(strcmp(outputfiles,cleanoutlist{j}));
        end
        nonemptyvals              = find(cellfun(@(x) ~isempty(x), ind2delete));
        ind2delete                = [ind2delete{nonemptyvals}];
        outputfiles(ind2delete) = []; 
        clear ind2delete;
    end
    for nfiles = 1:size(outputfiles,1)
        STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).model(procOK(i)).filename{nfiles,1}  = fullfile(fileparts(LIMO_files.mat{procOK(i)}),outputfiles{nfiles}); 
        [trash,tmp] = fileparts(outputfiles{nfiles});
        STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).model(procOK(i)).guiname{nfiles,1} = tmp;
    end
end

if ~isempty(STUDY.filepath)
    cd(STUDY.filepath);
end;
if isempty(find(procstatus==0))
    % Computing one sample t-test for each parameters
    % -------------------------------------------------------------------------
    % NOTE:
    % if a Group lvel variable exist, we need to split data accordingly (TO DO)
    
    % 1 - Computing ttest for each parameter
    nparams = 0;
    if exist('categ','var'),   nparams = max(categ); end;
    if exist('contvar','var'), nparams = nparams + size(contvar,2); end;
    
    foldername = [STUDY.filename(1:end-6) '_GLM' num2str(STUDY.currentdesign) '_' model.defaults.type '_' model.defaults.analysis '_'];
    nbootval   = 0;
    tfceval    = 0;
    
    try
        filesout = limo_random_select(1,LIMO_files.expected_chanlocs,'nboot'         ,nbootval...
                                                                    ,'type'          ,model.defaults.type ...
                                                                    ,'tfce'          ,tfceval...
                                                                    ,'analysis_type' ,'fullchan'...
                                                                    ,'parameters'    ,{[1:nparams]}...
                                                                    ,'limofiles'     ,{LIMO_files.Beta}...;
                                                                    ,'folderprefix'  ,foldername...
                                                                    ,'folderpath'    ,LIMO_files.LIMO);
        testype    = cell([length(filesout),1]);
        testype(:) = {'ttest'};
        for i = 1: nparams
            testname{i,1}   = [testype{i} '_parameter' num2str(i)];
        end

        % 2- Computing Paired ttest for each combination of parameters
        ncomb        = combnk(1:nparams,2);
        limofiles{1} = LIMO_files.Beta;
        limofiles{2} = LIMO_files.Beta;

        for i = 1:size(ncomb,1)
            tmpname          = [foldername 'par_' num2str(ncomb(i,1)) '_' num2str(ncomb(i,2))];
            mkdir(LIMO_files.LIMO,tmpname);
            folderpath       = fullfile(LIMO_files.LIMO,tmpname);
            filesout{end+1}  = limo_random_select(2,LIMO_files.expected_chanlocs,'nboot'         ,nbootval...
                                                                                ,'type'          ,model.defaults.type ...
                                                                                ,'tfce'          ,tfceval...
                                                                                ,'analysis_type' ,'fullchan'...
                                                                                ,'parameters'    ,{[ncomb(i,1)] [ncomb(i,2)]}...
                                                                                ,'limofiles'     ,limofiles...
                                                                                ,'folderpath'    ,folderpath);
            testype{end+1}      = 'p_ttest';
            testname{end+1,1}   = [testype{end} '_parameter_' num2str(ncomb(i,1)) '-' num2str(ncomb(i,2))];
        end
        % Assigning Level 2 info to STUDY
        STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).l2files = [testype testname filesout' ];

        for i = 1:size(filesout,2)
            STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).groupmodel(i).filename = filesout{1,i};
            STUDY.design(STUDY.currentdesign).limo(stdlimo_indx).groupmodel(i).guiname  = testname{i,1};
        end
    catch,
        disp('2nd level LIMO function failed - call from the command line');
    end
 end

% Saving STUDY
STUDY = pop_savestudy( STUDY, [],'filepath', STUDY.filepath,'savemode','resave');
cd(STUDY.filepath);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function file_fullpath = rel2fullpath(studypath,filepath)
% Return full path if 'filepath' is a relative path. The output format will
% fit the one of 'filepath'. That means that if 'filepath' is a cell array,
% then the output will a cell array too, and the same if is a string.

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
end