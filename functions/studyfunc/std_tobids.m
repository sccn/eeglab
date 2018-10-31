function STUDY = std_tobids(varargin)

% From a STUDY rename files if needed, sort folders, and export metadata
% following the BIDS specification
%
% FORMAT STUDY = std_tobids(STUDY,options)
%        options are: 'export_dir' followed by the path
%                     'keep_files' being 'yes' (default) or 'no'
%                     'update_study' being 'yes' or 'no'
%
% Author: Cyril Pernet - LIMO Team, University of Edinurgh
%         Arnaud Delorme, EEGLAB, SCCN, 2018
%

%% options
if nargin == 1
    export_dir   = STUDY.filepath; % where to make the bids directory
    keep_files   = 'no';          % copy files to move and rename
    update_study = 'yes';          % return study with updated names and paths
else
    for in = 2:nargin
        if strcmpi(varargin{in},'export_dir') || strcmpi(varargin{in},'exportdir')
            export_dir = varargin{in+1};
        elseif strcmpi(varargin{in},'keep_files') || strcmpi(varargin{in},'keepfiles')
            keep_files =  varargin{in+1};
        elseif strcmpi(varargin{in},'update_study') || strcmpi(varargin{in},'updatestudy')
            update_study = varargin{in+1};
        end
    end
end

%% get subjects info
N = length(STUDY.subject);
if isempty(STUDY.task)
    STUDY.task = inputdlg2('input task name','task info missing'); % @Arno inputdlg2 fails
end
% task name with no space
task = lower(STUDY.task);
task(1) = upper(task(1));
task(find(isspace(STUDY.task))+1) = upper(task(find(isspace(STUDY.task))+1));
task(isspace(task)) = [];


% let's get the index of which files belong to the same subject
datasetindex = cell(N,length(STUDY.session));
if length(STUDY.session) > 1
    for subject =1:N
        file_index = 1;
        for files = 1:length(STUDY.datasetinfo)
            if strcmp(STUDY.subject{subject},STUDY.datasetinfo(files).subject)
                datasetindex{subject,file_index} = files;
                file_index = file_index+1;
            end
        end
    end
else
    for subject =1:N
        datasetindex{subject} = STUDY.subject{subject};
    end
end


% % check if subjects are all together or in different folders
% for files = 1:length(STUDY.datasetinfo)
%     paths{files} = STUDY.datasetinfo(files).filepath;
% end

% for each subject, move/rename and export metadata
for subject =1:N
    if strncmp(STUDY.subject{subject},'sub-',4)
        newname = STUDY.subject{subject};
    else
        newname = ['sub-' STUDY.subject{subject}];
    end
    
    % create sub- folder
    if strcmp(export_dir(length(export_dir)),filesep)
        export_dir = export_dir(1:end-1);
    end
    subj_dir             = [export_dir filesep newname filesep 'eeg'];
    subj_derivatives_dir = [export_dir filesep 'derivatives' filesep newname filesep 'eeg'];
    mkdir(subj_dir)
    
    % rename/move files
    for run = 1:length(STUDY.session)
        if ~isempty(datasetindex{subject,run})
            if ischar(datasetindex{subject,run})
                tmp = STUDY.datasetinfo(str2num(datasetindex{subject,run}));
            else
                tmp = STUDY.datasetinfo(datasetindex{subject,run});
            end
            [~,name,ext]=fileparts(tmp.filename);
            
            if strcmpi(keep_files,'no')
                movefile([tmp.filepath filesep tmp.filename],[subj_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_eeg' ext]);
                movefile([tmp.filepath filesep name '.fdt'], [subj_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_eeg.fdt']);
            elseif strcmpi(keep_files,'yes')
                copyfile([tmp.filepath filesep tmp.filename],[subj_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_eeg' ext]);
                copyfile([tmp.filepath filesep name '.fdt'], [subj_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_eeg.fdt']);
            end
            EEG              = pop_loadset([subj_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_eeg' ext]);
            EEG.setname      = [subj_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_eeg'];
            EEG.datfile      = [newname '_task-' task '_sess-1_run-' num2str(run) '_eeg.fdt'];
            % metadata
            SF(subject,run)  = EEG.srate;
            CC(subject,run)  = size(EEG.chanlocs,2);
            Ref{subject,run} = EEG.ref;
            
            % check derivatives and update .set
            if exist([tmp.filepath filesep STUDY.subject{subject} '.daterp'],'file') || ...
                    exist([tmp.filepath filesep STUDY.subject{subject} '.daterpim'],'file') || ...
                    exist([tmp.filepath filesep STUDY.subject{subject} '.datspec'],'file') || ...
                    exist([tmp.filepath filesep STUDY.subject{subject} '.dattimef'],'file') 
                name = STUDY.subject{subject}; % otherwise comes from STUDY.datasetinfo, depends if computed from study or not
            end
            
            if exist([tmp.filepath filesep name '.daterp'],'file')
                if ~exist(subj_derivatives_dir,'dir'); mkdir(subj_derivatives_dir); end
                if strcmpi(keep_files,'no')
                    movefile([tmp.filepath filesep name '.daterp'], [subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-erp_eeg.daterp']);                    
                elseif strcmpi(keep_files,'yes')
                    copyfile([tmp.filepath filesep name '.daterp'], [subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-erp_eeg.daterp']);
                end
                EEG.etc.datafiles.daterp = [subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-erp_eeg.daterp'];
            end
            
            if exist([tmp.filepath filesep name '.daterpim'],'file')
                if ~exist(subj_derivatives_dir,'dir')' mkdir(subj_derivatives_dir); end
                if strcmpi(keep_files,'no')
                    movefile([tmp.filepath filesep name '.daterpim'],[subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-erpimg_eeg.daterpim']);
                elseif strcmpi(keep_files,'yes')
                    copyfile([tmp.filepath filesep name '.daterpim'],[subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-erpimg_eeg.daterpim']);
                end
                EEG.etc.datafiles.daterpim = [subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-erpimg_eeg.daterpim'];
            end
            
            if exist([tmp.filepath filesep name '.datspec'],'file')
                if ~exist(subj_derivatives_dir,'dir')' mkdir(subj_derivatives_dir); end
                if strcmpi(keep_files,'no')
                    movefile([tmp.filepath filesep name '.datspec'],[subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-spectrum_eeg.datspec']);
                elseif strcmpi(keep_files,'yes')
                    copyfile([tmp.filepath filesep name '.datspec'],[subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-spectrum_eeg.datspec']);
                EEG.etc.datafiles.datspec = [subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-spectrum_eeg.datspec'];
                end
            end
            
            if exist([tmp.filepath filesep name '.dattimef'],'file')
                if ~exist(subj_derivatives_dir,'dir')' mkdir(subj_derivatives_dir); end
                if strcmpi(keep_files,'no')
                    movefile([tmp.filepath filesep name '.dattimef'],[subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-timefrequency_eeg.dattimef']);
                elseif strcmpi(keep_files,'yes')
                    copyfile([tmp.filepath filesep name '.dattimef'],[subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-timefrequency_eeg.dattimef']);
                EEG.etc.datafiles.dattimef = [subj_derivatives_dir filesep newname '_task-' task '_sess-1_run-' num2str(run) '_desc-timefrequency_eeg.dattimef'];
                end
            end
            
            % resave the .set with updated info
            pop_saveset(EEG,'savemode','resave');
            
            % export metadata 
            if run == 1
                electrodes_to_tsv(EEG);
            end
           channelloc_to_tsv(EEG);
           events_to_tsv(EEG);
        end
        
        if strcmpi(update_study,'yes')
            STUDY.datasetinfo(subject).filepath = subj_dir;
            STUDY.datasetinfo(subject).filename = [newname '_task-' task '_sess-1_run-' num2str(run) '_eeg' ext];
            STUDY.datasetinfo(subject).subject  = newname;
            STUDY.subject{subject}              = newname;
        end
    end
end

% update task name
if strcmpi(update_study,'yes')
    STUDY.task = task;
    STUDY = pop_savestudy(STUDY,'savemode','resave');
end
        
% finish of with root metadata
% ----------------------------
% make *_eeg.json
PLF = str2num(cell2mat(inputdlg('What was the Power Line Frequency?','BIDS spec requirement')));
if size(unique(Ref),1) > 1
    warning('difference referencing scheme found of this dataset')
    Ref = cell2mat(inputdlg('What was the reference?','BIDS spec requirement'));
end

json = struct('TaskName',task, ...
    'SamplingFrequency', mean(SF(:)), ...
    'EEGChannelCount', max(CC(:)), ...
    'EEGReference', unique(Ref), ...
    'PowerLineFrequency', PLF, ...
    'SoftwareFilters', ' ');
jsonwrite([export_dir filesep task '_eeg.json'],json,struct('indent','  '))

% make a participants table and save 
age = zeros(N,1);
sex = repmat(' ',[N 1]); 
t = table(STUDY.subject',age,sex,'VariableNames',{'participant_id','age','sex'});
writetable(t,[export_dir filesep 'Participants.tsv'],'FileType','text','Delimiter','\t');
warndlg('Job done, but metadata created need editing','BIDS spec','modal')




