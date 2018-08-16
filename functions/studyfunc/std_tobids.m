% std2bids(study)
%
% from a STUDY rename files if needed, sort folders, and export metadata
% following the BIDS specification
%
% 
% Author: Cyril Pernet - LIMO Team, University of Edinurgh
%         Arnaud Delorme, EEGLAB, SCCN, 2018

function std_tobids.m(study)

%% get subjects info
N = length(STUDY.subject); 
if isempty(STUDY.task)
    STUDY.task = inputdlg('input task name','task info missing'); % @Arno inputdlg2 fails
end
task = STUDY.task(~isspace(STUDY.task));

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


% check if subjects are all together or in different folders
for files = 1:length(STUDY.datasetinfo)
      paths{files} = STUDY.datasetinfo(files).filepath;
end
        
% if altogether 
if length(unique(paths)) == 1
    root = paths{1};
    for subject =1:N
        if strncmp(STUDY.subject{subject},'sub-',4)
            newname = STUDY.subject{subject};
        else
            newname = ['sub-' STUDY.subject{subject}];
        end
        % create sub- folder
        mkdir(newname);
        % rename/move files
        for run = 1:length(STUDY.session) % by default assumes a single session - should update to distinguish this
            if ~isempty(datasetindex{subject,sess})
                tmp = STUDY.datasetinfo(datasetindex{subject,sess});
                [~,name,ext]=fileparts(tmp.filename);
                movefile([tmp.filepath filesep tmp.filename],[root filesep newname filesep newname '_task-' task '_run-' num2str(run) '_eeg' ext]);
                movefile([tmp.filepath filesep name '.fdt'],[root filesep newname filesep newname '_task-' task '_run-' num2str(run) '_eeg.fdt']);
                % check for derivatives
                if exist([STUDY.subject{subject} '.daterp'],'file')
                    if ~exist([root filesep 'derivatives'],'dir')
                        mkdir([root filesep 'derivatives'])
                    end
                    
                    if ~exist([root filesep 'derivatives' filesep newname],'dir')
                        mkdir([root filesep 'derivatives' filesep newname]);
                    end
                    
                    movefile([STUDY.subject{subject} '.daterp'],[root filesep 'derivatives' filesep newname filesep newname '_task-' task '_run-' num2str(run) '_eeg.daterp']);
                    
                    % add other cases here
                end
                
                % export metadata
                events_to_tsv
                channelloc_to_tvs
                electrodes_to_tsv
            end
        end
    end
end

STUDY.datasetinfo
STUDY.session

% rename 


%% export metadata



 