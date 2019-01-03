% bids_format('/Users/arno/temp/rishikesh', 24, {
% Optional inputs:
%  'targetdir' - [string] target directory. Default is 'bidsexport' in the
%                current folder.
%
%  'taskname'  - [string] name of the task. No space are allowed and no
%                special characters. Default is ''.
%
%  'README'    - [string] content of the README file. If the string points
%                to a file that exists on the path, the file is copied.
%                Otherwise a new README file is created and the content of
%                this variable is copied to the new file.
%
%  'CHANGES'   - [string] content of the README file. If the string points
%                to a file that exists on the path, the file is copied.
%
%  'codefiles' - [cell] cell array of file names containing code related
%                to importing or processing the data.
%
%  'stimuli'   - [cell] cell array of type and corresponding file names.   XXXXXXXXXXXX
%                For example: { 'sound1' '/Users/xxxx/sounds/sound1.mp3';  XXXXXXXXXXXX
%                               'img1'   '/Users/xxxx/sounds/img1.jpg' }   XXXXXXXXXXXX
%                (the semicolumn above is optional)
%
%  'gInfo'     - [struct] general information fields. See BIDS specifications.
%                For example.
%                info.ReferencesAndLinks = { 'Pubmed 122032020' };
%                info.Name = 'This is a custom task';
%                info.License = 'Creative commons';
%
%  'tInfo'     - [struct] task information fields. See BIDS specifications.
%                For example.
%                tInfo.InstitutionAddress = '9500 Gilman Drive, CA92093-0559 La Jolla, USA';
%                tInfo.InstitutionName = 'Univesity of California, San Diego';
%                tInfo.InstitutionalDepartmentName = 'Institute of Neural Computation';
%                tInfo.PowerLineFrequency = 50;
%
%  'pInfo'     - [cell] cell array of participant values, with one row
%                per participants. The first row contains columns names.
%                For example for 2 participants:
%                { 'sex'     'age';
%                  'male'    20;
%                  'female'  25 };
%
%  'pInfoDesc' - [struct] structure describing participant file columns.
%                The fields are described in the BIDS format.
%                For example
%                pInfo.participant_id.LongName = 'Participant ID';
%                pInfo.participant_id.Description = 'Event onset';
%                pInfo.sex.LongName = 'Gender';
%                pInfo.sex.Description = 'Gender';
%                pInfo.age.LongName = 'Age';
%                pInfo.age.Description = 'Age in years';
%
%  'eInfo'     - [cell] additional event information columns and their corresponding
%                event fields in the EEGLAB event structure. Note that
%                EEGLAB event latency, duration, and type are inserted
%                automatically as columns "onset" (latency in sec), "duration"
%                (duration in sec), "event_sample" (latency), "event_type"
%                (time locking event), "event_value" (type). For example
%                { 'HED' 'reaction_time';
%                  'HED' 'rt' }
%
%  'eInfoDesc' - [struct] structure describing additional or/and original
%                event fields if you wish to redefine these.
%                These are EEGLAB event fields listed above.
%                See the BIDS format for more information
%                eInfo.onset.LongName = 'Event onset'; % change default value
%                eInfo.onset.Description = 'Event onset';
%                eInfo.onset.Units = 'seconds';
%                eInfo.HED.LongName = 'Hierarchival Event Descriptors';
%                eInfo.reaction_time.LongName = 'Reaction time';
%                eInfo.reaction_time.Units = 'seconds';
%
%  'cInfo'     - [cell] cell array containing the field names for the
%                channel file in addition to the channel "name", "type" and
%                "unit" which are extracted from the channel location
%                structure. For example, to add the reference as an
%                additional column:
%                {'reference';
%                 'ref' }
%
%  'cInfoDesc' - [struct] structure describing additional or/and original
%                channel fields if you wish to redefine these.
%                cInfo.name.LongName = 'Channel name'; % change default description
%                cInfo.name.Description = 'Channel name';
%                cInfo.reference.LongName = 'Channel reference';
%                cInfo.reference.Description = 'Channel reference montage 10/20';
%

%  The following are automatically
%                populated using channel structure info ('eeg', 'ecg', 'emg', 'eog', 'trigger')
%                tInfo.ECGChannelCount = xxx;
%                tInfo.EEGChannelCount = xxx;
%                tInfo.EMGChannelCount = xxx;
%                tInfo.EOGChannelCount = xxx;
%                tInfo.MiscChannelCount = xxx;
%                tInfo.TriggerChannelCount = xxx;
%                tInfo.EEGReference = xxx;
%                tInfo.EpochLength = xxx;
%                tInfo.RecordingDuration = xxx;
%                tInfo.SamplingFrequency = xxx;
%                tInfo.TaskName = 'meditation';
%                However, they may be overwritten man.

function bids_format_eeglab(files, varargin)

if nargin < 1
    help bids_format_eeglab;
    return
end

opt = finputcheck(varargin, { 'ReferencesAndLinks' 'cell'   {}   { 'n/a' };
    'Name'               'string' {}   'n/a';
    'License'            'string' {}   'n/a';
    'targetdir' 'string'  {}   fullfile(pwd, 'bidsexport');
    'taskName'  'string'  {}   'meditation'; % TRY EMPTY TO SEE IF IT WORKS ***********************************
    'codefiles' 'cell'    {}   {};
    'stimuli'   'cell'    {}   {};
    'pInfo'     'cell'    {}   {};
    'eInfo'     'cell'    {}   {};
    'cInfo'     'cell'    {}   {};
    'gInfo'     'struct'  {}   struct([]);
    'tInfo'     'struct'  {}   struct([]);
    'pInfoDesc' 'struct'  {}   struct([]);
    'eInfoDesc' 'struct'  {}   struct([]);
    'cInfoDesc' 'struct'  {}   struct([]);
    'README'    'string'  {}   '';
    'CHANGES'   'string'  {}   '' }, 'bids_format_eeglab');
if isstr(opt), error(opt); end
if size(opt.stimuli,1) == 1 || size(opt.stimuli,1) == 1
    opt.stimuli = reshape(opt.stimuli, [2 length(opt.stimuli)/2])';
end

% deleting folder
fprintf('Exporting data to %s...\n', opt.targetdir);
if exist(opt.targetdir)
    disp('Deleting folder...')
    rmdir(opt.targetdir, 's');
end

disp('Creating sub-directories...')
mkdir( fullfile(opt.targetdir, 'code'));
mkdir( fullfile(opt.targetdir, 'stimuli'));
mkdir( fullfile(opt.targetdir, 'sourcedata'));

% write dataset info
% ------------------
gInfoFields = { 'ReferencesAndLinks' 'required' 'cell' { 'n/a' };
    'Name'               'required' 'char' '';
    'License'            'required' 'char' '';
    'BIDSVersion'        'required' 'char' '1.1.1' };
opt.gInfo = checkfields(opt.gInfo, gInfoFields, 'gInfo');
writejson(fullfile(opt.targetdir, 'dataset_description.json'), opt.gInfo);

% write participant information
% -----------------------------
participants = { 'participant_id' };
for iSubj=1:length(files)
    participants{iSubj+1, 1} = sprintf('sub-%3.3d', iSubj);
end
if ~isempty(opt.pInfo)
    if size(opt.pInfo,1) ~= length(participants)
        error(sprintf('Wrong number of participant (%d) in tInfo structure, should be %d based on the number of files', size(opt.pInfo,1)-1, length(files)));
    end
    participants(:,2:size(opt.pInfo,2)+1) = opt.pInfo;
end
writetsv(fullfile(opt.targetdir, 'participants.tsv'), participants);

% write participation field description
% -------------------------------------
descFields = { 'LongName'     'required' 'char' '';
    'Levels'       'optional' 'struct' {};
    'Description'  'optional' 'char' '';
    'Units'        'optional' 'char' '';
    'TermURL'      'optional' 'char' '' };
fields = fieldnames(opt.pInfoDesc);
if ~isempty(setdiff(fields, participants(1,:)))
    error('Some field names in the pInfoDec structure do not have a corresponding column name in pInfo');
end
fields = participants(1,:);
for iField = 1:length(fields)
    descFields{1,4} = fields{iField};
    if ~isfield(opt.pInfoDesc, fields{iField}), opt.pInfoDesc(1).(fields{iField}) = struct([]); end
    opt.pInfoDesc.(fields{iField}) = checkfields(opt.pInfoDesc.(fields{iField}), descFields, 'pInfoDesc');
end
writejson(fullfile(opt.targetdir, 'participants.json'), opt.pInfoDesc);

% write event file information
% ----------------------------
events = { 'onset' 'duration' 'event_sample' 'event_type' 'event_value' };std_2
fields = fieldnames(opt.eInfoDesc);
if ~isempty(setdiff(fields, events(1,:)))
    error('Some field names in the eInfoDesc structure do not have a corresponding column name in eInfo');
end
fields = events(1,:);
for iField = 1:length(fields)
    descFields{1,4} = fields{iField};
    if ~isfield(opt.eInfoDesc, fields{iField}), opt.eInfoDesc(1).(fields{iField}) = struct([]); end
    opt.eInfoDesc.(fields{iField}) = checkfields(opt.eInfoDesc.(fields{iField}), descFields, 'eInfoDesc');
end
writejson(fullfile(opt.targetdir, [ 'task-' opt.taskName '_events.json' ]), opt.eInfoDesc);

% Write README and CHANGES files
% ------------------------------
if ~isempty(opt.README)
    if ~exist(opt.README)
        fid = fopen(fullfile(opt.targetdir, 'README'), 'w');
        if fid == -1, error('Cannot write README file'); end
        fprintf(fid, '%s', opt.README);
        fclose(fid);
    else
        copyfile(opt.README, fullfile(opt.targetdir, 'README'));
    end
end
if ~isempty(opt.CHANGES)
    if ~exist(opt.CHANGES)
        fid = fopen(fullfile(opt.targetdir, 'CHANGES'), 'w');
        if fid == -1, error('Cannot write README file'); end
        fprintf(fid, '%s', opt.CHANGES);
        fclose(fid);
    else
        copyfile(opt.CHANGES, fullfile(opt.targetdir, 'CHANGES'));
    end
end

% write code files
% ----------------
if ~isempty(opt.codefiles)
    for iFile = 1:length(opt.codefiles)
        [~,fileName,Ext] = fileparts(opt.codefiles{iFile});
        if ~isempty(dir(opt.codefiles{iFile}))
            copyfile(opt.codefiles{iFile}, fullfile(opt.targetdir, 'code', [ fileName Ext ]));
        else
            fprintf('Warning: cannot find code file %s\n', opt.codefiles{iFile})
        end
    end
end

% write stimulus files
% --------------------
if ~isempty(opt.stimuli)
    if size(opt.stimuli,1) == 1, opt.stimuli = opt.stimuli'; end
    for iStim = 1:size(opt.stimuli,1)
        [~,fileName,Ext] = fileparts(opt.stimuli{iStim,1});
        if ~isempty(dir(opt.stimuli{iStim,1}))
            copyfile(opt.stimuli{iStim,1}, fullfile(opt.targetdir, 'stimuli', [ fileName Ext ]));
        else
            fprintf('Warning: cannot find stimulus file %s\n', opt.codefiles{iFile});
        end
    end
end

% check task info
% ---------------
opt.tInfo.TaskName = opt.taskName;
tInfoFields = {...
    'TaskName' 'REQUIRED' '' '';
    'TaskDescription' 'RECOMMENDED' '' '';
    'Instructions' 'RECOMMENDED' 'char' '';
    'CogAtlasID' 'RECOMMENDED' 'char' '';
    'CogPOID' 'RECOMMENDED' 'char' '';
    'InstitutionName' 'RECOMMENDED' 'char' '';
    'InstitutionAddress' 'RECOMMENDED' 'char' '';
    'InstitutionalDepartmentName' ' RECOMMENDED' 'char' '';
    'DeviceSerialNumber' 'RECOMMENDED' 'char' '';
    'SamplingFrequency' 'REQUIRED' '' '';
    'EEGChannelCount' 'REQUIRED' '' '';
    'EOGChannelCount' 'REQUIRED' '' 0;
    'ECGChannelCount' 'REQUIRED' '' 0;
    'EMGChannelCount' 'REQUIRED' '' 0;
    'EEGReference' 'REQUIRED' 'char' 'Unknown';
    'PowerLineFrequency' 'REQUIRED' '' 'Unknown';
    'EEGGround' 'RECOMMENDED ' 'char' '';
    'MiscChannelCount' ' OPTIONAL' 'char' '';
    'TriggerChannelCount' 'RECOMMENDED' 'char' '';
    'EEGPlacementScheme' 'RECOMMENDED' 'char' '';
    'Manufacturer' 'RECOMMENDED' 'char' '';
    'ManufacturersModelName' 'OPTIONAL' 'char' '';
    'CapManufacturer' 'RECOMMENDED' 'char' 'Unknown';
    'CapManufacturersModelName' 'OPTIONAL' 'char' '';
    'HardwareFilters' 'OPTIONAL' 'char' '';
    'SoftwareFilters' 'RECOMMENDED' 'char' '';
    'RecordingDuration' 'RECOMMENDED' 'char' 'n/a';
    'RecordingType' 'RECOMMENDED' 'char' '';
    'EpochLength' 'RECOMMENDED' '' 'n/a';
    'SoftwareVersions' 'RECOMMENDED' 'char' '';
    'SubjectArtefactDescription' 'OPTIONAL' 'char' '' };
opt.tInfo = checkfields(opt.tInfo, tInfoFields, 'tInfo');

% copy EEG files
% --------------
disp('Copying EEG files...')
for iSubj = 1:length(files)
    subjectStr    = sprintf('sub-%3.3d', iSubj);
    if iscell(files{iSubj})
        for iSess = 1:length(files{iSubj})
            copy_data_bids( files{iSubj}{iSess}, fullfile(opt.targetdir, subjectStr, sprintf('ses-%2.2d', iSess), 'eeg', [ subjectStr sprintf('_ses-%2.2d', iSess) '_task-' opt.taskName '_eeg' files{iSubj}{iSess}(end-3:end)]),opt.tInfo);
        end
    else
        copy_data_bids( files{iSubj}, fullfile(opt.targetdir, subjectStr, 'eeg', [ subjectStr '_task-' opt.taskName '_eeg' files{iSubj}{iSess}(end-3:end) ]),opt.tInfo);
    end
end

function copy_data_bids(fileIn, fileOut, tInfo)

    folderOut = fileparts(fileOut);
    if ~exist(folderOut)
        mkdir(folderOut);
    end
    if ~exist(fileOut)
    end
    % if BDF file anonymize records
    [~,~,ext] = fileparts(fileOut);
    if strcmpi(ext, '.bdf')
        fileIDIn  = fopen(fileIn,'rb','ieee-le');  % see sopen
        fileIDOut = fopen(fileOut,'wb','ieee-le');  % see sopen
        data = fread(fileIDIn, Inf);
        data(9:9+160-1) = ' ';
        fwrite(fileIDOut, data);
        fclose(fileIDIn);
        fclose(fileIDOut);
        tInfo.EEGReference = 'CMS/DRL';
        tInfo.Manufacturer = 'BIOSEMI';
    else
        copyfile(fileIn, fileOut);
    end

    EEG = pop_biosig(fileOut);
        
    % old conversion using data2bids
    %     addpath(fullfile(p, 'fieldtrip'));
    %     cfg = [];
    %     cfg.dataset                     = fileOut;
    %
    %     cfg.eeg.writesidecar            = 'replace';
    %     cfg.channels.writesidecar       = 'replace';
    %     cfg.events.writesidecar         = 'replace';
    %
    %     cfg.InstitutionName             = 'Paul Sabatier University';
    %     cfg.InstitutionalDepartmentName = 'Centre de Recherche Cerveau et Cognition';
    %     cfg.InstitutionAddress          = 'ISCT - Centre de Recherche Cerveau et Cognition, Place du Docteur Baylac, Pavillon Baudot, 31059 Toulouse';
    %
    %     % provide the long rescription of the task
    %     cfg.TaskName = taskName;
    %
    %     % these are EEG specific
    %     cfg.eeg.PowerLineFrequency      = 50;
    %     cfg.eeg.EEGReference            ='Na';
    %     %data2bids(cfg)

    insertEpoch = false;
    if EEG.trials > 1
        % get TLE events
        insertEpoch = true;
        eventlat = abs(eeg_point2lat( [ EEG.event.latency ], [ EEG.event.epoch ], EEG.srate, [EEG.xmin EEG.xmax]));
        indtle    = find(eventlat == 0);
        if length(indtle) < EEG.trials
            indtle    = find(eventlat < 0.02);
        end
        if length(indtle) ~= EEG.trials
            insertEpoch = false;
        end
    end

    % write event file information
    fid = fopen( [ fileOut(1:end-7) 'events.tsv' ], 'w');
    fprintf(fid, 'onset\tduration\tevent_sample\tevent_type\tevent_value\n');
    for iEvent = 1:length(EEG.event)
        onset = (EEG.event(iEvent).latency-1)/EEG.srate;

        % duration
        if isfield(EEG.event, 'duration') && ~isempty(EEG.event(iEvent).duration)
            duration = num2str(numEEG.event(iEvent).duration, '%1.10f') ;
        else
            duration = 'n/a';
        end

        % event type (which is the type of event - not the same as EEGLAB)
        eventType = 'STATUS';
        if insertEpoch
            if any(indtle == iEvent)
                eventType = 'Epoch';
            end
        end

        % event value
        if isstr(EEG.event(iEvent).type)
            eventValue = EEG.event(iEvent).type;
        else
            eventValue = num2str(EEG.event(iEvent).type);
        end

        fprintf(fid, '%1.10f\t%s\t%1.10f\t%s\t%s\n', onset, duration, EEG.event(iEvent).latency-1, eventType, eventValue);
    end
    fclose(fid);

    % write channel file information
    fid = fopen( [ fileOut(1:end-7) 'channels.tsv' ], 'w');
    if isempty(EEG.chanlocs)
        fprintf(fid, 'name\n');
        for iChan = 1:EEG.nbchan, printf(fid, 'E%d\n', iChan); end
    else
        fprintf(fid, 'name\ttype\tunits\n');

        for iChan = 1:EEG.nbchan
            if ~isfield(EEG.chanlocs, 'type') || isempty(EEG.chanlocs(iChan).type)
                type = 'n/a';
            else
                type = EEG.chanlocs(iChan).type;
            end
            if strcmpi(type, 'eeg')
                unit = 'uV';
            else
                unit = 'n/a';
            end

            fprintf(fid, '%s\t%s\t%s\n', EEG.chanlocs(iChan).labels, type, unit);
        end
    end
    fclose(fid);

    % write task information
    tInfo.EEGChannelCount = EEG.nbchan;
    if ~isfield(tInfo, 'EEGReference')
        tInfo.EEGReference    = EEG.ref;
    end
    if EEG.trials == 1
        tInfo.RecordingType = 'continuous';
    else
        tInfo.RecordingType = 'epoched';
        tInfo.EpochLength = EEG.pnts/EEG.srate;
    end
    tInfo.RecordingDuration = EEG.pnts/EEG.srate;
    tInfo.SamplingFrequency = EEG.srate;
    jsonStr = jsonencode(tInfo);
    fid = fopen( [fileOut(1:end-7) 'eeg.json' ], 'w');
    fprintf(fid, '%s', jsonStr);
    fclose(fid);

    % write channel information
%     cInfo.name.LongName = 'Channel name';
%     cInfo.name.Description = 'Channel name';
%     cInfo.type.LongName = 'Channel type';
%     cInfo.type.Description = 'Channel type';
%     cInfo.units.LongName = 'Channel unit';
%     cInfo.units.Description = 'Channel unit';
%     jsonStr = jsonencode(cInfo);
%     fid = fopen( [fileOut(1:end-7) 'channels.json' ], 'w');
%     fprintf(fid, '%s', jsonStr);
%     fclose(fid);

% check the fields for the structures
% -----------------------------------
function s = checkfields(s, f, structName)

    fields = fieldnames(s);
    diffFields = setdiff(fields, f(:,1)');
    if ~isempty(diffFields)
        error(sprintf('Invalid field name(s) %sfor structure %s', sprintf('%s ',diffFields{:}), structName));
    end
    for iRow = 1:size(f,1)
        if isempty(s) || ~isfield(s, f{iRow,1})
            if strcmpi(f{iRow,2}, 'required') % required or optional
                s = setfield(s, {1}, f{iRow,1}, f{iRow,4});
            end
        elseif ~isempty(f{iRow,3}) && ~isa(s.(f{iRow,1}), f{iRow,3})
            error(sprintf('Parameter %s.%s must be a %s', structName, f{iRow,1}, f{iRow,3}));
        end
    end

% write JSON file
% ---------------
function writejson(fileName, matlabStruct)
    jsonStr = jsonencode(matlabStruct);

    fid = fopen(fileName, 'w');
    if fid == -1, error('Cannot write file - make sure you have writing permission'); end
    fprintf(fid, '%s', jsonStr);
    fclose(fid);

% write TSV file
% --------------
function writetsv(fileName, matlabArray)
    fid = fopen(fileName, 'w');
    if fid == -1, error('Cannot write file - make sure you have writing permission'); end
    for iRow=1:size(matlabArray,1)
        for iCol=1:size(matlabArray,2)
            if isempty(matlabArray{iRow,iCol})
                disp('Empty value detected, replacing by n/a');
                fprintf(fid, 'n/a');
            elseif ischar(matlabArray{iRow,iCol})
                fprintf(fid, '%s', matlabArray{iRow,iCol});
            elseif isnumeric(matlabArray{iRow,iCol}) && rem(matlabArray{iRow,iCol},1) == 0
                fprintf(fid, '%d', matlabArray{iRow,iCol});
            elseif isnumeric(matlabArray{iRow,iCol})
                fprintf(fid, '%1.10f', matlabArray{iRow,iCol});
            else
                error('Table values can only be string or numerical values');
            end
            if iCol ~= size(matlabArray,2)
                fprintf(fid, '\t');
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
