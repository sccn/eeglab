function channelloc_to_tvs

% start at the root and convert the channel_loc.txt files for BIDS
% SPM12 dependency

folders = dir('sub*')
for f=7:length(folders)
    cd(folders(f).name)
    file = dir('*.txt')
    fileID = fopen(file.name);
    C = textscan(fileID,'%d %f32 %f32 %f32 %s %s');
    fclose(fileID);
    
    % do the electrodes.tsv
    clear name
    for i=1:74
        name(i,:) = cell2mat(C{5}(i));
    end
    x = C{2}(1:74);
    y = C{3}(1:74);
    z = C{4}(1:74);
    T = table(name,x,y,z);
    writetable(T,[pwd filesep 'eeg' filesep folders(f).name '_task-facerecognition_electrodes.txt'],'Delimiter','\t')
    movefile([pwd filesep 'eeg' filesep folders(f).name '_task-facerecognition_electrodes.txt'],[pwd filesep 'eeg' filesep folders(f).name '_task-facerecognition_electrodes.tsv'])
    
    % do the electrodes.json
    json = struct('EEGCoordinateSystem','T1w', ...
        'EEGCoordinateUnits','mm','IntendedFor',['anat' filesep folders(f).name '_T1w.nii.gz']);
    spm_jsonwrite(['eeg' filesep folders(f).name '_task-facerecognition_electrodes.json'],json,struct('indent','  '))
    
    % do the fid.json
    json = struct(...
        'EEGCoordinateSystem','T1w', ...
        'EEGCoordinateUnits', 'mm', ...
        'AnatomicalMRICoordinateSystem', 'T1w', ...
        'AnatomicalMRICoordinateUnits', 'mm', ...
        'LandmarkCoordinates', ' ', ...
        'LandmarkCoordinateSystem','T1w', ...
        'LandmarkCoordinateUnits','mm');
    
    clear name
    for i=75:77
        tmp = cell2mat(C{5}(i));
        name(i-74,:) = tmp(1:3);
    end
        json.LandmarkCoordinates = ['{"' name(1,:) '":[' num2str(double([C{2}(75),C{3}(75),C{4}(75)])) '], "' ...
       name(2,:) '":[' num2str(double([C{2}(76),C{3}(76),C{4}(76)])) '], "' ...
        name(3,:) '":[' num2str(double([C{2}(77),C{3}(77),C{4}(77)])) ']}'];

    
    spm_jsonwrite(['eeg' filesep folders(f).name '_task-facerecognition_fid.json'],json,struct('indent','  '))
    
    % clean up
    delete(file.name)
    cd ..
end