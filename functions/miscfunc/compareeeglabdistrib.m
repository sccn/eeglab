% compare EEGLAB distributions

function compareeeglabdistrib(folder1, folder2)

fileList1 = getmfilelist(folder1);
fileList2 = getmfilelist(folder2);

% strip list from folders
for iFile = 1:length(fileList1), fileList1{iFile} = fileList1{iFile}(length(folder1)+1:end); end
for iFile = 1:length(fileList2), fileList2{iFile} = fileList2{iFile}(length(folder2)+1:end); end

for iFile1 = 1:length(fileList1)
    iFile2 = strmatch(fileList1{iFile1}, fileList2, 'exact');
    if isempty(iFile2)
        fprintf('Missing file %s in folder %s\n', fileList1{iFile1}, folder2);
    else
        dir1 = dir(fullfile(folder1, fileList1{iFile1}));
        dir2 = dir(fullfile(folder2, fileList2{iFile2}));
        if length(dir1) ~= 1 || length(dir2) ~= 1 
            error('Size error');
        end
        if dir1.bytes ~= dir2.bytes
            fprintf('Files %s differ\n', fileList2{iFile2});
            %fprintf('---------------------\n');
            %system([ '/usr/bin/diff ' fullfile(folder1, fileList1{iFile1}) ' ' fullfile(folder2, fileList2{iFile2}) ]);
        end
        fileList2(iFile2) = [];
    end
end

for iFile2 = 1:length(fileList2)
    fprintf('Missing file %s in folder %s\n', fileList2{iFile2}, folder1);
end
    


            
            
        
    
