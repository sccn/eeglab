% this function detects potential malware in the current folder and subfolders
%
% Author: A. Delorme, Cotober 2013

function detectmalware(currentFolder);

if nargin < 1
    currentFolder = pwd;
end;

folderContent = dir(currentFolder);
folderContent = { folderContent.name };

malwareStrings = { 'eval(' 'evalin(' 'evalc(' 'delete(' 'movefile(' 'rmdir' 'mkdir' 'copyfile' 'system(' '!' };

for iFile = 1:length(folderContent)
    
    currentFile = folderContent{iFile};
    if length(currentFile) > 2 && strcmpi(currentFile(end-1:end), '.m')
        
        fid = fopen(fullfile(currentFolder, currentFile), 'r');
        countLine = 0;
        prevstr = '';
        while ~feof(fid)
            str = fgetl(fid);
            countLine = countLine+1;
        
            if length(str) > 1 && str(1) ~= '%'
                res = cellfun(@(x)~isempty(findstr(x, str)), malwareStrings(1:end-1));
                if str(1) == '!', res(end+1) = 1; end;
                if any(res)
                    pos = find(res); pos = pos(1);
                    disp('************************************')
                    fprintf('Potential malware command detected containing "%s" in\n %s line %d\n', malwareStrings{pos}, fullfile(currentFolder, currentFile), countLine);
                    fprintf('%d: %s\n%d: %s\n%d: %s\n', countLine-1, prevstr, countLine, str, countLine+1, fgetl(fid));
                    countLine = countLine+1;
                end;
            end;
            
            prevstr = str;
        end;
        fclose(fid);
        
    elseif exist(currentFile) == 7 && ~strcmpi(currentFile, '..') && ~strcmpi(currentFile, '.')
        detectmalware(fullfile(currentFolder, currentFile));
    end;
end;
