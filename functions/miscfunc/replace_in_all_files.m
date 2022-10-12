% Replace XXXX by XXXXX for better documentation in MATLAB
% allFiles = dir('*.m')
% for iFile = 1:length(allFiles)
%    replace_in_all_files(allFiles(iFile).name);
%    disp(allFiles(iFile).name);
% end

function replace_in_all_files(fileName)

fid = fopen(fileName, 'r');
if fid == -1
    error('File not found');
end

allStrs = {};
count   = 1;
expression1 = '[_A-Za-z0-9]+\(\)'; % + means at least 1
expression2 = '[A-Za-z0-9]\(\)';
while ~feof(fid)
    strTmp = fgetl(fid);

    if length(strTmp) > 0 && strTmp(1) == '%'
        pos1 = regexp(strTmp, expression1);
        pos2 = regexp(strTmp, expression2)+1;

        if length(pos1) ~= length(pos2) 
            error('Issue with length of regular expressions')
        end
        for iPos = length(pos1):-1:1
            strTmp(pos2(iPos):pos2(iPos)+1) = [];
            strTmp(pos1(iPos):pos2(iPos)-1) = upper(strTmp(pos1(iPos):pos2(iPos)-1));
        end
    end
    
    allStrs{count} = strTmp;
    count = count+1;
end
fclose(fid);    

% write new file
fid = fopen(fileName, 'w');
for iRow = 1:length(allStrs)
    fprintf(fid, '%s\n', allStrs{iRow});
end
fclose(fid);    

function [twoLines, content] = replace_txt(twoLines, pattern1, pattern2)

indB1 = strfind(twoLines, pattern1{1}); indPat1 = 1;
indB2 = strfind(twoLines, pattern2{1}); indPat2 = 1;
content = '';
if isempty(indB1), indB1 = strfind(twoLines, pattern1{2}); indPat1 = 2; end
if isempty(indB2), indB2 = strfind(twoLines, pattern2{2}); indPat2 = 2; end
if ~isempty(indB1) && ~isempty(indB2)
    content = twoLines(indB1+length(pattern1{indPat1}):indB2-1);
    fprintf('Found: %s\n', content);
    twoLines = [ twoLines(1:indB1-1) '[' content '](http://sccn.ucsd.edu/eeglab/locatefile.php?file=' content ')' twoLines(indB2+length(pattern2{indPat2}):end) ];
end
