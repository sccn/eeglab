function result = msi_file_find_keyword(fid, keyword)
%	MSI_FILE_FIND_KEYWORD(fid, keyword)	Skip to a keyword in a file FID
%
%

result = 0;
frewind(fid);
while feof(fid) == 0
    	line = fgetl(fid);
	matches = findstr(line, keyword);
	num = length(matches);
	if num > 0
	    result = result + 1;
	    % fprintf(1, 'I found a %d\n', result);
	end
end

