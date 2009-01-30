function result = msi_file_get_long(fid, keyword)
%	MSI_FILE_GET_LONG(fid, keyword)	Skip to a keyword in a file FID
%	and return the single value
%

frewind(fid);
test_line = sprintf('%s\t%%d', keyword);
while feof(fid) == 0
    	line = fgetl(fid);
	matches = findstr(line, keyword);
	num = length(matches);
	if num > 0
	    result = sscanf(line, test_line, 1);
	    % fprintf(1, 'I found a %d\n', result);
	end
end

