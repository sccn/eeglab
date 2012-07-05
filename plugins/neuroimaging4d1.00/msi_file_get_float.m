function result = msi_file_get_float(fid, keyword)
%	MSI_FILE_GET_LONG(fid, keyword)	Skip to a keyword in a file FID
%	and return the single float value
%

frewind(fid);
test_line = sprintf('%s\t%%f', keyword);
while feof(fid) == 0
    	line = fgetl(fid);
	matches = findstr(line, keyword);
	num = length(matches);
	if num > 0
	    result = sscanf(line, test_line, 1);
	    % fprintf(1, 'I found a %f\n', result);
	end
end

