function result = msi_file_get_index(fid, keyword)
%	MSI_FILE_GET_INDEX(fid, keyword)  Skip to a keyword in a file FID
%	and return the single valuea one based index array from a zero
%	based index array
%

result(1) = 0;
frewind(fid);
index = 1;
test_line = sprintf('%s\t%%d', keyword);
while feof(fid) == 0
    	line = fgetl(fid);
	matches = findstr(line, keyword);
	num = length(matches);
	if num > 0
	    % eat up keyword
	    [chopped,remainder] = strtok(line);
	    while(any(remainder));
		[chopped,remainder] = strtok(remainder, ',');
		temp_value = sscanf(chopped, '%d');
		result(index) = temp_value;
		index = index + 1;
	    end
	end
end

