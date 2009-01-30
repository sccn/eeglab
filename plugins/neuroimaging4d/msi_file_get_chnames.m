function msi_file_get_index(fid, keyword, names)
%	MSI_FILE_GET_INDEX(fid, keyword)  Skip to a keyword in a file FID
%	and return the single valuea one based index array from a zero
%	based index array
%

frewind(fid);
test_line = sprintf('%s\t%%d', keyword);
index=0;
while feof(fid) == 0
    	line = fgetl(fid);
	matches = findstr(line, keyword);
	num = length(matches);
	if num > 0
	    % eat up keyword
	    [chopped,remainder] = strtok(line);
	    while(any(remainder));
		[chopped,remainder] = strtok(remainder, ',');
		names(index) = sscanf(chopped, '%s');
                index = index + 1;
	    end
	end
end

