function [eeg_loc, eeg_dir] = msi_file_EEG_loc(set_fid)
%function [eeg_loc, eeg_dir] = msi_file_EEG_loc(set_fid, chans)
%	MSI_FILE_GET_LONG(set_fid, chans)	Skip to a keyword in a file FID
%	and return the single float value
%

key_count = msi_file_find_keyword(set_fid, 'MSI.EegChanCount:');
if (key_count == 1)
    EegChanCount = msi_file_get_long(set_fid, 'MSI.EegChanCount:');
end

frewind(set_fid);
keyword = 'MSI.Eeg_Position_Information.Begin:';

while feof(set_fid) == 0
    	line = fgetl(set_fid);
	matches = findstr(line, keyword);
	num = length(matches);
	if num > 0
	    fprintf(1, 'I found a %s\n', line);
	    for (ii=1:EegChanCount)
		line = fgetl(set_fid);
		[name, count, errmsg, nextindex] = sscanf(line, '%s\t', 1);
		[token, nextline] = strtok(line, '	');
		[bigresult, count, errmsg, nextindex] = sscanf(nextline, '\t%f\t%f\t%f\t%f\t%f\t%f', 6);
		loc(1) = bigresult(1);
		loc(2) = bigresult(2);
		loc(3) = bigresult(3);
		dir(1) = bigresult(4);
		dir(2) = bigresult(5);
		dir(3) = bigresult(6);
		%fprintf(1, 'chan: %s\n',  name);
		%fprintf(1, '  %f %f %f   ',  loc(1), loc(2), loc(3));
		%fprintf(1, '  %f %f %f \n',  dir(1), dir(2), dir(3));
		%fprintf(1, 'chan: %s (%f,%f,%f) (%f,%f,%f)\n',  name, loc[1], loc[2], loc[3], dir[1], dir[2], dir[3]);
		outloc(ii,1) = loc(1);
		outloc(ii,2) = loc(2);
		outloc(ii,3) = loc(3);

		outdir(ii,1) = dir(1);
		outdir(ii,2) = dir(2);
		outdir(ii,3) = dir(3);
	    end
eeg_loc = outloc;
eeg_dir = outdir;
	end
end

