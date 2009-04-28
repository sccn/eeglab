function eeg_write(file,data,PRECISION)

% eeg_write - Write matrix of EEG data into binary file
%
% Useage: eeg_write(file,data,PRECISION)
%
% where:    'file' is the path + filename, eg "c:\data.dat"
%           'data' is a matlab variable
%           'PRECISION' is a binary format ('double' default, see fread)
%
% comment:  Uses the fwrite method to output both 'size(data)' and
%           data. The size(data) integers are written first and used
%           by 'eeg_load' to determine the size of the data matrix.
%           All files are 'ieee-le'.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('PRECISION','var'), PRECISION = 'double'; end

[path,name,ext] = fileparts(file);
file = fullfile(path,[name ext]);

fid = fopen(file,'w','ieee-le');

% Write magicN data file size
magicN = size(data);
count = fwrite(fid,magicN,'int');
    
% Write the rest of the data
count = fwrite(fid,data,PRECISION);  fclose(fid);

fprintf('EEG_WRITE: Wrote %d bytes to:\t%s\n', count, file);

return
