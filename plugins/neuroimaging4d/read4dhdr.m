% read4dhdr() - read header information from 4d data file.
%
% Usage:
%   >> [head] = read4dhdr(fid)
%
% Input:
%   fid - file identifier of 4D set file
%
% Output:
%   head - structure containing header information.
%          Structure fields are:
%
% Author: Christian Wienbruch, 14 Mar 2004
%
%
% See also: read4d()

function head = read4dhdr(filename)

disp(filename);
if nargin < 1
    help read4dhdr;
    return;
end;
set_fid=fopen(filename,'rt');

head.fileformat=0;
key_count = msi_file_find_keyword(set_fid, 'MSI.Format: SHORT');
if (key_count == 1)
  head.fileformat=1;
  disp('File Format is SHORT');
end
key_count = msi_file_find_keyword(set_fid, 'MSI.Format: LONG');
if (key_count == 1)
  head.fileformat=2;
  disp('File Format is LONG');
end
key_count = msi_file_find_keyword(set_fid, 'MSI.Format: FLOAT');
if (key_count == 1)
  head.fileformat=3;
  disp('File Format is FLOAT');
end
key_count = msi_file_find_keyword(set_fid, 'MSI.Format: DOUBLE');
if (key_count == 1)
  head.fileformat=4;
  disp('File Format is DOUBLE');
end

key_count = msi_file_find_keyword(set_fid, 'MSI.TotalChannels:');
if (key_count == 1)
    head.TotalChannels = msi_file_get_long(set_fid, 'MSI.TotalChannels:');
end

key_count = msi_file_find_keyword(set_fid, 'MSI.TotalEpochs:');
if (key_count == 1)
    head.TotalEpochs = msi_file_get_long(set_fid, 'MSI.TotalEpochs:');
end

key_count = msi_file_find_keyword(set_fid, 'MSI.SlicesPerEpoch:');
if (key_count == 1)
    head.TotalSlices = msi_file_get_long(set_fid, 'MSI.SlicesPerEpoch:');
end

key_count = msi_file_find_keyword(set_fid, 'MSI.MegChanCount:');
if (key_count == 1)
    head.MegChanCount = msi_file_get_long(set_fid, 'MSI.MegChanCount:');
end

key_count = msi_file_find_keyword(set_fid, 'MSI.SamplePeriod:');
if (key_count == 1)
    head.SamplePeriod = msi_file_get_float(set_fid, 'MSI.SamplePeriod:');
end

key_count = msi_file_find_keyword(set_fid, 'MSI.FirstLatency:');
if (key_count == 1)
    head.FirstLatency = msi_file_get_float(set_fid, 'MSI.FirstLatency:');
end
key_count = msi_file_find_keyword(set_fid, 'MSI.TotalEvents:');
 if (key_count == 1)
     head.TotalEvents = msi_file_get_long(set_fid, 'MSI.TotalEvents:');
 end
key_count = msi_file_find_keyword(set_fid, 'MSI.TriggerIndex:');
if (key_count == 1)
    head.TriggerIndex = msi_file_get_long(set_fid, 'MSI.TriggerIndex:');
end
key_count = msi_file_find_keyword(set_fid, 'MSI.ResponseIndex:');
if (key_count == 1)
    head.ResponseIndex = msi_file_get_long(set_fid, 'MSI.ResponseIndex:');
end

fprintf(1, 'TotalChannels = %d\n', head.TotalChannels);
fprintf(1, ' = %d\n', head.TotalSlices);
fprintf(1, 'SamplePeriod is %f\n', head.SamplePeriod);
fprintf(1, 'FirstLatency is %f\n', head.FirstLatency);
head.Duration = head.SamplePeriod * head.TotalSlices;
head.samp_rate = 1./head.SamplePeriod;
head.LastSample = head.Duration + head.FirstLatency;
fprintf(1, 'Duration is %f\n', head.Duration);
fprintf(1, 'LastSample is %f\n', head.LastSample);
head.segsamps = head.TotalSlices ;

head.MegIndexArray = msi_file_get_index(set_fid, 'MSI.MegChanIndex:');
head.EegIndexArray = msi_file_get_index(set_fid, 'MSI.EegChanIndex:');
head.Events        = msi_file_get_index(set_fid, 'MSI.Events:');
head.segments      = head.Events;
head.EventCodes    = msi_file_get_index(set_fid, 'MSI.EventCodes:');

key_count = msi_file_find_keyword(set_fid, 'MSI.Meg_Position_Information.Begin:');
if (key_count == 1)
    fprintf(1, 'reading coil locations\n');
    [head.Meg_pos, head.Meg_dir] = msi_file_MEG_loc(set_fid);
end

key_count = msi_file_find_keyword(set_fid, 'MSI.Eeg_Position_Information.Begin:');
if (key_count == 1)
    fprintf(1, 'reading coil locations\n');
    [head.Eeg_pos, head.Eeg_dir] = msi_file_EEG_loc(set_fid);
end

key_count = msi_file_find_keyword(set_fid, 'MSI.Der_Position_Information.Begin:');
if (key_count == 1)
    fprintf(1, 'reading coil locations\n');
    [head.Der_pos, head.Der_dir] = msi_file_DER_loc(set_fid);
end

fclose(set_fid);
