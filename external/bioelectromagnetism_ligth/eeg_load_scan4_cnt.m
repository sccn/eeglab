function cnt = eeg_load_scan4_cnt(filename,channels,range)

% eeg_load_scan4_cnt - Load a scan4.1+ CNT file
% 
% cnt = eeg_load_scan4_cnt('filename',channels,range)
% 
% filename   -  Filename string
% channels   -  see eeg_load_scan4_cnt_data,
%               default is 'all' channels (sensors/electrodes)
% range      -  see eeg_load_scan4_cnt_data,
%               1x2 start/stop sample points
%               eg, [1 1000] loads first 1000 points
%
%               NOTE: default range is [1 1000], not
%               all of the CNT file!  If you want it all,
%               specify 'all'
%
% cnt        -  struct, see eeg_load_scan4_cnt_data
%
% See also: eeg_load_scan4_cnt_data, eeg_load_scan4_cnt_event
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/2002, Darren.Weber_at_radiology.ucsf.edu
%                    created
%           10/2003, Darren.Weber_at_radiology.ucsf.edu
%                    modified default range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('filename','var'),
    msg = sprintf('EEG_LOAD_SCAN4_CNT: No filename provided\n');
    error(msg);
end

if ~exist('channels','var'), channels = 'all'; end

if ~exist('range','var'), range = [1 1000]; end

[path,name,ext] = fileparts(filename);
filename = fullfile(path,[name ext]);

if ~isequal(exist(filename),2),
    lookfile = which(filename);
    if isempty(lookfile),
        msg = sprintf('Cannot locate %s\n', filename);
        error(msg);
    else
        filename = lookfile;
    end
end

fprintf('EEG_LOAD_SCAN4_CNT: Reading file: %s\n',filename);

fid = fopen(filename,'r','ieee-le');

[cnt.volt,cnt.srate,cnt.numSamples,cnt.labels,cnt.events] = eeg_load_scan4_cnt_data(fid,channels,range);

cnt.path = [path,filesep];
cnt.file = [name,ext];

fclose(fid);

return
