function [volt] = eeg_load(FILENAME,SIZE,PRECISION)

% eeg_load - Load binary EEG data
%
% Can be used to easily load files written by 'eeg_write'.  Otherwise,
% can be used to load a binary matrix data file.
%
% Useage: [volt] = eeg_load(filename,SIZE,PRECISION)
%
% where:    'filename' is the full path + fileprefix + filextension
%           'SIZE' is [M,N] rows by columns; if omitted, SIZE can
%           be read from any data files written by 'eeg_write' but
%           otherwise all data is read into an 1D array.
%           'PRECISION' is the data type (default is 'double').
%           See matlab fread command for more info on size/precision.
%
% comment:  A simple wrapper for 'fread' matlab function, useful in 
%           other scripts. This is complementary to 'eeg_write'.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2000, Darren.Weber_at_radiology.ucsf.edu
%           04/2002, Darren.Weber_at_radiology.ucsf.edu
%                    added machineformat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('PRECISION','var'), PRECISION = 'double'; end

file = char(FILENAME);    

[path,name,ext] = fileparts(file);
file = fullfile(path,[name ext]);

if ~isequal(exist(file),2),
    lookfile = which(file);
    if isempty(lookfile),
        msg = sprintf('Cannot locate %s\n', file);
        error(msg);
    else
        file = lookfile;
    end
end

fprintf('EEG_LOAD: Reading:\n... %s : ', file);

fid = fopen(file,'r','ieee-le');

if ~exist('SIZE','var'),
    [SIZE, COUNT] = fread(fid,[1,2],'int');
end

[volt, COUNT] = fread(fid,SIZE,PRECISION);

S = size(volt);     fprintf('%d rows, %d cols\n', S(1), S(2) );
