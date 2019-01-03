function rename_brainvision_files(varargin)

% RENAME_BRAINVISION_FILES renames a BrainVision EEG dataset, which consists of a vhdr header
% file, vmrk marker file and a data file that usually has the extension dat, eeg or seg.
%
% Use as
%   rename_brainvision_files(oldname, newname, 'rmf', 'on')
% where both the old and the new filename should be strings corresponding to the
% header file, i.e. including the vhdr extension.
% 'rmf' option indicates to remove old files and can be turned 'on' of 'off' (default)
%
% See also http://www.fieldtriptoolbox.org/ and https://sccn.ucsd.edu/wiki/EEGLAB for
% open source software to process BrainVision EEG data.
%
% Robert Oostenveld https://gist.github.com/robertoostenveld/e31637a777c514bf1e86272e1092316e
% Cyril Pernet - fixed few bugs here and there https://gist.github.com/CPernet/e037df46e064ca83a49fb4c595d4566a

%% deal with inputs

rmf = 'off'; % by default do not delete files
if nargin >= 3
    if strcmpi(varargin{3},'rmf')
        if nargin == 3
            disp('no value associated to the remove file option, assumning off')
        else
            rmf = varargin{4};
        end
    else
        error(['unrecognized option argument ''' varargin{3} ''''])
    end
end

oldheaderfile = varargin{1};
newheaderfile = varargin{2};
clear varargin

% determine whether the file extensions should be in lower or upper case
if ~isempty(regexp(newheaderfile, 'VHDR$', 'once'))
  switchcase = @upper;
else
  switchcase = @lower;
end

% determine the filename without extension
[~, f, ~] = fileparts(newheaderfile);

%% do the renaming

% deal with the header file
assert(exist(oldheaderfile, 'file')~=0, 'the file %s does not exists', oldheaderfile);
assert(exist(newheaderfile, 'file')==0, 'the file %s already exists', newheaderfile);
fid1 = fopen(oldheaderfile, 'r'); % read old
fid2 = fopen(newheaderfile, 'w'); % write new

while ~feof(fid1)
  line = fgetl(fid1);
  if ~isempty(regexp(line, '^MarkerFile', 'once'))
    [~, rem] = strtok(line, '=');
    oldmarkerfile = rem(2:end);
    [~, ~, x] = fileparts(oldmarkerfile);
    newmarkerfile = [f switchcase(x)];
    line = sprintf('MarkerFile=%s', newmarkerfile);
  elseif ~isempty(regexp(line, '^DataFile', 'once'))
    [~, rem] = strtok(line, '=');
    olddatafile = rem(2:end);
    [~, ~, x] = fileparts(olddatafile);
    newdatafile = [f switchcase(x)];
    line = sprintf('DataFile=%s', newdatafile);
  end
  fprintf(fid2, '%s\r\n', line);
end
fclose(fid1);
fclose(fid2);

% deal with the marker file
if exist(oldmarkerfile, 'file') == 0
    [~,~,ext] = fileparts(oldmarkerfile);
    [~, ff, ~] = fileparts(oldheaderfile); % re-reading as sometimes weird names comes up
    if exist([ff ext], 'file')
        oldmarkerfile = [ff ext];
    else
        error('the file %s does not exists', oldmarkerfile);
    end
end
assert(exist(newmarkerfile, 'file')==0, 'the file %s already exists', newmarkerfile);
fid1 = fopen(oldmarkerfile, 'r');
fid2 = fopen(newmarkerfile, 'w');

while ~feof(fid1)
  line = fgetl(fid1);
  if ~isempty(regexp(line, '^HeaderFile', 'once'))
    [~, rem] = strtok(line, '=');
    oldheaderfile = rem(2:end);
    [~, ~, x] = fileparts(oldheaderfile);
    newheaderfile = [f switchcase(x)];
    line = sprintf('HeaderFile=%s', newheaderfile);
  elseif ~isempty(regexp(line, '^DataFile', 'once'))
    [~, rem] = strtok(line, '=');
    olddatafile = rem(2:end);
    [~, ~, x] = fileparts(olddatafile);
    newdatafile = [f switchcase(x)];
    line = sprintf('DataFile=%s', newdatafile);
  end
  fprintf(fid2, '%s\r\n', line);
end
fclose(fid1);
fclose(fid2);

% deal with the data file
if exist(olddatafile, 'file') == 0
    [~,~,ext] = fileparts(olddatafile);
    [~, ff, ~] = fileparts(oldheaderfile); % re-reading as sometimes weird names comes up
    if exist([ff ext], 'file')
        olddatafile = [ff ext];
    else
        error('the file %s does not exists', oldmarkerfile);
    end
end
assert(exist(newdatafile, 'file')==0, 'the file %s already exists', newdatafile);
status = copyfile(olddatafile, newdatafile);
if ~status
  error('failed to copy data from %s to %s', olddatafile, newdatafile);
end

%% delete old files *try* in case of user restriction
if strcmpi(rmf,'on')
    try delete(oldheaderfile); end
    try delete(oldmarkerfile); end
    try delete(olddatafile); end
end