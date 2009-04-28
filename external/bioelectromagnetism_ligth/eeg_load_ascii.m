function [volt,var] = eeg_load_ascii(filename)

% eeg_load_ascii - Load ascii EEG data matrix
%
% Assumes a simple matrix format of data, compatible with matlab
% 'load' command, with no strings (eg, labels) in the data file.
% Looks for a voltage file and an associated .var variance file.
%
% Useage: [volt,var] = eeg_load_ascii(filename)
%
% where:    filename is the full path + fileprefix + filextension
%
% comment:  A very simple function, useful in other scripts.
%           Could be developed further to include matrix orientation
%           checks and/or reorientation options.  Currently a simple
%           wrapper for the 'load' matlab function.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no implied or express warranties
% History:  07/2000, Darren.Weber_at_radiology.ucsf.edu
%           10/2003, Darren.Weber_at_radiology.ucsf.edu
%                    added variance file handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegversion = '$Revision: 1.1 $';
fprintf('\nEEG_LOAD_ASCII [v %s]\n',eegversion(11:15)); tic;

file = char(filename);

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

if ~isequal(exist(file),2),
  lookfile = which(file);
  if isempty(lookfile),
    msg = sprintf('...cannot locate %s\n', filename);
    error(msg);
  else
    file = lookfile;
  end
end

fprintf('...reading %s : ', [name ext]);
volt = load(file);
s = size(volt);
fprintf('%d rows, %d cols\n', s(1), s(2));


% Attempt to load associated variance file
varfile = fullfile(path,strcat(name,'.var'));
if isequal(exist(varfile),2),
  fprintf('...reading %s : ', [name ext]);
  var = load(file);
  s = size(var);
  fprintf('%d rows, %d cols\n', s(1), s(2));
else
  lookfile = which(varfile);
  if isempty(lookfile),
    fprintf('...cannot locate ASCII variance : %s\n', strcat(name,'.var'));
    var = [];
  else
    varfile = lookfile;
    fprintf('...reading %s : ', [name ext]);
    var = load(file);
    s = size(var);
    fprintf('%d rows, %d cols\n', s(1), s(2));
  end
end

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
