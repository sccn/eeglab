function [file] = read_ctf_ascii(filename);

% READ_CTF_ASCII reads general data from an CTF configuration file
%
% The file should be formatted like
%    Group
%    {
%      item1 : value1a value1b value1c 
%      item2 : value2a value2b value2c 
%      item3 : value3a value3b value3c 
%      item4 : value4a value4b value4c 
%    }
%
% This fileformat structure is used in 
%   params.avg
%   default.hdm
%   multiSphere.hdm
%   processing.cfg
% and maybe for other files as well.

% Copyright (C) 2003, Robert Oostenveld
% 
% $Log: not supported by cvs2svn $
% Revision 1.1  2005/12/06 06:24:20  psdlw
% Alternative functions from the FieldTrip package, which is now released under GPL (so I assume these functions can be committed to the sourceforge cvs)
%
% Revision 1.2  2003/04/17 12:38:08  roberto
% *** empty log message ***
%
% Revision 1.1  2003/03/24 12:30:42  roberto
% new implementation
%

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

line = '';
while ischar(line)
  line = cleanline(fgetl(fid));
  if isempty(line) | line==-1
    continue
  end

  % the line is not empty, which means that we have encountered a chunck of information
  subline = cleanline(fgetl(fid));	% read the {
  subline = cleanline(fgetl(fid));	% read the first item
  while isempty(findstr(subline, '}'))
    if ~isempty(subline)
      [item, value] = strtok(subline, ':');
      value(1) = ' ';			% remove the :
      value  = deblank2(value);
      item   = deblank2(item);
      warning off
      if isempty(str2num(value))
        % the value appears to be a string
        eval(sprintf('file.%s.%s = [ ''%s'' ];', line, item, value));
      else
        % the value appears to be a number or a list of numbers
        eval(sprintf('file.%s.%s = [ %s ];', line, item, value));
      end
      warning on
    end
    subline = cleanline(fgetl(fid));	% read the first item
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line = cleanline(line)
  if isempty(line) | line==-1
    return
  end
  comment = findstr(line, '//');
  if ~isempty(comment)
    line(min(comment):end) = ' ';
  end
  line = deblank2(line);

