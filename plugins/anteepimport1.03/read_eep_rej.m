function [rej] = read_eep_rej(fn);

% READ_EEP_REJ reads rejection marks from an EEProbe *.rej file
%
% This function returns a Nx2 matrix with the begin and end latency
% of N rejection marks. The latency is in miliseconds.
%
% rej = read_eep_rej(filename)
%
% An EEProbe rejection file is formatted like
%   0.0000-0.3640
%   2.4373-3.5471
%   ... 
% where rejection begin and end are given in seconds. This function 
% converts the latency in miliseconds.
%
% Author: Robert Oostenveld, Aalborg University, Denmark, 11 March 2003
%
% See also READ_EEP_CNT, READ_EEP_TRG, READ_EEP_AVR
%

% Copyright (C) 2002, Robert Oostenveld
%                     Aalborg University, Denmark
%                     http://www.smi.auc.dk/~roberto/
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

rej = [];

fid = fopen(fn, 'rb');
if fid<0
   return 
end
while ~feof(fid)
  tmp = fscanf(fid, '%f-%f', 2);
  if ~isempty(tmp)
    rej = [rej; tmp'];
  end
end

% convert to ms
rej = 1000*rej;

fclose(fid);  
