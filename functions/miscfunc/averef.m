% averef() - convert common-reference EEG data to average reference
%
% Usage:
%   >> data = averef(data);
%
% Notes:
%   Input and output data are (chans,frames*epochs) matrices of EEG 
%   or MEG data
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1999 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: not supported by cvs2svn $

% 12/16/99 Corrected denomiator on the suggestion of Ian Nimmo-Smith, Cambridge UK
% 01-25-02 reformated help & license -ad 

function data = averef(data)

if nargin<1
  help averef
  return
end
chans = size(data,1);
if chans < 2 
  help averef
  return
end

% avematrix = eye(chans)-ones(chans)*1/chans;
% data = avematrix*data; % implement as a matrix multiply
% else (faster?)

  data = data - ones(chans,1)*sum(data)/chans;

