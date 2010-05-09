% abspeak() - find peak amps/latencies in each row of a single-epoch data matrix 
%
% Usage:
%        >> [amps,frames,signs] = abspeak(data);
%        >> [amps,frames,signs] = abspeak(data,frames);
%
% Inputs:
%   data   - single-epoch data matrix
%   frames - frames per epoch in data {default|0 -> whole data}
%
% Outputs:
%   amps   - amplitude array
%   frames - array of indexes
%   sign   - sign array
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1998 

% Copyright (C) 1998 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 3-9-98 added frames arg -sm
% 01-25-02 reformated help & license -ad 

function [amps,frames,signs]= abspeak(data,fepoch);

if nargin < 1
  help abspeak
  return
end

[chans,ftot] = size(data);
if nargin < 2
  fepoch = 0;
end
if fepoch == 0
  fepoch = ftot;
end
epochs = floor(ftot/fepoch);
if fepoch*epochs ~= ftot
   fprintf('abspeak(): frames arg does not divide data length.\n')
   return
end

amps   = zeros(chans,epochs);
frames = zeros(chans,epochs);
signs  = zeros(chans,epochs);

for e=1:epochs
  for c=1:chans
    dat = abs(matsel(data,fepoch,0,c,e))';
    [sdat,si]   = sort(dat);
    amps(c,e)   = dat(si(fepoch));           % abs value at peak
    frames(c,e) = si(fepoch);                % peak frame
    signs(c,e)  = sign(data(c,frames(c,e))); % sign at peak
  end
end
