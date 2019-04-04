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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
