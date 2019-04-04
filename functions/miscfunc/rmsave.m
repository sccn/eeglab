% rmsave() - return the RMS in each channel, epoch
%
% Usage:
%         >> ave = rmsave(data,frames);

% Copyright (C) Scott Makeig, CNL/Salk Institute, La Jolla, 9/98
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

function ave = rmsave(data,frames)

if nargin<1
  help rmsave
  return
end
if nargin<2
	frames = size(data,2); 
	data = reshape(data, size(data,1), size(data,2)*size(data,3));
end

chans = size(data,1);
datalength = size(data,2);
if rem(datalength,frames)
   fprintf('frames should divide data length.\n');
   return
end
if frames < 1
   fprintf('frames should be > 1.\n');
   return
end

epochs = datalength/frames;
ave = zeros(chans,epochs);
i=1;
while i<= epochs
  dat = matsel(data,frames,0,0,i);
  dat = dat.*dat;
  ave(:,i) = sqrt(mean(dat'))';
  i = i+1;
end

