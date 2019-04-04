% rmsave() - return the RMS in each channel, epoch
%
% Usage:
%         >> ave = rmsave(data,frames);

% Copyright (C) Scott Makeig, CNL/Salk Institute, La Jolla, 9/98
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

