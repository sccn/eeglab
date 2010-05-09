% blockave() - make block average of concatenated data sets of same size 
%              Each data set is assumed to be of size (chans,frames).  
% Usage: 
%        >> aveout = blockave(data,frames) 
%        >> aveout = blockave(data,frames,epochs,weights) 
% Inputs: 
%         data    = data matrix of size (chans,frames*epochs) 
%         frames  = columns per data epoch
%         epochs  = vector of epochs to average {default|0 -> all}
%         weights = vector of epoch weights {default|0 -> equal}
%                   non-equal weights are normalized internally to sum=1
% Output:
%         aveout = averaged data of size (chans, frames)
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

% 5-14-98 improved error checking, allowed 2 args -sm
% 01-25-02 reformated help & license -ad 

function ave = blockave(data,frames,epochs,weights)

if nargin < 2
  help blockave
  return
end
if nargin<3
  epochs = 0;
end
if nargin<4
  weights = nan;
end
if nargin>4
  help blockave
end

[chans,framestot]=size(data);
nepochs = floor(framestot/frames);
if nepochs*frames ~= framestot
   fprintf('blockave(): frames arg does not divide data length,\n');
   return
end
if max(epochs) > nepochs
   help blockave
   return
end
if min(epochs) < 1
   if size(epochs,1)>1 | size(epochs,2) > 1
     help blockave
     return
   else
     epochs = 1:nepochs;
   end
end
if ~isnan(weights)
  if length(weights) ~= nepochs
      fprintf(...
  'blockave(): number of weights (%d) must match number of epochs (%d).\n',...
                               length(weights),                  nepochs);
  end
  weights = weights/sum(weights); % normalize
end

ave = zeros(chans,frames);
epochs = epochs(:)'; % make epochs a row vector
n = 1;
for e=epochs
   if isnan(weights)
     ave = ave + data(:,(e-1)*frames+1:e*frames);
   else
     ave = ave + weights(n)*data(:,(e-1)*frames+1:e*frames);
     n = n+1;
   end
end
ave = ave/length(epochs);
