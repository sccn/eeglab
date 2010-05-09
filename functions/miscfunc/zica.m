% zica() - Z-transform of ICA activations; useful for studying component SNR
%
% Usage: >> [zact,basesd,maz,mazc,mazf] = zica(activations,frames,baseframes)
%
% Inputs:
%   activations - activations matrix produced by runica()
%   frames      - frames per epoch {0|default ->  length(activations)}
%   baseframes  - vector of frames in z-defining baseline period {default frames}
% 
% Outputs:
%   zact        - activations z-scaled and reorded in reverse order of max abs
%   basesd      - standard deviations in each activation row (reverse ordered)
%   maz         - maximum absolute z-value for each activation row (rev ordered)
%   mazc        - component indices of the reverse-sorted max abs z-values
%                 (this is the act -> zact reordering)
%   mazf        - frame indices of the max abs z-values (reverse ordered)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 2-25-98 
%
% See also: runica()

function [zact,basesd,maxabsz,maxc,maxabszf] = zica(activations,frames,baseframes)

% Copyright (C) 2-25-98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 3-2-98 added frames variable -sm
% 1-25-01 put revsort subfunction in the core of the zica program -ad
% 01-25-02 reformated help & license, added link -ad 

if nargin<1
   help zica
   return
end

[chans,framestot] = size(activations);

if nargin < 3
   baseframes = 0;
end
if nargin < 2
   frames = 0;
end
if frames == 0
   frames = framestot;
end
if baseframes == 0
   baseframes = 1:frames
end
epochs = floor(framestot/frames);
if frames*epochs ~= framestot
   fprintf('zica(): indicated frames does not divide data length.\n');
   return
end

if length(baseframes) < 3
  fprintf('\n  zica() - too few baseframes (%d).\n',length(baseframes));
  help zica
  return
end

if min(baseframes) < 1 | max(baseframes) > frames
  fprintf('\n  zica() - baseframes out of range.\n');
  help zica
  return
end

baselength = length(baseframes);
baseact = zeros(epochs*baselength,chans);
for e=1:epochs
 baseact((e-1)*baselength+1:e*baselength,:) = ...
                        matsel(activations,frames,baseframes,0,e)';
end
basesd = sqrt(covary(baseact));
zact = activations./(basesd'*ones(1,framestot));
[maxabsz,maxabszf] = sort(abs(zact'));
maxabsz  = maxabsz(frames,:);
maxabszf = maxabszf(frames,:);

%
% reorder outputs in reverse order of max abs z
[maxabsz,maxc] = revsort(maxabsz);
zact = zact(maxc,:);
basesd = basesd(maxc);
maxabszf = maxabszf(maxc);

% revsort  - reverse sort columns (biggest 1st, ...)
function [out,i] = revsort(in)

if size(in,1) == 1
   in = in'; % make column vector
end

[out,i] = sort(in);
out = out(size(in,1):-1:1,:);
  i = i(size(in,1):-1:1,:);
return;
