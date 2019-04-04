% posact() - Make runica() activations all RMS-positive.
%            Adjust weights and inverse weight matrix accordingly.
%
% Usage: >> [actout,winvout,weightsout] = posact(data,weights,sphere) 
%
% Inputs:
%    data        = runica() input data
%    weights     = runica() weights
%    sphere      = runica() sphere {default|0 -> eye()}
%
% Outputs:
%    actout      = activations reoriented to be RMS-positive
%    winvout     = inv(weights*sphere) reoriented to match actout
%    weightsout  = weights reoriented to match actout (sphere unchanged)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 11/97 

% Copyright (C) 11/97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license, added links -ad 

function [actout,winvout,weightsout] = posact(data,weights,sphere)

if nargin < 2
   help posact
   return
end
if nargin < 3
   sphere = 0;
end

[chans,frames]=size(data);
[r,c]=size(weights);
if sphere == 0
  sphere = eye(chans);
end
[sr,sc] = size(sphere);
if sc~= chans
   fprintf('posact(): Sizes of sphere and data do not agree.\n')
   return
elseif c~=sr
   fprintf('posact(): Sizes of weights and sphere do not agree.\n')
   return
end

activations = weights*sphere*data;

if r==c
  winv = inv(weights*sphere);
else
  winv = pinv(weights*sphere);
end

[rows,cols] = size(activations);

actout = activations;
winvout = winv;

fprintf('Inverting negative activations: ');
for r=1:rows,
        pos = find(activations(r,:)>=0);
        posrms = sqrt(sum(activations(r,pos).*activations(r,pos))/length(pos));
        neg = find(activations(r,:)<0);
        negrms = sqrt(sum(activations(r,neg).*activations(r,neg))/length(neg));
        if negrms>posrms
            fprintf('-');   
            actout(r,:) = -1*activations(r,:);
            winvout(:,r) = -1*winv(:,r);
        end
        fprintf('%d ',r);
end
fprintf('\n');

if nargout>2
  if r==c,
    weightsout = inv(winvout);
  else
    weightsout = pinv(winvout);
  end
  if nargin>2 % if sphere submitted
    weightsout = weightsout*inv(sphere); % separate out the sphering
  end
end
