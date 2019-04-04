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
