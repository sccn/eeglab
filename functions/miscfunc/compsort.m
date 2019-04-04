% compsort() - reorder ICA components, first largest to smallest by the size 
%             of their maximum variance in the single-component projections, 
%             then (if specified) the nlargest component projections are 
%             reordered by the (within-epoch) time point at which they reach 
%             their max variance.
%
% Usage:
%   >> [windex,maxvar,maxframe,maxepoch,maxmap] ...
%                          = compsort(data,weights,sphere,datamean, ...
%                                                 frames,nlargest,chanlist);
% Inputs:
%   data     = (chans,frames*epochs) the input data set decomposed by runica()
%   weights  = ica weight matrix returned by runica() 
%   sphere   = sphering matrix returned by runica() 
%
% Optional:
%   datamean = means removed from each row*epoch in runica() 
%              (Note: 0 -> input data means are distributed among components
%              : 1 -> input data means are removed from the components (default))
%   frames   = frames per epoch in data (0 -> all)
%   nlargest = number of largest ICA components to order by latency to max var
%              (other returned in reverse order of variance) (0 -> none)
%   chanlist = list of channel numbers to sort on (0 -> all)
%
% Outputs:
%   windex   = permuted order of rows in output weights (1-chans)
%   maxvar   = maximum variance of projected components (in perm. order)
%   maxframe = frame of maximum variance (1-frames) (in perm. order)
%   maxepoch = epoch number of frame of maxvar (1-nepochs) (in perm. order)
%   maxmap   = projected scalp map at max (in perm. order)
%
% Authors: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1996 
%
% See also: runica()

% Copyright (C) 1996 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 11-30-96 as compsort.m
% 12-19-96 moved def of epochs below default frames def -sm
% 02-18-97 added chanlist -sm
% 03-11-97 added default chanlist=0 -> all chans, fixed datamean and var() -sm
% 03-17-97 made nlargest default -> nlargest=0 -sm
% 04-03-97 fixed problems and removed permweights, permcomps outputs -sm 
% 04-04-97 shortened name to compsort() -sm
% 06-05-97 corrected variance computation -sm
% 06-07-97 changed order of args to conform to runica -sm
% 06-10-97 fixed recent bug in maxvar order -sm
% 07-23-97  made datamean==0 distribute means among components -sm
% 08-04-97  added datamean=1 option -sm
% 01-25-02 reformated help & license, added links -ad 

function [windex,maxvar,maxframe,maxepoch,maxmap] = compsort(data,weights,sphere,datamean,frames,nlargest,chanlist)

if nargin<3,
    fprintf('compsort(): needs at least three arguments.\n\n');
    help compsort
    return
end

[chans,framestot] = size(data);
if framestot==0,
    fprintf('Gcompsort(): cannot process an empty data array.\n');
    return
end
[srows,scols] = size(sphere);
[wrows,wcols] = size(weights);

if nargin<7,
    chanlist = [1:chans];
end
if chanlist==0,
    chanlist = [1:chans];
end
if length(chanlist)~=chans && wrows<chans,  % if non-square weight matrix
    fprintf('compsort(): chanlist not allowed with non-square weights.\n');
    return
end
if size(chanlist,1)>1
    chanlist = chanlist'; % make a row vector
end
if size(chanlist,1)>1
    fprintf('compsort(): chanlist must be a vector.\n');
    return
end
if nargin<6,
    nlargest=0;
end
if nargin<5,
    frames = 0;
end
if nargin<4,
   datamean = 1;
end
if frames ==0,
    frames = framestot;
end
epochs = framestot/frames;

% activations = (wrows,wcols)x(srows,scols)x(chans,framestot)
if chans ~= scols || srows ~= wcols,
   fprintf('compsort(): input data dimensions do not match.\n');
   fprintf(...
 ' i.e., Either chans %d ~= sphere cols %d or sphere rows %d ~= weights cols %d\n',...
                             chans,scols,srows,wcols);
   return
end
if wrows ~= chans && nlargest ~= 0 && nlargest ~= wrows,
   fprintf(...
 'compsort(): cannot project components back to order by size - nchans ~= ncomponents.\n');
   return
end
if floor(framestot/frames)*frames ~= framestot
   fprintf(...
 'compsort(): input data frames does not divide data length.\n');
   return
end

if nlargest > wrows,
   fprintf(...
 'compsort(): there are only %d rows in the weight matrix.\n',wrows);
   return
end

if epochs ~= floor(epochs),
   fprintf(...
 'compsort(): input frames does not subdivide data length.\n');
   return
end

%
%%%%%%%%%%%%%%%%%%%% Reorder weight matrix %%%%%%%%%%%%%%%%%%%%%
%
if datamean == 1,
    data = data - mean(data')'*ones(1,framestot); % remove channel means

elseif datamean~=0,                               % remove given means
    if size(datamean,2) ~= epochs || size(datamean,1) ~= chans,
        fprintf('compsort(): datamean must be 0, 1, or (chans,epochs)\n');
        return
    end
    for e=1:epochs
        data(:,(e-1)*frames+1:e*frames) = data(:,(e-1)*frames+1:e*frames)...
             - datamean(:,e)*ones(1,frames);
    end
end     % compute mean data matrix inherited from runica()

comps   = weights*sphere*data; % Note: distributes means if datamean==0

maxvar   = zeros(wrows,1);     % size of the projections
maxframe = zeros(wrows,1);     % frame of the abs(max) projection
maxepoch = zeros(wrows,1);     % epoch of the abs(max) projection

maxmap = zeros(wrows,chans);   % leave 0s unless weights is square
if chans==wrows,
  icainv = inv(weights*sphere);
  fprintf('Computing projected variance for all %d components:\n',wrows);
  for s=1:wrows
    fprintf('%d ',s);          % construct single-component data matrix
                               % project to scalp
    compproj   = icainv(:,s)*comps(s,:); 
    compv = zeros(frames*epochs,1);
    compv = (sum(compproj(chanlist,:).*compproj(chanlist,:)))' ...
                      /(length(chanlist)-1);
    [m,mi] = max(compv);       % find max variance
    maxvar(s)  = m;
    maxframe(s)=rem(mi,frames);    % count from beginning of each epoch!
    maxepoch(s)=floor(mi/frames)+1;% record epoch number of max
    if maxframe(s)==0,             % if max var is in last frame . . .
        maxframe(s) = frames;
        maxepoch(s) = maxepoch(s)-1;
    end
    maxmap(:,s) = compproj(:,mi);  % record scalp projection at max(var(proj))
  end

else % weight matrix is non-square, sort components by latency only
  fprintf('compsort() - non-square weights - finding max latencies.\n');
  for s=1:wrows
    compv = comps(s).*comps(s,:)'; % get variance analogues at each time point
    [m,mi] = max(abs(compv)');     % find abs(max)'s
    maxvar(s) = m;
    maxframe(s)=rem(mi,frames);    % count from beginning of each epoch!
    maxepoch(s)=floor(mi/frames)+1;% record epoch number of max
    if maxframe(s)==0,             % if max var is in last frame . . .
        maxframe(s) = frames;
        maxepoch(s) = maxepoch(s)-1;
    end
  end
end
fprintf('\n');

[maxvar windex] = sort(maxvar');% sort components by relative size
windex   = windex(wrows:-1:1)';    % reverse order
maxvar   = maxvar(wrows:-1:1)';    % make each returned vector a column vector

% Note: maxvar reordered by sort() above
maxframe = maxframe(windex);       % reorder maxframes
maxepoch = maxepoch(windex);       % reorder maxepoch
maxmap   = maxmap(:,windex);       % reorder maxmap columns

if nlargest>0,
  fprintf('Ordering largest %d components by frame of abs(max): ',nlargest);
  [m mfi] = sort(maxframe(1:nlargest)); % sort by frame order 
  windex(1:nlargest) = windex(mfi);
  maxframe(1:nlargest) = m;
  maxvar(1:nlargest) = maxvar(mfi);
  maxepoch(1:nlargest) = maxepoch(mfi);
  maxmap(:,1:nlargest) = maxmap(:,mfi); % reorder largest maxmap columns
  fprintf('\n');
end
