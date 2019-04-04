% icavar() - project ICA component activations through the ICA weight matrices 
%            to reconstitute the observed data using selected ICA components. 
%            Returns time course of variance on scalp for each component.
%
% Usage: >> [srcvar] = icavar(data,weights,sphere,compnums);
%
% Inputs:
%   data     - data matrix returned by runica()
%   weights  - weight matrix returned by runica()
%   sphere   - sphere matrix returned by runica() (default|0 -> eye(ncomps))
%   compnums - list of component numbers (default|0 -> all)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 11-30-1996 
%
% See also: runica()

% Copyright (C) 11-30-1996 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 04-03-97  made more efficient -sm
% 05-20-97  debugged variance calculation and made it more efficient -sm
% 06-07-97  changed order of args to conform to runica -sm
% 07-23-97  do not add back datamean to projections -sm
% 12-23-97  added compnums -sm
% 02-25-98  changed activations input to data -sm
% 07-20-98  use pinv() for inverting non-square weights -sm
% 01-25-02 reformated help & license, added links -ad 

function [srcvar] = icavar(data,weights,sphere,compnums)

if nargin<2
   help icavar
   return
end

if size(weights,2) ~= size(sphere,1) || size(sphere,2) ~= size(data,1)
   fprintf('icavar() - sizes of weights, sphere, and data incompatible.\n')
   whos data weights sphere
   return
end
activations = weights*sphere*data;
[ncomps,frames] = size(activations);
[wr,chans]      = size(weights);    % Note: ncomps may < data chans

if nargin<4
    compnums = 0;
end
if compnums == 0,
    compnums = 1:ncomps;
end
srcvar = zeros(length(compnums),frames);

if nargin < 3
   sphere = 0;
end
if sphere == 0,
   sphere = eye(ncomps);
end

project = weights*sphere;

if wr~=ncomps,
  fprintf('icavar: Number of rows in activations and weights must be equal.\n');
  return
end

% if wr<chans,
  % fprintf('Filling out a square projection matrix with small noise.\n');
  % for r=1:chans-wr,
    % project = [project;0.00001*randn(1,chans)];
  % end
% end
% project = inv(project);                        % invert projection matrix

if wr<chans
  project = pinv(project);
else
  project = inv(project);
end

fprintf('Projecting ICA component ');
nout = length(compnums);
for i=1:nout                                 % for each component 
  comp = compnums(i);
  fprintf('%d ',comp); 
  projdata = project(:,comp)*activations(comp,:);       % reproject eeg data 
  srcvar(i,:) = sum(projdata.*projdata)/(chans-1);% compute variance at each time point
end
fprintf('\n');
