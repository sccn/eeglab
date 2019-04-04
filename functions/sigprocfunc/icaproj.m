% icaproj() - project ICA component activations through the
%             associated weight matrices to reconstitute the
%             observed data using only the selected ICA components.
% Usage:
%   >> [icaprojdata] = icaproj(data,weights,compindex,[datameans],chansout);
%
% Inputs:
%   data        - data matrix (chans, frames*epochs)
%   weights     - unmixing weight matrix (e.g., weights*sphere from runica())
%   compindex   - vector of ICA component indices to project
%
% Optional inputs:
%   datamean    - Optional ICA row means (for each epoch) from runica() 
%                 {default 0 -> distribute data offsets among the ICA components}
%   chansout    - Optional vector of channel indices to output {default: all}
%
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 11-30-96 
%
% See also: icavar(), runica()

% Copyright (C) 11-30-96 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 11-30-96 Scott Makeig CNL / Salk Institute, La Jolla as icaproject.m
% 12-24-96 added define for verbose -sm (V1.3)
% 2-11-97  use outer product math when only one component in compindex -sm
% 3-11-97  remove row means instead of grand mean -sm
% 3-19-97  used datamean argument instead of frames/baseframes -sm
% 4-03-97  changed name to icaproj() -sm
% 6-07-97  changed order of args to conform to runica -sm
% 6-10-97  fixed data mean handling -sm
% 6-23-97  trying pseudo-inverse for non-square weights -sm
% 7-23-97  distributed baseline offset (if any) among the activations -sm
% 10-31-97 removed errcode -sm
% 12-19-00 removed sphere, shifted order of args -sm
% 05-29-01 added chansout, made more efficient -sm
% 01-25-02 reformated help & license, added links -ad 

function  [icaprojdata] = icaproj(data,weights,compindex,datamean,chansout)

verbose = 0;   % default no-verbose

if nargin<3   % need 3 args
    help icaproj
    return
end
if nargin<4
  datamean = 0;  % default
end
if isempty(datamean)
  datamean = 0;  % default
end

[chans,framestot] = size(data);
if nargin<5
  chansout = [];
end
if isempty(chansout) || chansout(1) == 0
  chansout = 1:chans;
end
if min(chansout)<1 || max(chansout)> chans
  fprintf('icaproj(): chansout variable out of 1:chans range.\n')
  return
end

[mchans,epochs] = size(datamean);
frames = floor(framestot/epochs);
if epochs*frames ~= framestot || frames < 1,
    fprintf(...
        'icaproj(): frames (%d) does not divide data length (%d)\n.',...
                frames,framestot);
    return
end

[ncomps,cols] = size(compindex);
if cols>1,
  if ncomps==1,    % if row vector, 
      compindex = compindex';    % make col vector
      ncomps = cols;
  else
      fprintf('icaproj(): compindex must be a vector\n');
      return
  end
end
if ncomps>chans,
  fprintf('icaproj(): compindex must have <= %d entries\n',chans);
  return
end
for a=1:ncomps-1
  for b=a+1:ncomps
      if compindex(a,1)==compindex(b,1),
            fprintf('icaproj(): component index repeated in compindex\n.');
          return
      end
  end
end
for a=1:ncomps
    if compindex(a)>chans || compindex(a)<1
      fprintf('icaproj(): component index %d out of range!\n',compindex(a));
      return
      break
    end
end

if nargin<4
  datamean = 0; % default
end

if datamean ~= 0,
  %
  % Remove row means, e.g. those subtracted prior to ICA training by runica()
  %
  if verbose==1,
     fprintf('Removing data means of each channel and data epoch...\n');
  end
  for e=1:epochs
      data(:,(e-1)*frames+1:e*frames) = ...
         data(:,(e-1)*frames+1:e*frames) - datamean(:,e)*ones(1,frames);
  end
end
if verbose == 1
  fprintf('Final input data range: %g to %g\n', ...
                min(min(data)),max(max(data)));
end

if size(weights,1) == size(weights,2)
  iweights    = inv(weights);             % inverse weight matrix
else
  iweights    = pinv(weights);            % pseudo-inverse weight matrix
end

activations = weights(compindex,:)*data;  % activation waveforms
                                          % at desired components
if ncomps==1,  % compute outer product only for single component projection
  if verbose==1,
    fprintf('icaproj(): Projecting data for ICA component %d\n',...
                     compindex(1));
  end
else % if ncomps > 1
  if verbose==1,
    fprintf('icaproj(): Projecting data for ICA components ');
    if ncomps<32
     for n=1:ncomps % for each included component 
      fprintf('%d ',compindex(n));
     end
    else
      fprintf('specified.');
    end
    fprintf('\n');         % copy selected activations
  end
end

icaprojdata = iweights(chansout,compindex)*activations; 
% reconstitute selected scalp data channels from selected ICA components
