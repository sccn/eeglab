% matsel() - select rows, columns, and epochs from given multi-epoch data matrix
%
% Usage:
%       >> [dataout] = matsel(data,frames,framelist);
%       >> [dataout] = matsel(data,frames,framelist,chanlist);
%       >> [dataout] = matsel(data,frames,framelist,chanlist,epochlist);
%
% Inputs:
%   data      - input data matrix (chans,frames*epochs) 
%   frames    - frames (data columns) per epoch (0 -> frames*epochs)
%   framelist - list of frames per epoch to select (0 -> 1:frames)
%   chanlist  - list of chans to select (0 -> 1:chans)
%   epochlist - list of epochs to select (0 -> 1:epochs)
%
% Note: The size of dataout is (length(chanlist), length(framelist)*length(epochlist))
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 5-21-96 

% Copyright (C) 5-21-96 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 5-25-96 added chanlist, epochlist -sm
% 10-05-97 added out of bounds tests for chanlist and framelist -sm
% 02-04-00 truncate to epochs*frames if necessary -sm
% 01-25-02 reformated help & license -ad 

function [dataout] = matsel(data,frames,framelist,chanlist,epochlist)

if nargin<1
  help matsel
  return
end

[chans framestot] = size(data);
if isempty(data)
  fprintf('matsel(): empty data matrix!?\n')
  help matsel
  return
end

if nargin < 5,
  epochlist = 0;
end
if nargin < 4,
  chanlist = 0;
end
if nargin < 3,
  fprintf('matsel(): needs at least 3 arguments.\n\n');
  return
end

if frames == 0,
  frames = framestot;
end

if framelist == 0,
  framelist = [1:frames];
end

framesout = length(framelist);

if isempty(chanlist) | chanlist == 0,
  chanlist = [1:chans];
end

chansout = length(chanlist);
epochs = floor(framestot/frames);

if epochs*frames ~= framestot
    fprintf('matsel(): data length %d was not a multiple of %d frames.\n',...
                          framestot,frames);
    data = data(:,1:epochs*frames);
end

if isempty(epochlist) | epochlist == 0,
    epochlist = [1:epochs];
end
epochsout = length(epochlist);

if max(epochlist)>epochs 
      fprintf('matsel() error: max index in epochlist (%d) > epochs in data (%d)\n',...
                        max(epochlist),epochs);
      return
end

if max(framelist)>frames 
      fprintf('matsel() error: max index in framelist (%d) > frames per epoch (%d)\n',...
                        max(framelist),frames);
      return
end
    
if min(framelist)<1
      fprintf('matsel() error: framelist min (%d) < 1\n', min(framelist));
      return
end

if min(epochlist)<1
      fprintf('matsel() error: epochlist min (%d) < 1\n', min(epochlist));
      return
end
   
if max(chanlist)>chans
      fprintf('matsel() error: chanlist max (%d) > chans (%d)\n',...
                        max(chanlist),chans);
      return
end
    
if min(chanlist)<1
      fprintf('matsel() error: chanlist min (%d) <1\n',...
                        min(chanlist));
      return
end
    
dataout = zeros(chansout,framesout*epochsout);
for e=1:epochsout
    dataout(:,framesout*(e-1)+1:framesout*e) = ...
                data(chanlist,framelist+(epochlist(e)-1)*frames); 
end

