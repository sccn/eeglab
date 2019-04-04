% averef() - convert common-reference EEG data to average reference
%            Note that this old function is not being used in EEGLAB. The
%            function used by EEGLAB is reref().
%
% Usage:
%   >> data = averef(data);
%   >> [data_out W_out S_out meandata] = averef(data,W);
%
% Inputs:
%   data - 2D data matrix (chans,frames*epochs) 
%   W    - ICA weight matrix
%
% Outputs:
%   data_out - Input data converted to average reference.
%   W_out    - ICA weight matrix converted to average reference
%   S_out    - ICA sphere matrix converted to eye()
%   meandata - (1,dataframes) mean removed from each data frame (point)
%
% Note: If 2 args, also converts the weight matrix W to average reference:
%         If ica_act = W*data, then data = inv(W)*ica_act; 
%         If R*data is the average-referenced data, 
%         R*data=(R*inv(W))*ica_act and W_out = inv(R*inv(W));
%         The average-reference ICA maps are the columns of inv(W_out).
%
% Authors: Scott Makeig and Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 
%
% See also: reref()

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 12/16/99 Corrected denomiator on the suggestion of Ian Nimmo-Smith, Cambridge UK
% 01-25-02 reformated help & license -ad 

function [data, W, S, meandata] = averef(data, W, S)

if nargin<1
  help averef
  return
end
chans = size(data,1);
if chans < 2 
  help averef
  return
end

% avematrix = eye(chans)-ones(chans)*1/chans;
% data = avematrix*data; % implement as a matrix multiply
% else (faster?)

meandata = sum(data)/chans;
data = data - ones(chans,1)*meandata;

% treat optional ica parameters
if nargin == 2
	winv  = pinv(W);
    size1 = size(winv,1);
	avematrix = eye(size1)-ones(size1)*1/size1;
	W = pinv(avematrix*winv);
end
if nargin >= 3
	winv = pinv(W*S);
    size1 = size(winv,1);
	avematrix = eye(size1)-ones(size1)*1/size1;
	W = pinv(avematrix*winv);
	S = eye(chans);
end
