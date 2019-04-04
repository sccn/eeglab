% pcsquash() - compress data using Principal Component Analysis (PCA)
%              into a principal component subspace.  To project back 
%              into the original channel space, use pcexpand()
%
% Usage: 
%        >> [eigenvectors,eigenvalues] = pcsquash(data,ncomps);
%        >> [eigenvectors,eigenvalues,compressed,datamean] ...
%                                                    = pcsquash(data,ncomps);
%
% Inputs:
%   data    = (chans,frames) each row is a channel, each column a time point
%   ncomps  = numbers of components to retain
%
% Outputs: 
%   eigenvectors = square matrix of (column) eigenvectors 
%   eigenvalues  = vector of associated eigenvalues 
%   compressed   = data compressed into space of the ncomps eigenvectors
%                  with largest eigenvalues (ncomps,frames)
%                  Note that >> compressed = eigenvectors(:,1:ncomps)'*data;
%   datamean     = input data channel (row) means (used internally)
%
% Author: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 6-97 
%
% See also: pcexpand(), svd()

% Copyright (C) 2000 Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, 
% scott@sccn.ucsd.edu
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

function [EigenVectors,EigenValues,Compressed,Datamean]=pcsquash(matrix,ncomps)

if nargin < 1 
   help pcsquash
   return
end
if nargin < 2
   ncomps = 0;
end
if ncomps == 0
  ncomps = size(matrix,1);
end
if ncomps < 1
   help pcsquash
   return
end

data = matrix';                    % transpose data
[n,p]=size(data);                  % now p chans,n time points
if ncomps > p
   fprintf('pcsquash(): components must be <= number of data rows (%d).\n',p);
   return;
end

Datamean = mean(data,1);  % remove column (channel) means
data = data-ones(n,1)*Datamean;    % remove column (channel) means
out=data'*data/n;
[V,D] = eig(out);                  % get eigenvectors/eigenvalues
diag(D);
[eigenval,index] = sort(diag(D));
index=rot90(rot90(index));
EigenValues=rot90(rot90(eigenval))';
EigenVectors=V(:,index);

if nargout >= 3
   Compressed = EigenVectors(:,1:ncomps)'*data';
end
