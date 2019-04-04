% pcexpand() - expand data using Principal Component Analysis (PCA)
%              returns data expanded from a principal component subspace 
%                              [compare pcsquash()]
% Usage: 
%        After  >> [eigenvectors,eigenvalues,projections] = pcsquash(data,ncomps);
%        then   >> [expanded_data] = pcexpand(projections,eigenvectors,mean(data'));
%
% Inputs:
%    projections  = (comps,frames) each row is a component, each column a time point
%    eigenvectors = square matrix of (column) eigenvectors 
%    datameans    = vector of original data channel means
%
% Outputs: 
%    projections  = data projected back into the original data space
%                   size (chans=eigenvector_rows,frames)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 2000 
%
% See also: pcsquash(), svd()

% Copyright (C) 6-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 4-15-98  debugged -sm & t-pj
% 01-25-02 reformated help & license, added links -ad 

function [expanded_data]=pcexpand(PCAproj,EigenVectors,Datameans)

if nargin < 2
   help pcexpand
   return;
end

[ncomps,frames]=size(PCAproj);
[j,k]=size(EigenVectors);

if j ~= k 
   error('Wrong array input size (eigenvectors matrix not square');
end

if j < ncomps
   error('Wrong array input size (eigenvectors rows must be equal to projection matrix rows');
end

if size(Datameans,1) == 1,
    Datameans = Datameans';   % make a column vector
end
expanded_data = EigenVectors(:,1:ncomps)*PCAproj + Datameans*ones(1,frames);

