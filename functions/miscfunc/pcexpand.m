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

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 6-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: not supported by cvs2svn $

% 4-15-98  debugged -sm & t-pj
% 01-25-02 reformated help & license, added links -ad 

function [expanded_data]=pcexpand(PCAproj,EigenVectors,Datameans)

if nargin < 2
   help pcexpand
end

[ncomps,frames]=size(PCAproj);
[j,k]=size(EigenVectors);

if j < ncomps | nargin < 2 | j ~= k 
   help pcexpand
end

if size(Datameans,1) == 1,
    Datameans = Datameans';   % make a column vector
end
expanded_data = EigenVectors(:,1:ncomps)*PCAproj; + Datameans*ones(1,frames);

