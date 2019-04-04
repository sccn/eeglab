% matperm() - transpose and sign rows of x to match y (run after matcorr() )
%
% Usage: >> [permx indperm] = matperm(x,y,indx,indy,corr);
%
% Inputs:
%   x     = first input matrix 
%   y     = matrix with same number of columns as x
%   indx  = column containing row indices for x (from matcorr())
%   indy  = column containing row indices for y (from matcorr())
%   corr  = column of correlations between indexed rows of x,y (from matcorr())
%           (used only for its signs, +/-) 
% Outputs:
%   permx   = the matrix x permuted and signed according to (indx, indy,corr) 
%             to best match y. Rows of 0s added to x to match size of y if nec.
%   indperm = permutation index turning x into y;
%
% Authors: Scott Makeig, Sigurd Enghoff & Tzyy-Ping Jung 
%          SCCN/INC/UCSD, La Jolla, 2000 

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

% 04-22-99  Adjusted for fixes and speed by Sigurd Enghoff & Tzyy-Ping Jung
% 01-25-02  Reformated help & license, added links -ad 

function [permx,indperm]= matperm(x,y,indx,indy,corr)

[m,n] = size(x);
[p,q] = size(y);
[ix,z] = size(indx);
[iy,z] = size(indy);
oldm = m;

errcode=0;
if  ix ~= iy || p ~= iy,
 fprintf('matperm: indx and indy must be column vectors, same height as y.\n');
 errcode=1
end

if n~=q,
   fprintf('matperm(): two matrices must be same number of columns.\n');
   errcode=2;
else
  if m<p,
  		x = [x;zeros(p-m,n)];	% add rows to x to match height of y
  		p=m;
  elseif p<m,
  		y = [y;zeros(m-p,n)];	% add rows to y to match height of x
  		m=p;
  end
end
if errcode==0,
%
% Return the row permutation of matrix x most correlated with matrix y:
%  plus the resulting permutation index
%
  indperm = [1:length(indx)]';	% column vector [1 2 ...nrows]   
  permx  = x(indx,:); 
  indperm = indperm(indx,:);
  ydni(indy) = 1:length(indy);
  permx = permx(ydni,:);% put x in y row-order
  indperm = indperm(ydni,:);
  permx = permx.*(sgn(corr(ydni))*ones(1,size(permx,2))); 
								 % make x signs agree with y
  permx = permx(1:oldm,:);		 % throw out bottom rows if 
								 % they were added to match y
  indperm = indperm(1:oldm,:);
end

return

function vals=sgn(data)

 vals = 2*(data>=0)-1;

return
