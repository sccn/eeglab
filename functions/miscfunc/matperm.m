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

% 04-22-99  Adjusted for fixes and speed by Sigurd Enghoff & Tzyy-Ping Jung
% 01-25-02  Reformated help & license, added links -ad 

function [permx,indperm]= matperm(x,y,indx,indy,corr)

[m,n] = size(x);
[p,q] = size(y);
[ix,z] = size(indx);
[iy,z] = size(indy);
oldm = m;

errcode=0;
if  ix ~= iy | p ~= iy,
 fprintf('matperm: indx and indy must be column vectors, same height as y.\n');
 errcode=1
end;

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
  end;
end;
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
end;

return

function vals=sgn(data)

 vals = 2*(data>=0)-1;

return
