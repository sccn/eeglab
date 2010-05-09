% fastranc.m -  adaptive noise cancellation for fmrib_fastr.m
%   [out,y]=fmribanc(refs,d,N,mu)
%
%   Inputs:
%       refs: Reference data
%       d:   Input data
%       N:   filter order (length)
%       mu:  Step size
%
%   Outputs:
%       out: error signal (signal after noise cancellation)
%       y: filtered out noise
%
%   If a binary version of this program is present in same directory
%   it will be run instead.  The binary version is MUCH fastr.  However,
%   if it is not present or unavailable for your platform, you can use
%   the matlab compiler to compile the C Mex-File  fastranc.c by typing
%   >> mex fastranc.c
%   from your command line after you CD into the fmrib1.0 plugin directory.  
%   This will create the appropriate binary file for your platform 
%   (if you have the MATLAB compiler). You can check what binaries are
%   available by looking at the extention of the fastranc.* files:
%   fastranc.mexa64   -  Linux 64 bit (compiled on Xeon EM64T but might
%                          be the same for athalon 64)
%   fastranc.mexaxp   -  Alpha 64 processor,
%   fastranc.mexglx   -  Linux 32 bit processor.
%   fastranc.dll      -  Windows.  
%
%   Check your matlab help for more info on compiling MEX files.
%
%
%   Author:  Rami K. Niazy
%   
%   Copyright (c) 2004 University of Oxford

% Copyright (C) 2004 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
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

function [out,y]=fastranc(refs,d,N,mu)

nargchk(4,4,nargin);
refs=refs(:);
d=d(:);
mANC=length(d);
if length(d)~=length(refs)
error('Reference and Input data must be of the same length','fmribanc() error!');
end

W=zeros(N+1,1);
r=flipud([0;refs(1:N)]);
out=zeros(mANC,1);
y=zeros(mANC,1);

for E=N+1:mANC

    %-----------------calc------------------------------------------

    r=[refs(E);r(1:end-1)];
    y(E)=sum(W.*r);
    out(E)=d(E)-y(E);
    W=W+2*mu*out(E)*r;

    %---------------------------------------------------------------
end
return;
