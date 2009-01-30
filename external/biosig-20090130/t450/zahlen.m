% M1-1 Zufallszahlen aus dem Intervall [1,N]
% The function [g] = zahlen(M,N) is a subfunction in the main-function fdp.m.
%For help use the main-function.
%
%
%Copyright (C) 2006 by Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
%Institute of Medical Statistics, Computer Sciences and Documantation
%University of Jena
%This work was supported by DFG Project VO 683/2-1
%This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%%***
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.
%
%--------------------------------------------------------------------------



function [g]=zahlen(M,N);
if  M<N/2

g=ceil(N*rand(1,M));

g=unique(g);

l=size(g,2);
while l<M
g1=ceil(N*rand(1,M-l));
g=union(g1,g);l=size(g,2);
end;

else
M=N-M;
g=ceil(N*rand(1,M));

g=unique(g);

l=size(g,2);
while l<M
g1=ceil(N*rand(1,M-l));
g=union(g1,g);l=size(g,2);
end;
g=setdiff([1:N],g);
end;
