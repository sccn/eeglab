% neural_net() - computes clusters using Matlab Neural Net toolbox. 
%        Alternative clustering algorithm to kmeans().
%        This is a helper function called from pop_clust(). 

% Copyright (C) 2006 UCSD
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

function [IDX,C] = neural_net(clustdata,clus_num)

nmin = min(clustdata);
nmax = max(clustdata);
net = newc([nmin ;nmax].',clus_num);
net = train(net,(clustdata).');
Y = sim(net,(clustdata).');
IDX = vec2ind(Y);
C = zeros(clus_num,size(clustdata,2));
for k = 1:clus_num
    C(k,:) = sum(clustdata(find(IDX == k),:));
end
