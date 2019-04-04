% std_findoutlierclust() - determine whether an outlier cluster already exists 
%                 for a specified cluster. If so, return the outlier cluster index.  
%                 If not, return zero. This helper function is called by 
%                 pop_clustedit(), std_moveoutlier(), std_renameclust(). 

% Copyright (C) 2006 UCDS
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

function outlier_clust = std_findoutlierclust(STUDY,clust)
outlier_clust = 0;
outlier_name = [ 'Outliers ' STUDY.cluster(clust).name ' ' ];
for k = 1:length(STUDY.cluster)
    if strncmpi(outlier_name,STUDY.cluster(k).name,length(outlier_name))
        outlier_clust = k;
    end
end
