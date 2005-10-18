% cls_findoutlierclust() - find if an outlier cluster already exist for a 
%       specific cluster. If so return the outlier cluster index, if not 
%       return zero.
%
%       This is a helper function called from pop_clustedit(), 
%       cls_moveoutlier() & cls_renameclust(). 

function outlier_clust = cls_findoutlierclust(STUDY,clust)
outlier_clust = 0;
outlier_name = [ 'Outliers ' STUDY.cluster(clust).name ' ' ];
for k = 1:length(STUDY.cluster)
    if strncmpi(outlier_name,STUDY.cluster(k).name,length(outlier_name))
        outlier_clust = k;
    end
end
