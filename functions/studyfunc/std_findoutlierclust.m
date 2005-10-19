% cls_findoutlierclust() - determine whether an outlier cluster already exist 
%       for a specified cluster. If so return the outlier cluster index;  
%       if not return zero.
%
%       This helper function is called by  pop_clustedit(), 
%       cls_moveoutlier(), & cls_renameclust(). 

function outlier_clust = cls_findoutlierclust(STUDY,clust)
outlier_clust = 0;
outlier_name = [ 'Outliers ' STUDY.cluster(clust).name ' ' ];
for k = 1:length(STUDY.cluster)
    if strncmpi(outlier_name,STUDY.cluster(k).name,length(outlier_name))
        outlier_clust = k;
    end
end
