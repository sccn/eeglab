function outlier_clust = cls_findoutlierclust(STUDY,clust)
% Finds if Outlier cluster already exist for a specific cluster
% if yes returns the outlier cluster index, if not returns zero
outlier_clust = 0;
outlier_name = [ 'Outliers ' STUDY.cluster(clust).name ' ' ];
for k = 1:length(STUDY.cluster)
    if strncmpi(outlier_name,STUDY.cluster(k).name,length(outlier_name))
        outlier_clust = k;
    end
end
