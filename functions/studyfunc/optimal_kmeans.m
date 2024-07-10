function [IDX, C, sumd, D] = optimal_kmeans(clustdata, clusnum)
% THe aim of this function is to find the optimal number of clusters based
% on the "clusdata" and the "clusnum" range.
% INPUT: 
%       clustdata: the data matrix for the componenets. It is the same
%       varuable used in EEGLAB pop_clus.m
%       clusnum: number of cluster is a [1 x 2] vector specifying the lower
%       and upper bound of the optimal kmeans search.
%
% NOTE: This fucntion uses MALTAB's Statistics and Machine Learning
% Toolbox, and can't run without that.
%
% Created by Seyed Yahya Shirazi @ UCF 2019

%% intialize
if size(clusnum,2) ~= 2
    error('when using OPTIMAL_KMEANS, clusnum shoub be a 1x2 vector')
end

% try
% clust = parcluster;
% clust.NumWorkers = clust.NumWorkers - 2;
% parpool(clust);
% catch
% end

%% evaluate
disp('Finding optimal number of clusters. please be patient \n')
tic
eva = evalclusters(clustdata,@par_kmeans,'silhouette','KList',clusnum(1):clusnum(2));
toc
[IDX,C, sumd, D] = kmeans(clustdata,eva.OptimalK,'Replicates',10);  

end % end of the main funtion

function C = par_kmeans(DATA,K)
%     if gpuDeviceCount
%         disp('using GPU for kmeans')
%         X = gpuArray(DATA);
%         gC = kmeans(X,K,'Replicates',10);
%         C = gather(gC);
%     else
%         disp('using parpool for kmeans')
        C = kmeans(DATA,K,'Replicates',5); % ,'Options',statset('UseParallel',1)     
%     end
end % end of function
