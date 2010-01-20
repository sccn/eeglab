function [STUDY, clusters] = std_createclusters(STUDY,IDX, algorithm, parentClusterNumber)
% 

if nargin<4
    parentClusterNumber = 1;
    STUDY.cluster(2:end) = [];
    STUDY.cluster(1).child = [];
end;
    



% IDX - index of cluster for each component. Ex: 63 components and 2
% clusters: IDX will be a 61x1 vector of 1 and 2 (and 0=outlisers)
% C - centroid for clusters. If 2 clusters, size will be 2 x
%     width of the preclustering matrix


clusters = [];
sets = [];
comp = [];
nsets = length(STUDY.etc.preclust.preclustcomps); % number of sets from which components are originated, maybe less than actual number of sets in the STUDY since some sets may have not any selected componnets
for k = 1: nsets
    sets = [sets k*ones(1,length(STUDY.etc.preclust.preclustcomps{k}))];
    comp = [comp STUDY.etc.preclust.preclustcomps{k}];
end

% Find the next available cluster index
% -------------------------------------
 cls = length(unique(IDX)); % number of cluster = number of row of centroid matrix
nc  = 0; % index of last cluster 
for k =  1:length(STUDY.cluster)
    ti = strfind(STUDY.cluster(k).name, ' ');
    tmp = STUDY.cluster(k).name(ti(end) + 1:end);
    nc = max(nc,str2num(tmp));
    % check if there is a cluster of Notclust components
    if strcmp(STUDY.cluster(k).parent,STUDY.cluster(STUDY.etc.preclust.clustlevel).name) 
        STUDY.cluster(k).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
        clusters = [clusters k];
    end
end
len = length(STUDY.cluster);

if ~isempty(find(IDX==0)) %outliers exist
    firstind = 0;
    nc  = nc + 1;
    len = len + 1;
else
    firstind = 1;
end

% create all clusters
% -------------------
for k = firstind:cls 
    
    % cluster name
    % ------------
    if k == 0
         STUDY.cluster(len).name   = [ 'outlier ' num2str(k+nc)];
    else STUDY.cluster(k+len).name = [ 'Cls ' num2str(k+nc)];
    end

    % find indices
    % ------------
    tmp = find(IDX==k);
    STUDY.cluster(k+len).sets  = STUDY.cluster(STUDY.etc.preclust.clustlevel).sets(:,tmp);
    STUDY.cluster(k+len).comps = STUDY.cluster(STUDY.etc.preclust.clustlevel).comps(tmp);
    STUDY.cluster(k+len).algorithm = algorithm;
    STUDY.cluster(k+len).parent{end+1} = STUDY.cluster(STUDY.etc.preclust.clustlevel).name;
    STUDY.cluster(k+len).child = [];
   % STUDY.cluster(k+len).preclust.preclustdata = STUDY.etc.preclust.preclustdata(tmp,:);
   % STUDY.cluster(k+len).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
    STUDY.cluster(k+len).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;

    %update parents clusters with cluster child indices
    % -------------------------------------------------
    STUDY.cluster(STUDY.etc.preclust.clustlevel).child{end+1} = STUDY.cluster(k+nc).name;
end

clusters = [ clusters firstind+len:cls+len];%the new created clusters indices.

















% clusters = [];
% sets = [];
% comp = [];
% nsets = length(STUDY.etc.preclust.preclustcomps); % number of sets from which components are originated, maybe less than actual number of sets in the STUDY since some sets may have not any selected componnets
% for k = 1: nsets
%     sets = [sets k*ones(1,length(STUDY.etc.preclust.preclustcomps{k}))];
%     comp = [comp STUDY.etc.preclust.preclustcomps{k}];
% end
% 
% % Find the next available cluster index
% % -------------------------------------
% cls = length(unique(IDX)); % number of cluster = number of row of centroid matrix
% nc  = 0; % index of last cluster 
% for k =  1:length(STUDY.cluster)
%     ti = strfind(STUDY.cluster(k).name, ' ');
%     tmp = STUDY.cluster(k).name(ti(end) + 1:end);
%     nc = max(nc,str2num(tmp));
%     % check if there is a cluster of Notclust components
%     if strcmp(STUDY.cluster(k).parent,STUDY.cluster(STUDY.etc.preclust.clustlevel).name) 
%         STUDY.cluster(k).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
%         clusters = [clusters k];
%     end
% end
% len = length(STUDY.cluster);
% 
% if ~isempty(find(IDX==0)) %outliers exist
%     firstind = 0;
%     nc  = nc + 1;
%     len = len + 1;
% else
%     firstind = 1;
% end
% 
% % create all clusters
% % -------------------
% for k = firstind:cls 
%     
%     % cluster name
%     % ------------
%     if k == 0
%          STUDY.cluster(len).name   = [ 'outlier ' num2str(k+nc)];
%     else STUDY.cluster(k+len).name = [ 'Cls ' num2str(k+nc)];
%     end
% 
%     % find indices
%     % ------------
%     tmp = find(IDX==k);
%     STUDY.cluster(k+len).sets  = STUDY.cluster(STUDY.etc.preclust.clustlevel).sets(:,tmp);
%     STUDY.cluster(k+len).comps = STUDY.cluster(STUDY.etc.preclust.clustlevel).comps(tmp);
%     STUDY.cluster(k+len).algorithm = algorithm;
%     STUDY.cluster(k+len).parent{end+1} = STUDY.cluster(STUDY.etc.preclust.clustlevel).name;
%     STUDY.cluster(k+len).child = [];
%     
%     if strcmp(algorithm,'Measure Product')
%         STUDY.cluster(k+len).preclust.preclustdata = []; % just for now
%         STUDY.cluster(k+len).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;  
%     else
%         STUDY.cluster(k+len).preclust.preclustdata = STUDY.etc.preclust.preclustdata(tmp,:);
%         STUDY.cluster(k+len).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;
%     end;
% 
%     STUDY.cluster(k+len).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
% 
%     %update parents clusters with cluster child indices
%     % -------------------------------------------------
%     STUDY.cluster(STUDY.etc.preclust.clustlevel).child{end+1} = STUDY.cluster(k+nc).name;
% end
% 
% clusters = [ clusters firstind+len:cls+len];%the new created clusters indices.
