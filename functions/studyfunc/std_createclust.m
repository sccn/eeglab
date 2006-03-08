function [STUDY, clusters] = std_createclust(STUDY,IDX,C, algorithm)

clusters = [];
sets = [];
comp = [];
Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end
nsets = length(STUDY.etc.preclust.preclustcomps); 
for k = 1: nsets
    sets = [sets k*ones(1,length(STUDY.etc.preclust.preclustcomps{k}))];
    comp = [comp STUDY.etc.preclust.preclustcomps{k}];
end
if Ncond > 1
    tmp = ones(Ncond, length(sets));
    for l = 1:Ncond
        tmp(l,:) = sets + (l-1)*nsets;
    end
    sets = tmp;
    clear tmp
end

cls = size(C,1);
% Find the next available cluster index
nc = 0; 
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
    l = 0;
    nc = nc + 1;
    len = len +1;
else
    l = 1;
end

for k = l:cls 
    tmp = find(IDX==k);
    STUDY.cluster(k+len).sets = sets(:,tmp);
    STUDY.cluster(k+len).comps = comp(tmp);
    if k == 0
        STUDY.cluster(len).name = ['outlier ' num2str(k+nc)];
    else
        STUDY.cluster(k+len).name = [ 'Cls ' num2str(k+nc)];
    end
    STUDY.cluster(k+len).algorithm = algorithm;
    if STUDY.etc.preclust.clustlevel == 0
        STUDY.cluster(k+len).parent = [];
    else
        STUDY.cluster(k+len).parent{end+1} = STUDY.cluster(STUDY.etc.preclust.clustlevel).name;
        %update parents clusters with cluster child indices
        STUDY.cluster(STUDY.etc.preclust.clustlevel).child{end+1} = STUDY.cluster(k+nc).name;
    end
    STUDY.cluster(k+len).child = [];
    STUDY.cluster(k+len).preclust.preclustdata = STUDY.etc.preclust.preclustdata(tmp,:);
    STUDY.cluster(k+len).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
    STUDY.cluster(k+len).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;
end

clusters = [ clusters l+len:cls+len];%the new created clusters indices.
