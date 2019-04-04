% robust_kmeans() - an extension of Matlab kmeans() that removes outlier 
%        components from all clusters. 
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

function  [IDX,C,sumd,D,outliers] = robust_kmeans(data,N,STD,MAXiter,method)
% data - pre-clustering data matrix.
% N - number of wanted clusters.

if nargin < 5
    method = 'kmeans';
end

flag  = 1;
not_outliers = 1:size(data,1);
old_outliers = [];
if strcmpi(method, 'kmeans')
    [IDX,C,sumd,D] = kmeans(data,N,'replicates',30,'emptyaction','drop'); % Cluster using K-means algorithm
else
    [IDX,C,sumd,D] = kmeanscluster(data,N); % Cluster using K-means algorithm
end;    
if STD >= 2 % STD for returned outlier
    rSTD = STD -1;
else
    rSTD = STD;
end
loop = 0;

while flag
     loop =  loop + 1;
	std_all = [];
    ref_D = 0;
	for k = 1:N
        tmp = ['cls' num2str(k) ' = find(IDX=='  num2str(k) ')''; ' ]; %find the component indices belonging to each cluster (cls1 = ...).
        eval(tmp);
        tmp = ['std' num2str(k) ' = std(D(cls'  num2str(k) ' ,' num2str(k) ')); ' ]; %compute the std of each cluster
        eval(tmp);
        std_all = [std_all ['std' num2str(k)  '  ']];
        tmp = [ 'ref_D = ' num2str(ref_D) ' + mean(D(cls'  num2str(k) ' ,' num2str(k) '));' ];
        eval(tmp);
	end
	std_all = [ '[ ' std_all ' ]' ];
    std_all = eval(std_all);
    
	% Find the outliers
    % Outlier definition - its distance from its cluster center is bigger
    % than STD times the std of the cluster, as long as the distance is bigger
    % than the mean distance times STD (avoid problems where all points turn to be outliers).
	outliers = [];
    ref_D = ref_D/N;
	for k = 1:N
        tmp = ['cls' num2str(k) '(find(D(find(IDX=='  num2str(k) ')'' , ' num2str(k) ') > ' num2str(STD)  '*std' num2str(k) ')); ' ];
        optionalO = eval(tmp);
        Oind = find(D(optionalO,k) >  ref_D*STD);
        outliers = [outliers optionalO(Oind)];
	end
    if isempty(outliers) || (loop == MAXiter)
        flag = 0;
    end
    l = length(old_outliers);
    returned_outliers = [];
 
    
    for k = 1:l
        tmp = sum((C-ones(N,1)*data(old_outliers(k),:)).^2,2)'; % Find the distance of each former outlier to the current cluster
        if isempty(find(tmp <= std_all*rSTD))  %Check if the outlier is still an outlier (far from each cluster center more than STD-1 times its std).
            returned_outliers = [returned_outliers old_outliers(k)];
        end
    end
    outliers = not_outliers(outliers);
    outliers = [outliers returned_outliers ];
	tmp = ones(1,size(data,1));
	tmp(outliers) = 0;
	not_outliers = (find(tmp==1));
    
    if strcmpi(method, 'kmeans')
        [IDX,C,sumd,D] = kmeans(data(not_outliers,:),N,'replicates',30,'emptyaction','drop');
    else
        [IDX,C,sumd,D] = kmeanscluster(data(not_outliers,:),N);
    end
    old_outliers = outliers;
    old_IDX = zeros(size(data,1),1);
    old_IDX(sort(not_outliers)) = IDX;
    
end

IDX = old_IDX;
