function STUDY = std_mpcluster(STUDY,ALLEEG, numberOfClusters, outlierSTD, measuresToUseInClustering, methodParameter, conditionsToInclude )
% std_mpcluster() - Clusters STUDY ICs using Measure Product method. 
%                   In this method, IC measures, except equiv. dipoles, (ERP, ERSP...)  
%                   are compared for each IC pair and their dissimilarity is multiplied
%                   together to form a combined pairwise dissimilarity matrix. This matrix
%                   is then normalized, weighted and added to the normalized and weighted 
%                   IC equiv. dipole distance matrix. The final dissimilarity matrix is
%                   then clustered using affinity clustering  method. 
%                   You can control the effect of equiv. dipole distances in
%                   the clustering by setting the 'Relative dipole weight'
%                   parameter in the pop-uo GUI. For example, by setting
%                   this value to 0.8, the final dissimilarity matrix will consist of 80% 
%                   distance dissimilarity and 20% of other measures combined together.
%                   Please note that the number of returned clusters may slighty 
%                   differ from the number requested in the GUI.
%                   Notes:
%                    1) Currently only clustering the parent cluster (containing all components) is supported.
%                    2) You need to first calculate pairwise similarity matrices (for measures which
%                       you intend to used in the clustering) with std_mpreclust() or pop_mpreclust() before 
%                       running this function.
%                    3) The number of returned clusters may slightly (up to 5%) differ from the number requested.
% 
% Usage:   
%           >> STUDY = std_mpcluster(STUDY,ALLEEG, numberOfClusters, outlierSTD, .. 
%                                  measuresToUseInClustering, dipoleParameter, conditionsToInclude );
%              
%  Inputs:
%           STUDY                       - an EEGLAB STUDY set of loaded EEG structures
%           ALLEEG                      - vector of EEG datasets (can also be one EEG set).
%                                          must contain the dataset of interest.
%           numberOfClusters            - [integer] number of clusters to be created. 
%                                                   Note: The actual number of clusters that 
%                                                   are returned may be up to 5% different from this.
%           outlierSTD                  - [number] identify as outliers any components further than 
%                                                  the given number of standard deviations (of
%                                                  combined distance) from any cluster examplar.  
%                                                  Use Inf to prevent outlier identification.
%                                                  (otherwise a value between 2 and 5 is suggested)
%              
%           measuresToUseInClustering   - [cell array] measures to include in clustering.
%                                                 Any combination of ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%             
%           dipoleParameter             - (0 <= number <= 1) relative weighting of the equivalent
%                                         dipole locations.
%                                            
% Example:
%     
%     >> STUDY = std_mpcluster(STUDY,ALLEEG, 10, Inf, {'dipole' ,'erp' ,'ersp'}, 0.8);    
%        creates 10 clusters using component equivalent dipole locations, ERPs, and ERSPs.
%        The final dissimilarity matrix will consist of 80% distance dissimilarity and 20% of other 
%        measures combined together. It does not indentify outliers (outlierSTD -> Inf).
%              
%  Outputs:
%                 STUDY      - EEGLAB STUDY structure with new clusters added.
%                              All previous clusters except the parent cluster will have been removed.
%              
%  See also  apclusterK(),  std_mpreclust(), pop_mpcluster(),  pop_mpreclust()
%  Authors: Nima Bigdely-Shamlo and Brandon Burdge, SCCN/INC/UCSD, 2008

   
% Copyright (C) 2008 Nima Bigdely Shamlo, SCCN/INC/UCSD, nima(removethis)@sccn.ucsd.edu
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


% if nargin<4 || isempty(measuresToUseInClustering)
%     measuresToUseInClustering = {'erp','ersp', 'itc','spec','dipole'};
% end;

% always use hybrid clustering. effectively we here make raduis clustering disabled.
clusteringMethod = 'hybrid';

%  if no measure is specified in input arguments, use meaures specified in preclustering
if  nargin<5 && (isfield(STUDY, 'preclust') && isfield(STUDY.preclust, 'similarity') && isfield(STUDY.preclust.similarity, 'measuresToUseInClustering'))
    measuresToUseInClustering = STUDY.preclust.similarity.measuresToUseInClustering;
end;

% if precluster data does not exist, cluster the parent cluster (#1)
if isempty(STUDY.etc) || ~isfield(STUDY.etc, 'preclust') || isempty(STUDY.etc.preclust)
    STUDY.etc.preclust.clustlevel = 1;
end;

% if measure product matrices are not present, calculate them
% if ~isfield(STUDY, 'preclust') || ~isfield(STUDY.preclust, 'similarity') || ~isfield(STUDY.preclust.similarity, 'finalDistances')
%     STUDY = std_mpreclust(STUDY,ALLEEG, measuresToUseInClustering, true);
% else
%     finalDistances = STUDY.preclust.similarity.finalDistances;
% end;

numberOfComponents = length(STUDY.cluster(STUDY.etc.preclust.clustlevel).comps);
%numberOfComponents = size(finalDistances,1);

if nargin<4
    outlierSTD = Inf;
end;

if nargin<6
    clusteringMethod = 'hybrid';
    methodParameter = 0.5;
end;



fprintf('Clustering with %s Measure Product method into %d clusters. Measures to be used:\n', clusteringMethod, numberOfClusters);
for i=1:length(measuresToUseInClustering)
    fprintf('%d-%s\n',i,measuresToUseInClustering{i});
end;

%%


probabilityOfNotSameClass = ones(numberOfComponents, numberOfComponents);

if ismember('erp', measuresToUseInClustering)
    probabilityOfNotSameClass = probabilityOfNotSameClass .* (1-abs(STUDY.preclust.similarity.erpCorr));
end;

if ismember('ersp', measuresToUseInClustering)
    probabilityOfNotSameClass = probabilityOfNotSameClass .* (1-((1+STUDY.preclust.similarity.erspCorr)/2));
end;

if ismember('itc', measuresToUseInClustering)
    probabilityOfNotSameClass = probabilityOfNotSameClass .* (1-((1+STUDY.preclust.similarity.itcCorr)/2));
end;

if ismember('spec', measuresToUseInClustering)
    probabilityOfNotSameClass = probabilityOfNotSameClass .* (1-((1+STUDY.preclust.similarity.specCorr)/2));
end;

if ismember('map', measuresToUseInClustering)
    probabilityOfNotSameClass = probabilityOfNotSameClass .* (1-((1+STUDY.preclust.similarity.mapCorr)/2));
end;

% extra dissimilatity provided externally
if ismember('extra', measuresToUseInClustering) && isfield(STUDY.preclust.similarity, 'extraDissimilarity') && ~isempty(STUDY.preclust.similarity.extraDissimilarity)
    fprintf('Also using additional dissimilarity matrix. \n');
    probabilityOfNotSameClass = probabilityOfNotSameClass .* STUDY.preclust.similarity.extraDissimilarity;
end;


if ismember('dipole', measuresToUseInClustering)
    
    compDistance = STUDY.preclust.similarity.compDistance;

    if length(measuresToUseInClustering) == 1
        finalDistances = compDistance;
    else  % other measures other than dipole are also requested
        if strcmp(clusteringMethod, 'hybrid') % hybrid mehod (euclidean combination)
            dipoleWeight = methodParameter;
            finalDistances = ( ((1-dipoleWeight) * (probabilityOfNotSameClass ./ norm(probabilityOfNotSameClass,'fro') ) ) .^2 + (dipoleWeight * (compDistance ./ norm(compDistance,'fro')) ).^2 ) .^ 0.5;
        else % radius mehod
            icMaxDistance = maxDistanceToICallowed(compDistance, numberOfClusters);
            for i=1:size(compDistance,1)
                pairsWithLongerDistance = find(compDistance(i,:)>(icMaxDistance(i)* methodParameter)) ;

                probabilityOfNotSameClass(i, pairsWithLongerDistance) = 1;
                probabilityOfNotSameClass(pairsWithLongerDistance,i) = 1;

            end;
            finalDistances = compDistance .* probabilityOfNotSameClass;
        end;


    end;
else
    finalDistances = probabilityOfNotSameClass;
end
%%

[examplarIdx,netsim,dpsim,expref]=apclusterK(-finalDistances, numberOfClusters, 5);

examplar = unique(examplarIdx);

IDX = zeros(1,length(examplar));
for i=1:length(examplar)
    comps = find(examplarIdx == examplar(i));
    IDX(comps) = i;
end;
%% find outlier ICs

stdOfDistances = std(finalDistances(:));
for i=1:numberOfComponents
    minDistanceToAnExamplar(i) = min(finalDistances(i,examplar))/ stdOfDistances;
end;
outlierIDs = find(minDistanceToAnExamplar>outlierSTD);

if ~isempty(outlierIDs)  
    IDX(outlierIDs) = 0; % ID of zero is for outliers.
end;

%% create clusters in STUDY

STUDY.preclust.similarity.IDX = IDX;
STUDY.preclust.similarity.finalDistances = finalDistances;
STUDY  = std_createclust(STUDY,IDX,[],'Measure Product');


function icMaxDistance = maxDistanceToICallowed(compDistance, numberOfClusters)
fprintf('Clustering based on dipole locations...\n');
[examplarIdx]=apclusterK(-compDistance, numberOfClusters, 0.1);
examplar = unique(examplarIdx);
IDX = zeros(1,length(examplar));
for i=1:length(examplar)
    comps = find(examplarIdx == examplar(i));
    IDX(comps) = i;
end;

for comp = 1:length(IDX)
    ics = find(IDX == IDX(comp));
    icMaxDistance(comp) = max(max(compDistance(comp, ics)));
end;

