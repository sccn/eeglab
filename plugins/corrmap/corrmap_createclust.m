% corrmap_createclust() - creates a new cluster in the STUDY structure
%                         containing the ICs selected and stored using 
%                         corrmap() function.
%
% Usage:
%          >> STUDY = corrmap_createclust(STUDY,CORRMAP)
%
% Inputs:
%  STUDY - input STUDY structure
%  ALLEEG - input ALLEEG structure
%
% Outputs:
%  STUDY - the input STUDY set structure modified with the new cluster.
%
% See also:  pop_corrmap(), corrmap()
%
% Author: F. Campos Viola, MRC-IHR, 19/10/2007

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos Viola, MRC-IHR
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

% revised by F Campos-Viola - corrmap1.01 (30/01/2009)

function STUDY = corrmap_createclust(STUDY,CORRMAP)

clusters = [];

cls = length(STUDY.cluster)+1;% number of cluster
nc  = length(STUDY.cluster); % index of last cluster

%%%%% CREATING NEW CLUSTER %%%%%

STUDY.cluster(cls).name = CORRMAP.input.clname; %cluster's name

%auxiliary variables to create matrix with indices - used to access
%correct preclustdata

% saving info int the STUDY.cluster fields
STUDY.cluster(cls).sets  = CORRMAP.corr.sets{2}(1:CORRMAP.clust.ics(2))';
STUDY.cluster(cls).comps = CORRMAP.corr.ics{2}(1:CORRMAP.clust.ics(2))';
STUDY.cluster(cls).algorithm = 'correlation (CORRMAP)';
STUDY.cluster(cls).parent = 'ParentCluster 1';
STUDY.cluster(cls).child = [];
STUDY.cluster(cls).centroid = [];
%STUDY.cluster(cls).preclust.preclustdata = STUDY.etc.preclust.preclustdata(ind,:);
STUDY.cluster(cls).preclust.preclustdata = STUDY.etc.preclust.preclustdata
STUDY.cluster(cls).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
STUDY.cluster(cls).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;

%update parents clusters with cluster child indices
STUDY.cluster(1).child{cls-1} = CORRMAP.input.clname;

% clusters = [clusters 1:cls];%the new created clusters indices.

fprintf('>');
fprintf('a new cluster was created. \n'); %info for the user
