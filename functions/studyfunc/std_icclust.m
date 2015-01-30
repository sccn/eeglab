% pop_icclust() - Mixed Affinity Propagation Clustering method for eeglab STUDY
%                 Call the function limo_clusterica to compute the clustering 
%                 based on the inputs defined in this function
%                  
%
% Usage:
%   >>  [STUDY ALLEEG IDX] = std_icclust(STUDY, ALLEEG,1, 1, {'scalp' 'abso' 1},...
%                                                            {'spec' 'freqrange' [3 25] },...
%                                                            {'erp' 'timewindow' []});
% 
%   >>  [STUDY ALLEEG IDX] = std_icclust(STUDY, ALLEEG,1, 1, {'scalp' 'abso' 1},...
%                                                            {'spec' 'freqrange' [3 25] },...
%                                                            {'erp' 'timewindow' []},...
%                                                            {'dipoles'});
% 
%
% Inputs:
%   clustind     - a cluster index for further (hierarchical) clustering -
%                  for example to cluster a spectrum-based mu-rhythm cluster into 
%                  dipole location-based left mu and right mu sub-clusters. 
%                  Should be empty for first stage (whole-STUDY) clustering {default: []}
% useclust       - Index of the results from limo_clusterica.m to use to
%                  perform the clustering
%
%   preproc      - {'command' 'key1' val1 'key2' val2 ...} component clustering measures to prepare
%
%            * 'command' = component measure to compute:
%                    'erp'     = cluster on the component ERPs,
%                    'dipoles' = cluster on the component [X Y Z] dipole locations
%                    'spec'    = cluster on the component log activity spectra (in dB)
%                                  (with the baseline mean dB spectrum subtracted).
%                    'scalp'   = cluster on component (topoplot()) scalp maps 
%                                  (or on their absolute values),
%                    'scalpLaplac' = cluster on component (topoplot()) laplacian scalp maps
%                                  (or on their absolute values),
%                    'scalpGrad' = cluster on the (topoplot()) scalp map gradients 
%                                  (or on their absolute values),
%                    'ersp'    = cluster on components ERSP. (requires: 'cycles', 
%                                  'freqrange', 'padratio', 'timewindow', 'alpha').
%                    'itc'     = cluster on components ITC.(requires: 'cycles', 
%                                  'freqrange', 'padratio', 'timewindow', 'alpha').
% Optional inputs:
%                    'freqrange'  = [min max] frequency range (in Hz) to select
%                                   spectrum, 'ersp', and 'itc' measures.  
%                    'timewindow' = [min max] time window (in sec) to include in 'erp',
%                                   'ersp', and 'itc' measures.  
%                    'abso'    =  [0|1] 1 = use absolute values of topoplot(), gradient, or 
%                                   Laplacian maps {default: 1}
% 
% Outputs:
%
% See also:
%
%
% Author: % Cyril Pernet          , The University of Edinburgh
%           Arnaud Delorme        , SCCN, INC, UCSD
%           Ramon Martinez-Cancino, SCCN, INC, UCSD
%
% Copyright (C) 2014  Ramon Martinez-Cancino, SCCN, INC, UCSD
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

function [STUDY, ALLEEG, IDX, com] = std_icclust(STUDY, ALLEEG, cluster_ind,useclust, varargin)
% function [ STUDY, ALLEEG ] = std_preclust(STUDY, ALLEEG, cluster_ind, varargin)

if nargin < 2
    help std_icclust;
    return;
end;

if nargin == 2
    cluster_ind = 1; % default to clustering the whole STUDY
    useclust    = 3; % Default output selection from limo_clusterica
end

% check dataset length consistency before computing
% -------------------------------------------------
pnts  = ALLEEG(STUDY.datasetinfo(1).index).pnts;
srate = ALLEEG(STUDY.datasetinfo(1).index).srate;
xmin  = ALLEEG(STUDY.datasetinfo(1).index).xmin;
xmax  = ALLEEG(STUDY.datasetinfo(1).index).xmax;
for index = 1:length(STUDY.datasetinfo)
    ind = STUDY.datasetinfo(index).index;
    if srate ~= ALLEEG(ind).srate, error(sprintf('Dataset %d does not have the same sampling rate as dataset 1', ind)); end;
    if ~all([ ALLEEG.trials ] == 1)
        if abs(xmin-ALLEEG(ind).xmin) > 1e-7, warning(sprintf('Dataset %d does not have the same time limit as dataset 1', ind)); end;
        if abs(xmax-ALLEEG(ind).xmax) > 1e-7, warning(sprintf('Dataset %d does not have the same time limit as dataset 1', ind)); end;
        if pnts ~= ALLEEG(ind).pnts, error(sprintf('Dataset %d does not have the same number of point as dataset 1', ind)); end;
    end;
end;

% Get component indices that are part of the cluster
% --------------------------------------------------
if isempty(cluster_ind)
    cluster_ind = 1;
end;
if length(cluster_ind) ~= 1
    error('Only one cluster can be sub-clustered. To sub-cluster multiple clusters, first merge them.');
end

if isempty(useclust)
    useclust = 3;
end;
if length(useclust) ~= 1  || useclust > 3
    error('Checky third input. useclust must be 1,2 or 3');
end

% Decode input arguments
% ----------------------
update_flag  = 0;
rv           = 1;

indrm        = [];
for index = 1:length(varargin) % scan commands
    strcom = varargin{index}{1};
    if strcmpi(strcom,'dipselect')
        update_flag = 1;
        rv = varargin{index}{3};
        indrm = [indrm index]; %remove this command
    elseif strcmpi(strcom,'finaldim') % second level pca
        secondlevpca = varargin{index}{3};
        indrm = [indrm index]; %remove this command
    end
end
varargin(indrm) = []; %remove commands

% Scan which component to remove (no dipole info)
% -----------------------------------------------
if update_flag % dipole information is used to select components
    error('Update flag is obsolete');
end;

% scan all commands
% -----------------
clustdata = [];
erspquery = 0;
lmclust_vararg = '';

for index = 1:length(varargin)
    
    % decode inputs
    % -------------
    strcom = varargin{index}{1};
    if any(strcom == 'X'), disp('character ''X'' found in command'); end;
    %defult values
    freqrange = [];
    timewindow = [];
    abso = 1;
    fun_arg = [];
    
    for subind = 2:2:length(varargin{index})
        switch varargin{index}{subind}
            case 'norm'
                norm = varargin{index}{subind+1};
            case 'freqrange'
                freqrange = varargin{index}{subind+1};
            case 'timewindow'
                timewindow = varargin{index}{subind+1};
            case 'abso'
                abso = varargin{index}{subind+1};
            case 'savetrials'
                error('You may now use the function std_precomp to precompute measures');
            case 'cycles'
                error('You may now use the function std_precomp to precompute measures');
            case 'alpha'
                error('You may now use the function std_precomp to precompute measures');
            case 'padratio'
                error('You may now use the function std_precomp to precompute measures');
            otherwise
                fun_arg{length(fun_arg)+1} =  varargin{index}{subind+1};
        end
    end
    
    % scan datasets
    % -------------
    if strcmpi(strcom, 'scalp'),           scalpmodif = 'none';
    elseif strcmpi(strcom, 'scalpLaplac'), scalpmodif = 'laplacian';
    else                                   scalpmodif = 'gradient';
    end;
    
    % check that all datasets are in preclustering for current design
    % ---------------------------------------------------------------
    tmpstruct = std_setcomps2cell(STUDY, STUDY.cluster(cluster_ind).sets, STUDY.cluster(cluster_ind).comps);
    alldatasets = unique_bc(STUDY.cluster(cluster_ind).sets(:));
    
    if length(alldatasets) < length(STUDY.datasetinfo) && cluster_ind == 1
        error( [ 'Some datasets not included in preclustering' 10 ...
            'because of partial STUDY design. You need to' 10 ...
            'use a STUDY design that includes all datasets.' ]);
    end;
    
    for si = 1:size(STUDY.cluster(cluster_ind).sets,2)
        
        % test consistency of the .set structure
        % --------------------------------------
        if strcmpi(strcom, 'erp') || strcmpi(strcom, 'spec') || strcmpi(strcom, 'ersp') || strcmpi(strcom, 'itc')
            if any(isnan(STUDY.cluster(cluster_ind).sets(:)))
                error( [ 'std_preclust error: some datasets do not have ICA pairs.' 10 ...
                    'Look for NaN values in STUDY.cluster(1).sets which' 10 ...
                    'indicate missing datasets. FOR CLUSTERING, YOU MAY ONLY' 10 ...
                    'USE DIPOLE OR SCALP MAP CLUSTERING.' ]);
            end;
        end;
        
        switch strcom
            
            % select ica component ERPs
            % -------------------------
            case 'erp',
                % read and concatenate all cells for this specific set
                % of identical ICA decompositions
                STUDY.cluster = checkcentroidfield(STUDY.cluster, 'erp', 'erp_times');
                tmpstruct = std_setcomps2cell(STUDY, STUDY.cluster(cluster_ind).sets(:,si), STUDY.cluster(cluster_ind).comps(si));
                cellinds  = [ tmpstruct.setinds{:} ];
                compinds  = [ tmpstruct.allinds{:} ];
                cells = STUDY.design(STUDY.currentdesign).cell(cellinds);
                fprintf('Creating input array for Affinity Propagation Clustering array row %d, adding ERP for design %d cell(s) [%s] component %d ...\n', si, STUDY.currentdesign, int2str(cellinds), compinds(1));
                X = std_readfile( cells, 'components', compinds, 'timelimits', timewindow, 'measure', 'erp');
                X = abs(X(:)'); % take the absolute value of the ERP to avoid polarities issues
                
                % select ica scalp maps
                % --------------------------
            case { 'scalp' 'scalpLaplac' 'scalpGrad' }
                idat  = STUDY.datasetinfo(STUDY.cluster(cluster_ind).sets(:,si)).index;
                icomp = STUDY.cluster(cluster_ind).comps(si);
                fprintf('Creating input array for Affinity Propagation Clustering row %d, adding interpolated scalp maps for dataset %d component %d...\n', si, idat, icomp);
                X = std_readtopo(ALLEEG, idat, icomp, scalpmodif, 'preclust');
                
                % select ica comp spectra
                % -----------------------
            case 'spec',
                % read and concatenate all cells for this specific set
                % of identical ICA decompositions
                STUDY.cluster = checkcentroidfield(STUDY.cluster, 'spec', 'spec_freqs');
                tmpstruct = std_setcomps2cell(STUDY, STUDY.cluster(cluster_ind).sets(:,si), STUDY.cluster(cluster_ind).comps(si));
                cellinds  = [ tmpstruct.setinds{:} ];
                compinds  = [ tmpstruct.allinds{:} ];
                cells = STUDY.design(STUDY.currentdesign).cell(cellinds);
                fprintf('Creating input array for Affinity Propagation Clustering row %d, adding spectrum for design %d cell(s) [%s] component %d ...\n', si, STUDY.currentdesign, int2str(cellinds), compinds(1));
                X = std_readfile( cells, 'components', compinds, 'freqlimits', freqrange, 'measure', 'spec');
                if size(X,2) > 1, X = X - repmat(mean(X,2), [1 size(X,2)]); end;
                X = X - repmat(mean(X,1), [size(X,1) 1]);
                X = X(:)';
                
                % select dipole information
                % -------------------------
            case 'dipoles'
                idat  = STUDY.datasetinfo(STUDY.cluster(cluster_ind).sets(1,si)).index;
                icomp = STUDY.cluster(cluster_ind).comps(si);
                fprintf('Creating input array for Affinity Propagation Clustering row %d, adding dipole for dataset %d component %d...\n', si, idat, icomp);
                try
                    % select among 3 sub-options
                    % --------------------------
                    ldip = 1;
                    if size(ALLEEG(idat).dipfit.model(icomp).posxyz,1) == 2 % two dipoles model
                        if any(ALLEEG(idat).dipfit.model(icomp).posxyz(1,:)) ...
                                && any(ALLEEG(idat).dipfit.model(icomp).posxyz(2,:)) %both dipoles exist
                            % find the leftmost dipole
                            [garb ldip] = max(ALLEEG(idat).dipfit.model(icomp).posxyz(:,2));
                        elseif any(ALLEEG(idat).dipfit.model(icomp).posxyz(2,:))
                            ldip = 2; % the leftmost dipole is the only one that exists
                        end
                    end
                    X = ALLEEG(idat).dipfit.model(icomp).posxyz(ldip,:);
                catch
                    error([ sprintf('Some dipole information is missing (e.g. component %d of dataset %d)', icomp, idat) 10 ...
                        'Components are not assigned a dipole if residual variance is too high so' 10 ...
                        'in the STUDY info editor, remember to select component by residual' 10 ...
                        'variance (column "select by r.v.") prior to preclustering them.' ]);
                end
                
                % cluster on ica ersp / itc values
                % --------------------------------
            case  {'ersp', 'itc' }
                % read and concatenate all cells for this specific set
                % of identical ICA decompositions
                STUDY.cluster = checkcentroidfield(STUDY.cluster, 'ersp', 'ersp_times', 'ersp_freqs', 'itc', 'itc_times', 'itc_freqs');
                tmpstruct = std_setcomps2cell(STUDY, STUDY.cluster(cluster_ind).sets(:,si), STUDY.cluster(cluster_ind).comps(si));
                cellinds  = [ tmpstruct.setinds{:} ];
                compinds  = [ tmpstruct.allinds{:} ];
                cells = STUDY.design(STUDY.currentdesign).cell(cellinds);
                fprintf('Pre-clustering array row %d, adding %s for design %d cell(s) [%s] component %d ...\n', si, upper(strcom), STUDY.currentdesign, int2str(cellinds), compinds(1));
                X = std_readfile( cells, 'components', compinds, 'timelimits', timewindow, 'measure', strcom);
        end;
        
        % copy data in the array
        % ----------------------
        if ~isreal(X) X = abs(X); end; % for ITC data
        X = reshape(X, 1, numel(X));
        if si == 1, data = zeros(size(STUDY.cluster(cluster_ind).sets,2),length(X)); end;
        data(si,:) = X;
        try
            data(si,:) = X;
        catch,
            error([ 'This type of  clustering requires that all subjects' 10 ...
                'be represented for all combination of selected independent' 10 ...
                'variables in the current STUDY design. In addition, for each' 10 ...
                'different ICA decomposition included in the STUDY (some' 10 ...
                'datasets may have the same decomposition), at least one' 10 ...
                'dataset must be represented.' ]);
        end;
    end; % end scan datasets
    
    % Storing data
    switch strcom
        case 'erp'
            apdata.erp = data;
            lmclust_vararg = [lmclust_vararg ',''erp'', apdata.erp'];
        case 'ersp'
            apdata.ersp = data;
            lmclust_vararg = [lmclust_vararg ',''ersp'', apdata.ersp'];
        case 'dipoles'
            apdata.dipoles = data;
            lmclust_vararg = [lmclust_vararg ',''dipoles'', apdata.dipoles'];
        case'spec'
            apdata.spec = data;
            lmclust_vararg = [lmclust_vararg ',''spect'', apdata.spec'];
        case'itc'
            apdata.itc = data;
            lmclust_vararg = [lmclust_vararg ',''itc'', apdata.itc'];
        case  { 'scalp' 'scalpLaplac' 'scalpGrad' }
            apdata.smap = data;
            lmclust_vararg = [lmclust_vararg ',''smap'', apdata.smap'];       
    end 
                
end

% Affinity Propagation Clustering
display('Computing clusters using Mixed Affinity Propagation Clustering method ...')
com = ['IDX = limo_clusterica(' lmclust_vararg(2:end) ')'];
eval(com);
if isempty(IDX{useclust})
    warning(['std_icclust warning: Unexpected output form clustering........ aborting clustering']);
    return
end

% Storing data
if ~isempty(STUDY.etc.apclust.data)
    display('Replacing clustering data with the new one computed ...');
    STUDY.etc.apclust = [];
end
    
STUDY.etc.apclust.data   = apdata;
STUDY.etc.apclust.preclustparams = varargin;
STUDY.etc.apclust.IDX = IDX{useclust};

if ~isempty(cluster_ind)
    STUDY.etc.preclust.clustlevel = cluster_ind;
else
    STUDY.etc.preclust.clustlevel = 1; % No parent cluster (cluster on all components in STUDY).
end

% Creating clusters in STUDY
[STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','multiapclust','clus_num',   length(unique(IDX{useclust})));

return

% the function below checks the precense of the centroid field
function cluster = checkcentroidfield(cluster, varargin)
for kk = 1:length(cluster)
    if ~isfield('centroid', cluster(kk)), cluster(kk).centroid = []; end;
    for vi = 1:length(varargin)
        if isfield(cluster(kk).centroid, varargin{vi})
            cluster(kk).centroid = rmfield(cluster(kk).centroid, varargin{vi});
        end;
    end;
end;
