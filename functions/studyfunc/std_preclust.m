% std_preclust() - select measures to be included in computation of a preclustering array.
%                    This array is used by pop_clust() to find component clusters from a
%                    specified parent cluster.
%                    Selected measures (dipole location, ERPs, spectra, scalp maps, ERSPs,
%                    and/or ITCs) should already be precomputed using pop-precomp(). Each 
%                    feature dimension is reduced by PCA decomposition. These PCA matrices 
%                    (one per measure) are concatenated and used as input to the clustering
%                    algorithm in pop_clust(). Follow with pop_clust(). 
%                    See Example below:
%
%  >> [STUDY,ALLEEG] = std_preclust(STUDY,ALLEEG); % prepare to cluster all comps 
%                                                                % in all sets on all measures
%
%  >> [STUDY,ALLEEG] = std_preclust(STUDY,ALLEEG, clustind, preproc1, preproc2...);
%                                                                % prepare to cluster specifed 
%                                                                % cluster on specified measures
% Required inputs:
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%
% Optional inputs:
%   clustind     - a cluster index for further (hierarchical) clustering -
%                  for example to cluster a spectrum-based mu-rhythm cluster into 
%                  dipole location-based left mu and right mu sub-clusters. 
%                  Should be empty for first stage (whole-STUDY) clustering {default: []}
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
%                    'finaldim' = final number of dimensions. Enables second-level PCA. 
%                                  By default this command is not used (see Example below).
%
%            * 'key'   optional keywords and [valuess] used to compute the above 'commands':
%                    'npca'    =  [integer] number of principal components (PCA dimension) of 
%                                   the selected measures to retain for clustering. {default: 5}
%                    'norm'    =  [0|1] 1 -> normalize the PCA components so the variance of 
%                                   first principal component is 1 (useful when using several 
%                                   clustering measures - 'ersp','scalp',...). {default: 1}
%                    'weight'  =  [integer] weight with respect to other clustering measures.
%                    'freqrange'  = [min max] frequency range (in Hz) to include in activity 
%                                   spectrum, 'ersp', and 'itc' measures.  
%                    'timewindow' = [min max] time window (in sec) to include in 'erp',
%                                   'ersp', and 'itc' measures.  
%                    'abso'    =  [0|1] 1 = use absolute values of topoplot(), gradient, or 
%                                   Laplacian maps {default: 1}
%                    'erpfilter' = [double] Perform low pass filter on ERPs
%                                  at the frequency provided. This is done ONLY for
%                                  preclustering purposes and do not act on the data. Input unit: Hz
%                    'funarg'  =  [cell array] optional function arguments for mean spectrum 
%                                   calculation (>> help spectopo) {default: none}
% Outputs:
%   STUDY        - the input STUDY set with pre-clustering data added, for use by pop_clust() 
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures, modified by adding preprocessing 
%                  data as pointers to Matlab files that hold the pre-clustering component measures.
%
% Example:
%   >> [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
%                        { 'spec'  'npca' 10 'norm' 1 'weight' 1 'freqrange'  [ 3 25 ] } , ...
%                        { 'erp'   'npca' 10 'norm' 1 'weight' 2 'timewindow' [ 350 500 ] } ,...
%                        { 'scalp' 'npca' 10 'norm' 1 'weight' 2 'abso' 1 } , ...
%                        { 'dipoles'         'norm' 1 'weight' 15 } , ...
%                        { 'ersp'  'npca' 10 'freqrange' [ 3 25 ] 'cycles' [ 3 0.5 ] 'alpha' 0.01 ....
%                                  'padratio' 4 'timewindow' [ -1600 1495 ] 'norm' 1 'weight' 1 } ,...
%                        { 'itc'   'npca' 10 'freqrange' [ 3 25 ] 'cycles' [ 3 0.5 ] 'alpha' 0.01 ...
%                                  'padratio' 4 'timewindow' [ -1600 1495 ] 'norm' 1 'weight' 1 }, ...
%                        { 'finaldim' 'npca' 10 });
%                          
%                   % This prepares, for initial clustering, all components in the STUDY datasets
%                   % except components with dipole model residual variance (see function 
%                   % std_editset() for how to select such components).
%                   % Clustering will be based on the components' mean spectra in the [3 25] Hz 
%                   % frequency range, on the components' ERPs in the [350 500] ms time window, 
%                   % on the (absolute-value) component scalp maps, on the equivalent dipole 
%                   % locations, and on the component mean ERSP and ITC images. 
%                   % The final keyword specifies final PCA dimension reduction to 10
%                   % principal dimensions. See the clustering tutorial for more details.
%
% Authors: Arnaud Delorme, Hilit Serby & Scott Makeig, SCCN, INC, UCSD, May 13, 2004

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, May 13,2004, arno@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [ STUDY, ALLEEG ] = std_preclust(STUDY, ALLEEG, cluster_ind, varargin)
    
    if nargin < 2
        help std_preclust;
        return;
    end
    
    if nargin == 2
        cluster_ind = 1; % default to clustering the whole STUDY 
    end    

    % check dataset length consistency before computing
    % -------------------------------------------------
    pnts  = ALLEEG(STUDY.datasetinfo(1).index).pnts;
    srate = ALLEEG(STUDY.datasetinfo(1).index).srate;
    xmin  = ALLEEG(STUDY.datasetinfo(1).index).xmin;
    xmax  = ALLEEG(STUDY.datasetinfo(1).index).xmax;
    for index = 1:length(STUDY.datasetinfo)
        ind = STUDY.datasetinfo(index).index;
        if srate ~= ALLEEG(ind).srate, error(sprintf('Dataset %d does not have the same sampling rate as dataset 1', ind)); end
        if ~all([ ALLEEG.trials ] == 1)
            if abs(xmin-ALLEEG(ind).xmin) > 1e-7, warning(sprintf('Dataset %d does not have the same time limit as dataset 1', ind)); end
            if abs(xmax-ALLEEG(ind).xmax) > 1e-7, warning(sprintf('Dataset %d does not have the same time limit as dataset 1', ind)); end
            if pnts ~= ALLEEG(ind).pnts, error(sprintf('Dataset %d does not have the same number of point as dataset 1', ind)); end
        end
    end
    
    % Get component indices that are part of the cluster 
    % --------------------------------------------------
    if isempty(cluster_ind)
        cluster_ind = 1;
    end
    if length(cluster_ind) ~= 1
        error('Only one cluster can be sub-clustered. To sub-cluster multiple clusters, first merge them.');
    end
    
    % Decode input arguments
    % ----------------------
    update_flag  = 0;
    rv           = 1;
    secondlevpca = Inf;
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
    end
    
    % scan all commands
    % -----------------
    clustdata = [];
    erspquery = 0;
    for index = 1:length(varargin)
        
        % decode inputs
        % -------------
        strcom = varargin{index}{1};
        if any(strcom == 'X'), disp('character ''X'' found in command'); end
        %defult values
        npca = NaN;
        norm = 1;
        weight = 1;
        freqrange = [];
        timewindow = [];
        erpfilter  = [];
        abso = 1;
        fun_arg = [];
        savetrials = 'off';
        recompute  = 'on';
        clear data;
        for subind = 2:2:length(varargin{index})
            switch varargin{index}{subind}
                case 'npca'
                    npca = varargin{index}{subind+1};
                case 'norm'
                    norm = varargin{index}{subind+1};
                case 'weight'
                    weight = varargin{index}{subind+1};
                case 'freqrange' 
                    freqrange = varargin{index}{subind+1};
                case 'timewindow' 
                    timewindow = varargin{index}{subind+1};
                case 'abso'
                    abso = varargin{index}{subind+1};
                case 'erpfilter'
                    erpfilter = varargin{index}{subind+1};
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
        end
        
        % check that all datasets are in preclustering for current design
        % ---------------------------------------------------------------
        alldatasets = unique_bc(STUDY.cluster(cluster_ind).sets(:));
        
        if length(alldatasets) < length(STUDY.datasetinfo) && cluster_ind == 1
            error( [ 'Some datasets not included in preclustering' 10 ... 
                     'because of partial STUDY design. You need to' 10 ...
                     'use a STUDY design that includes all datasets.' ]);
        end
        
        switch strcom
            
            % select ica component ERPs
            % -------------------------
            case 'erp',
                [STUDY data, datatime] = std_readerp( STUDY, ALLEEG, 'design', NaN, 'clusters', cluster_ind, 'timerange', timewindow, 'datatype', 'erp', 'componentpol', 'off');
                % Filtering data to be plotted
                if ~isempty(erpfilter), data = {myfilt(data{:}, 1000/(datatime(2)-datatime(1)), 0, erpfilter)}; end;       
                data = data{1}';
                
            % select ica component spectrum
            % -----------------------------
            case 'spec',
                [STUDY data] = std_readerp( STUDY, ALLEEG, 'design', NaN, 'clusters', cluster_ind, 'freqrange', freqrange, 'datatype', 'spec', 'componentpol', 'off');
                data = data{1}';
                
            % select ica scalp maps
            % ---------------------
            case { 'scalp' 'scalpLaplac' 'scalpGrad' }
                for si = 1:size(STUDY.cluster(cluster_ind).sets,2)
                    idat  = STUDY.datasetinfo(STUDY.cluster(cluster_ind).sets(:,si)).index;
                    icomp = STUDY.cluster(cluster_ind).comps(si);
                    fprintf('Pre-clustering array row %d, adding interpolated scalp maps for dataset %d component %d...\n', si, idat, icomp);
                    data(si,:) = std_readtopo(ALLEEG, idat, icomp, scalpmodif, 'preclust');
                end
                
            % select ica equivalent dipole locations
            % --------------------------------------
            case 'dipoles'
                for si = 1:size(STUDY.cluster(cluster_ind).sets,2)
                    idat  = STUDY.datasetinfo(STUDY.cluster(cluster_ind).sets(1,si)).index;
                    icomp = STUDY.cluster(cluster_ind).comps(si);
                    fprintf('Pre-clustering array row %d, adding dipole for dataset %d component %d...\n', si, idat, icomp);
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
                        data(si,:) = ALLEEG(idat).dipfit.model(icomp).posxyz(ldip,:);
                    catch
                        error([ sprintf('Some dipole information is missing (e.g. component %d of dataset %d)', icomp, idat) 10 ...
                            'Components are not assigned a dipole if residual variance is too high so' 10 ...
                            'in the STUDY info editor, remember to select component by residual' 10 ...
                            'variance (column "select by r.v.") prior to preclustering them.' ]);
                    end
                end

            % select ica equivalent dipole moments
            % ------------------------------------
            case 'moments'
                for si = 1:size(STUDY.cluster(cluster_ind).sets,2)
                    idat  = STUDY.datasetinfo(STUDY.cluster(cluster_ind).sets(1,si)).index;
                    icomp = STUDY.cluster(cluster_ind).comps(si);
                    fprintf('Pre-clustering array row %d, adding dipole moment for dataset %d component %d...\n', si, idat, icomp);
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
                        posxyz = ALLEEG(idat).dipfit.model(icomp).posxyz(ldip,:);
                        momxyz = ALLEEG(idat).dipfit.model(icomp).momxyz(ldip,:);
                        if sum(posxyz .* momxyz) < 0, momxyz = -momxyz; end % make dipole pointing outward
                        data(si,:) = momxyz;
                    catch
                        error([ sprintf('Some dipole information is missing (e.g. component %d of dataset %d)', icomp, idat) 10 ...
                            'Components are not assigned a dipole if residual variance is too high so' 10 ...
                            'in the STUDY info editor, remember to select component by residual' 10 ...
                            'variance (column "select by r.v.") prior to preclustering them.' ]);
                    end
                end
                
            % cluster on ica ersp / itc values
            % --------------------------------
            case  {'ersp', 'itc' }
                [STUDY data] = std_readersp( STUDY, ALLEEG, 'design', NaN, 'clusters', cluster_ind, 'timerange', timewindow, 'freqrange', freqrange, 'measure', strcom);
                data = reshape(data{1}, [size(data{1},1)*size(data{1},2) size(data{1},3)])';
                
        end

        % adjust number of PCA values
        % ---------------------------
        if isnan(npca), npca = 5; end; % default number of components
        if npca >= size(data,2)
            % no need to run PCA, just copy the data.
            % But still run it to "normalize" coordinates
            % --------------------------------------
            npca = size(data,2);
        end;        
        if npca >= size(data,1) % cannot be more than the number of components
            npca = size(data,1);            
        end;        
        if ~strcmp(strcom, 'dipoles')
            fprintf('PCA dimension reduction to %d for command ''%s'' (normalize:%s; weight:%d)\n', ...
                npca, strcom, fastif(norm, 'on', 'off'), weight);
        else
            fprintf('Retaining the three-dimensional dipole locations (normalize:%s; weight:%d)\n', ...
                fastif(norm, 'on', 'off'), weight);
        end
        
        % run PCA to reduce data dimension
        % --------------------------------
        switch strcom
            case {'ersp','itc'}
                dsflag = 1;
                while dsflag
                    try,
                        clustdatatmp = runpca( double(data.'), npca, 1);
                        dsflag = 0;
                    catch,
                        % downsample frequency by 2 and times by 2
                        % ----------------------------------------
                        data = data(:,1:2:end);
                        %idat = STUDY.datasetinfo(STUDY.setind(1)).index; 
                        %[ tmp freqs times ] = std_readersp( ALLEEG, idat, succompind{1}(1));
                        %[data freqs times ] = erspdownsample(data,4, freqs,times,Ncond); 
                        if strcmp(varargin{index}(end-1) , 'downsample')
                            varargin{index}(end) = {celltomat(varargin{index}(end)) + 4};
                        else
                            varargin{index}(end+1) = {'downsample'};
                            varargin{index}(end+1) = {4};
                        end
                    end
                end
                clustdatatmp = clustdatatmp.';
            case 'dipoles'
                % normalize each cordinate by the std dev of the radii
                normval = std(sqrt(data(:,1).^2 + data(:,2).^2 + data(:,3).^2));
                clustdatatmp = data./normval;
                norm = 0;
            case 'erp'
                clustdatatmp = runpca( double(data.'), npca, 1);
                clustdatatmp = abs(clustdatatmp.');
            otherwise
                clustdatatmp = runpca( double(data.'), npca, 1);
                clustdatatmp = clustdatatmp.';
        end
        
        if norm %normalize the first pc std to 1
            normval = std(clustdatatmp(:,1));
            for icol = 1:size(clustdatatmp,2)
                clustdatatmp(:,icol) = clustdatatmp(:,icol) /normval;
            end
        end
        if weight ~= 1
            clustdata(:,end+1:end+size(clustdatatmp,2)) = clustdatatmp * weight;
        else
            clustdata(:,end+1:end+size(clustdatatmp,2)) = clustdatatmp;
        end
        
        if strcmpi(strcom, 'itc') || strcmpi(strcom, 'ersp')
            erspmode = 'already_computed';
        end
    end
    
    % Compute a second PCA of the already PCA'd data if there are too many PCA dimensions.
    % ------------------------------------------------------------------------------------
    if size(clustdata,2) > secondlevpca
        fprintf('Performing second-level PCA: reducing dimension from %d to %d \n', ...
                size(clustdata,2), secondlevpca);
        clustdata = runpca( double(clustdata.'), secondlevpca, 1);
        clustdata = clustdata.';
    end
    
    STUDY.etc.preclust.preclustdata   = clustdata;
    STUDY.etc.preclust.preclustparams = varargin;
    if isfield(STUDY.etc.preclust, 'preclustcomps')
        STUDY.etc.preclust = rmfield(STUDY.etc.preclust, 'preclustcomps');
    end
    
    % The preclustering level is equal to the parent cluster that the components belong to.
    if ~isempty(cluster_ind)
        STUDY.etc.preclust.clustlevel = cluster_ind; 
    else
        STUDY.etc.preclust.clustlevel = 1; % No parent cluster (cluster on all components in STUDY).
    end

return

% erspdownsample() - down samples component ERSP/ITC images if the
%        PCA operation in the clustering feature reduction step fails.
%        This is a helper function called by eeg_preclust().

function [dsdata, freqs, times] = erspdownsample(data, n, freqs,times,cond)
    len = length(freqs)*length(times); %size of ERSP
    nlen = ceil(length(freqs)/2)*ceil(length(times)/2); %new size of ERSP
    dsdata = zeros(size(data,1),cond*nlen);
    for k = 1:cond
        tmpdata = data(:,1+(k-1)*len:k*len);
        for l = 1:size(data,1) % go over components
            tmpersp = reshape(tmpdata(l,:)',length(freqs),length(times));
            tmpersp = downsample(tmpersp.', n/2).'; %downsample times
            tmpersp = downsample(tmpersp, n/2); %downsample freqs
            dsdata(l,1+(k-1)*nlen:k*nlen)  = tmpersp(:)';
        end
    end

% the function below checks the precense of the centroid field
function cluster = checkcentroidfield(cluster, varargin);
    for kk = 1:length(cluster)
        if ~isfield('centroid', cluster(kk)), cluster(kk).centroid = []; end
        for vi = 1:length(varargin)
            if isfield(cluster(kk).centroid, varargin{vi})
                cluster(kk).centroid = rmfield(cluster(kk).centroid, varargin{vi});
            end
        end
    end
    
    % rapid filtering for ERP (from std_plotcurve)
% -----------------------
function tmpdata2 = myfilt(tmpdata, srate, lowpass, highpass);
if ischar(highpass)
    highpass = str2num(highpass);
end
    bscorrect = 1;
    if bscorrect
        % Getting initial baseline
        bs_val1  =  mean(tmpdata,1);
        bs1      = repmat(bs_val1, size(tmpdata,1), 1);
    end
    
    % Filtering
    tmpdata2 = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3)*size(tmpdata,4));
    tmpdata2 = eegfiltfft(tmpdata2',srate, lowpass, highpass)';
    tmpdata2 = reshape(tmpdata2, size(tmpdata,1), size(tmpdata,2), size(tmpdata,3), size(tmpdata,4));
    
    if bscorrect
        % Getting after-filter baseline
        bs_val2  =  mean(tmpdata2,1);
        bs2      = repmat(bs_val2, size(tmpdata2,1), 1);
        
        % Correcting the baseline
        realbs = bs1-bs2;
        tmpdata2 = tmpdata2 + realbs;
    end
