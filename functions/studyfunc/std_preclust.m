%   std_preclust() - select measures to be included in computation of a preclustering array.
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

function [ STUDY, ALLEEG ] = std_preclust(STUDY, ALLEEG, cluster_ind, varargin)
    
    if nargin < 2
        help std_preclust;
        return;
    end;
    
    if nargin == 2
        cluster_ind = 1; % default to clustering the whole STUDY 
    end    

    % check dataset length consistency before computing
    % -------------------------------------------------
    pnts  = ALLEEG(STUDY.datasetinfo(1).index).pnts;
    srate = ALLEEG(STUDY.datasetinfo(1).index).srate;
    xmin  = ALLEEG(STUDY.datasetinfo(1).index).xmin;
    for index = 1:length(STUDY.datasetinfo)
        ind = STUDY.datasetinfo(index).index;
        if pnts  ~= ALLEEG(ind).pnts, error(sprintf('Dataset %d does not have the same number of point as dataset 1', ind)); end;
        if srate ~= ALLEEG(ind).srate, error(sprintf('Dataset %d does not have the same sampling rate as dataset 1', ind)); end;
        if xmin  ~= ALLEEG(ind).xmin, warning(sprintf('Dataset %d does not have the same time limit as dataset 1', ind)); end;
    end;
    
    % Get component indices that are part of the cluster 
    % --------------------------------------------------
    if isempty(cluster_ind)
        cluster_ind = 1;
    end;
    if length(cluster_ind) ~= 1
        error('Only one cluster can be sub-clustered. To sub-cluster multiple clusters, first merge them.');
    end
    
    % the goal of this code below is to find the components in the cluster
    % of interest for each set of condition 
    for k = 1:size(STUDY.setind,2)
        % Find the first entry in STUDY.setind(:,k) that is non-NaN. We only need one since
        % components are the same across conditions.
        for ri = 1:size(STUDY.setind,1)
            if ~isnan(STUDY.setind(ri,k)), break; end
        end
        sind = find(STUDY.cluster(cluster_ind).sets(ri,:) == STUDY.setind(ri,k));
        succompind{k} = STUDY.cluster(cluster_ind).comps(sind);
    end;
    for ind = 1:size(STUDY.setind,2)
        succompind{ind} = succompind{ind}(find(succompind{ind})); % remove zeros 
                                                                  % (though there should not be any? -Arno)
        succompind{ind} = sort(succompind{ind}); % sort the components
    end;
    
    % scan all commands to see if file exist
    % --------------------------------------
    for index = 1:length(varargin) % scan commands
        strcom = varargin{index}{1};
        if strcmpi(strcom, 'scalp') | strcmpi(strcom, 'scalplaplac') | strcmpi(strcom, 'scalpgrad') 
            strcom = 'topo';
        end;
        if ~strcmpi(strcom, 'dipoles') & ~strcmpi(strcom, 'finaldim')
            tmpfile = fullfile( ALLEEG(1).filepath, [ ALLEEG(1).filename(1:end-3) 'ica' strcom ]); 
            if exist(tmpfile) ~= 2
                error([ 'Could not find "' upper(strcom) '" file, they must be precomputed' ]);
            end;
        end;
    end;
            
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
    end;
    
    % scan all commands
    % -----------------
    clustdata = [];
    erspquery = 0;
    for index = 1:length(varargin)
        
        % decode inputs
        % -------------
        strcom = varargin{index}{1};
        if any(strcom == 'X'), disp('character ''X'' found in command'); end;
        %defult values
        npca = NaN;
        norm = 1;
        weight = 1;
        freqrange = [];
        timewindow = [];
        abso = 1;
        fun_arg = [];
        savetrials = 'off';
        recompute  = 'on';
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
        data = [];
        if strcmpi(strcom, 'scalp'),           scalpmodif = 'none';
        elseif strcmpi(strcom, 'scalpLaplac'), scalpmodif = 'laplacian';
        else                                   scalpmodif = 'gradient';
        end;
        
        for si = 1:size(STUDY.setind,2)
            switch strcom
             
             % select ica component ERPs
             % -------------------------
             case 'erp', 
                  for kk = 1:length(STUDY.cluster)
                      if isfield(STUDY.cluster(kk).centroid, 'erp')
                          STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'erp');
                          STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'erp_times');
                      end;
                  end;
                  if ~isempty(succompind{si})
                      for cond = 1 : size(STUDY.setind,1)
                         if ~isnan(STUDY.setind(cond,si))
                            idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                            if ALLEEG(idat).trials == 1
                               error('No epochs in dataset: ERP information has no meaning');
                            end
                            con_t = []; con_data = [];
                            [X t] = std_readerp( ALLEEG, idat, succompind{si}, timewindow);

                            X = abs(X); % take the absolute value of the ERP
                            
                            STUDY.preclust.erpclusttimes = timewindow;
                            if cond == 1
                                con_data = X;
                                con_t = t;
                            else % concatenate across conditions
                                con_t = [con_t; t];
                                con_data = [con_data X];
                            end
                         end
                     end
                     if isempty(data)
                          times = con_t;
                     else
                          times = [ times ; con_t];
                     end
                     data = [ data; con_data ];
                     clear X t con_data con_t tmp;
                 end

              % select ica scalp maps
              % --------------------------
             case { 'scalp' 'scalpLaplac' 'scalpGrad' }

                 for cond = 1:size(STUDY.setind,1)   % Find first nonNaN index
                    if ~isnan(STUDY.setind(cond,si)), break; end
                 end       
                 idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;
                 fprintf('Computing/loading interpolated scalp maps for dataset %d...\n', idat);
                 if ~isempty(succompind{si})
                    X = std_readtopo(ALLEEG, idat, succompind{si}, scalpmodif, 'preclust');

                    if abso % absolute values
                       data = [ data; abs(X) ];
                    else
                       data = [ data; X ];
                    end
                    clear X tmp;
                 end
                                     
             % select ica comp spectra
             % -----------------------
             case 'spec', 
                 for kk = 1:length(STUDY.cluster)
                     if isfield(STUDY.cluster(kk).centroid, 'spec')
                         STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'spec');
                         STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'spec_freqs');
                     end;
                 end;
                 if ~isempty(succompind{si})
                     for cond = 1 : size(STUDY.setind,1)
                         if ~isnan(STUDY.setind(cond,si))
                            idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index; 
                            con_f = []; con_data = [];
                            [X, f] = std_readspec(ALLEEG, idat, succompind{si}, freqrange);
                            STUDY.preclust.specclustfreqs = freqrange;
                            if cond == 1
                                con_data = X;
                                con_f = f;
                            else % concatenate across conditions
                                con_f = [con_f; f];
                                con_data = [con_data X];
                            end
                         end
                     end
                     if isempty(data)
                          frequencies = con_f;
                     else
                         frequencies = [ frequencies; con_f];
                     end
                     con_data = con_data - repmat(mean(con_data,2), [1 size(con_data,2)]);
                     con_data = con_data - repmat(mean(con_data,1), [size(con_data,1) 1]);
                     data = [ data; con_data ];
                     clear f X con_f con_data tmp;  
                 end               
                 
             % select dipole information
             % -------------------------
             case 'dipoles' % NB: dipoles are identical across conditions (within session)
                            % (no need to scan across conditions)
              for cond = 1:size(STUDY.setind,1)  % Scan for first nonNaN index.
                 if ~isnan(STUDY.setind(cond,si)), break; end
              end
              idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
              count = size(data,1)+1;
              try
                 for icomp = succompind{si}
                 % select among 3 sub-options
                 % --------------------------
                 if ~isempty(succompind{si})
                    ldip = 1;
                    if size(ALLEEG(idat).dipfit.model(icomp).posxyz,1) == 2 % two dipoles model
                        if any(ALLEEG(idat).dipfit.model(icomp).posxyz(1,:)) ...
                            & any(ALLEEG(idat).dipfit.model(icomp).posxyz(2,:)) %both dipoles exist
                           % find the leftmost dipole
                           [garb ldip] = max(ALLEEG(idat).dipfit.model(icomp).posxyz(:,2)); 
                        elseif any(ALLEEG(idat).dipfit.model(icomp).posxyz(2,:)) 
                           ldip = 2; % the leftmost dipole is the only one that exists
                        end
                     end
                     data(count,:) = ALLEEG(idat).dipfit.model(icomp).posxyz(ldip,:);
                     count = count+1;
                 end
                 end 
              catch
                error([ sprintf('Some dipole information is missing (e.g. component %d of dataset %d)', icomp, idat) 10 ...
                              'Components are not assigned a dipole if residual variance is too high so' 10 ...
                              'in the STUDY info editor, remember to select component by residual' 10 ...
                              'variance (column "select by r.v.") prior to preclustering them.' ]);
              end
              
             % cluster on ica ersp / itc values
             % --------------------------------
             case  {'ersp', 'itc' }
                for kk = 1:length(STUDY.cluster)
                    if isfield(STUDY.cluster(kk).centroid, 'ersp')
                        STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'ersp');
                        STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'ersp_freqs');
                        STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'ersp_times');
                    end;
                    if isfield(STUDY.cluster(kk).centroid, 'itc')
                        STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'itc');
                        STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'itc_times');
                    end;
                end;
                if ~isempty(succompind{si})
                    idattot = [];
                    for cond = 1:size(STUDY.setind,1)
                        if ~isnan(STUDY.setind(cond,si))
                            idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                            idattot = [idattot idat];
                        end
                    end
                    STUDY.preclust.erspclustfreqs = freqrange;
                    STUDY.preclust.erspclusttimes = timewindow;
                    
                    % prepare ERSP / ITC data for clustering (select requested components, 
                    % mask and change to a common base if multiple conditions).
                    % --------------------------------------------------------------------
                    for k = 1:length(succompind{si})
                        if strcmpi(strcom, 'ersp')
                            tmp = std_readersp(ALLEEG, idattot, succompind{si}(k), timewindow, freqrange); 
                        else
                            tmp = real(std_readitc( ALLEEG, idattot, succompind{si}(k), timewindow, freqrange)); 
                        end;
                        if k == 1
                            X = zeros(length(succompind{si}), size(tmp,1), size(tmp,2), size(tmp,3)); 
                        end;
                        X(k,:,:,:) = tmp;
                    end;
                    data = [data; reshape(X, size(X,1), size(X,2)*size(X,3)*size(X,4)) ];
                    clear tmp X 
                end
                
             % execute custom command
             % ----------------------
             otherwise, 
                 if ~isempty(succompind{si})
                      if ~any(strcom == 'X'), strcom = [ 'X=' strcom ';' ]; end; 
                      EEG = ALLEEG(idat);
                      eval(strcom);
                      if length(findstr(strcom,'spectopo')) == 1
                          if isempty(data)
                              frequencies = f;
                          else
                              frequencies = [ frequencies; f];
                          end    
                      end
                      if size(X,1) == 1
                          data(end+1,:) = X;
                      else
                          if size(X,1) > length(succompind{si})
                              data(end+1:end+length(succompind{si}),:) = X(succompind{si},:);
                          elseif size(X,1) == length(succompind{si})
                              data(end+1:end+size(X,1),:) = X;
                          else
                              data(end+1:end+size(X,1),:) = X;
                              disp('Warning: wrong number of components (rows) generated by string command.');
                          end;
                      end;
                  end
            end;

            
        end; % end scan datasets

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
            end;
        end;
        if weight ~= 1
            clustdata(:,end+1:end+size(clustdatatmp,2)) = clustdatatmp * weight;
        else
            clustdata(:,end+1:end+size(clustdatatmp,2)) = clustdatatmp;
        end
        
        if strcmpi(strcom, 'itc') | strcmpi(strcom, 'ersp')
            erspmode = 'already_computed';
        end;
    end
    
    % Compute a second PCA of the already PCA'd data if there are too many PCA dimensions.
    % ------------------------------------------------------------------------------------
    if size(clustdata,2) > secondlevpca
        fprintf('Performing second-level PCA: reducing dimension from %d to %d \n', ...
                size(clustdata,2), secondlevpca);
        clustdata = runpca( double(clustdata.'), secondlevpca, 1);
        clustdata = clustdata.';
    end
    
    STUDY.etc.preclust.preclustdata = clustdata;
    STUDY.etc.preclust.preclustparams = varargin;
    STUDY.etc.preclust.preclustcomps = succompind;

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
