% eeg_preclust() - prepare STUDY components' location and activity measures for later clustering.
%                  Selected measures (one or more from options: ERP, dipole locations, spectra,
%                  scalp maps, ERSP, and ITC) are computed for each dataset in the STUDY 
%                  set, unless they already present. After all requested measures are computed 
%                  and saved in the STUDY datasets, a PCA  matrix (by runica() with 'pca' option) 
%                  is constructed (this is the feature reduction step). This matrix will be used 
%                  as input to the clustering  algorithm in pop_clust(). eeg_preclust() allows 
%                  selection of a subset of components to use in the clustering. This subset 
%                  may be a user-specified component subset, components with dipole model residual 
%                  variance lower than a defined threshold (see dipfit()), or components from 
%                  an already existing cluster (for hierarchical clustering). The EEG datasets
%                  in the ALLEEG structure are updated. If new measures are added, the updated 
%                  EEG sets are also saved to disk. Called by pop_preclust(). Follow with 
%                  eeg_clust() or pop_clust().
% Usage:    
%                >> [ALLEEG,STUDY] = eeg_preclust(ALLEEG,STUDY); % cluster all comps in all sets
%                >> [ALLEEG,STUDY] = eeg_preclust(ALLEEG,STUDY,clustind,compind, preproc1,...);
%
% Required inputs:
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%
% Optional inputs:
%   clustind     - a cluster index for further clustering (hierarchical clustering), for example 
%                  to cluster a mu component cluster into left mu and right mu sub-clusters. 
%                  Should be empty for first stage (whole STUDY) clustering {default: []}
%
%   compind      - [ cell array | filename | [] ] Indices of dataset components to cluster. 
%                  A cell array of vectors of component indices for each dataset: Empty brackets 
%                  ([]) means all components in a dataset; 0 means none of them.
%                  For example: {[2:5] [10:15] [] [0]} includes components 2 to 5 from the first 
%                  dataset, components 10 to 15 from the second dataset, excludes dataset three, 
%                  but includes all components of the fourth dataset.
%                  Else, a quoted string containing the name of a .mat file containing such a cell 
%                  array variable named 'clustcomp' {default: [] => include all components 
%                  from all datasets}.
%
%   preprocX     - {'comnd' 'key1' val1 'key2' val2 ...} component clustering measures to prepare
%
%                  'comnd' = component measures to compute:
%                    * 'erp'     = cluster on the component ERPs,
%                    * 'dipoles' = cluster on the component [X Y Z] dipole locations
%                    * 'dipselect' = select components to cluster on that have residual 
%                                  dipole variance less than a specified threshold. 
%                    * 'spec'    = cluster on the component activity spectra 
%                                  (with the baseline removed).
%                    * 'scalp'   = cluster on component topoplot() scalp maps 
%                                  (or their absolute values),
%                    * 'scalpLaplac' = cluster on component topoplot() laplacian scalp maps
%                                  (or their absolute values),
%                    * 'scalpGrad' = cluster on the topoplot() scalp map gradients 
%                                  (or their absolute values),
%                    * 'ersp'    = cluster on components ERSP. (requires: 'cycles', 
%                                  'freqrange', 'padratio', 'timewindow', 'alpha').
%                    * 'itc'     = cluster on components ITC.(requires: 'cycles', 
%                                  'freqrange', 'padratio', 'timewindow', 'alpha').
%
%                  'key'   optional inputs used in computing  the specified measures:
%                    * 'npca'    =  [integer] number of PCA components (PCA dimension) of the 
%                                   selected data to retain for clustering. {default: 5}
%                    * 'norm'    =  [0|1] 1 normalizes the PCA components so the variance of 
%                                   first principal component is 1 (useful when using several 
%                                   preprocessing measures - 'ersp','scalp',...). {default: 1}
%                    * 'weight'  =  [integer] weight with respect to other preprocessing measures.
%                    * 'freqrange'  = [min max] freq. range (in Hz) for spectrum, 'ersp', 
%                                   and 'itc' options.  
%                    * 'timewindow' = [min max] time window (in sec) for 'erp' 'ersp', 
%                                   and 'itc' options.  
%                    * 'abso'    =  [0|1] 1 = take absolute values of topoplot(), Gradient, or 
%                                   Laplacian maps {default: 1}
%                    * 'cycles'  =  [0| cycles_factor] for ERSP and ITC (see timef() for details) 
%                                   {default: 0 (=> FFT method)}
%                    * 'padratio'=  [integer] for ERSP and ITC (see timef() for details) {default:1}
%                    * 'alpha'   =  [integer] bootstrap probability significance trhshold 
%                                   for ERSP and ITC (see timef() for details). {default: 0.01}
%                    * 'funarg'  =  [cell array] optional function arguments for mean spectrum 
%                                   calculation (see spectopo() for details) {default: none}
%                    * 'rv'      =  [number < 1] for 'dipselect', a threshold on the component 
%                                   residual variance. Only components with a lower residual variance 
%                                   (rv) will be clustered {default: 0 (all components)}
% Outputs:
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures modified by adding preprocessing 
%                  data (pointers to float files that hold ERSP, spectrum, etc. information).
%   STUDY        - the input STUDY set with added pre-clustering data, for use by pop_clust() 
%
% Example:
%   >> [ALLEEG  STUDY] = eeg_preclust(ALLEEG, STUDY, [], [] , { 'dipselect'  'rv'  0.15  } ,...
%                        { 'spec'  'npca' 10 'norm' 1 'weight' 1 'freqrange'  [ 3 25 ] } , ...
%                        { 'erp'   'npca' 10 'norm' 1 'weight' 2 'timewindow' [ 350 500 ] } ,...
%                        { 'scalp' 'npca' 10 'norm' 1 'weight' 2 'abso' 1 } , ...
%                        { 'dipoles'         'norm' 1 'weight' 15 } , ...
%                        { 'ersp'  'npca' 10 'freqrange' [ 3 25 ] 'cycles' [ 3 0.5 ] 'alpha' 0.01 ....
%                                  'padratio' 4 'timewindow' [ -1600 1495 ] 'norm' 1 'weight' 1 } ,...
%                        { 'itc'   'npca' 10 'freqrange' [ 3 25 ] 'cycles' [ 3 0.5 ] 'alpha' 0.01 ...
%                                  'padratio' 4 'timewindow' [ -1600 1495 ] 'norm' 1 'weight' 1 });
%                          
%                        % This prepares for first-stage clustering all components in the STUDY
%                        % datasets except components with dipole model residual variance above 0.15.
%                        % Clustering will be based on the components' mean spectra in % the [3 25] Hz 
%                        % frequency range, the components' ERPs in the [350 500] ms % time window, 
%                        % the (absolute-value) component scalp maps, their equivalent dipole locations,
%                        % and their ERSP and ITC images. 
%
% Authors: Hilit Serby, Arnaud Delorme & Scott Makeig, SCCN, INC, UCSD, May 13, 2004

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, May 13,2004, hilit@sccn.ucsd.edu
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

function [ ALLEEG, STUDY ] = eeg_preclust(ALLEEG, STUDY, cluster_ind, components_ind, varargin)
    
    if nargin < 2
        help eeg_preclust;
        return;
    end;
    if nargin == 2
        cluster_ind = []; % default to clustering the whole STUDY 
    end    
    if nargin == 3
        components_ind = []; % default to clustering all components 
    end    

    Ncond = length(STUDY.condition);
    if Ncond == 0
        Ncond = 1;
    end
    
    % Get component indices that are part of the cluster 
    % --------------------------------------------------
    if ~isempty(cluster_ind)
        if length(cluster_ind) ~= 1
            error('Only one cluster can be sub-clustered. To sub-cluster multiple clusters, first merge them.');
        end
        for k = 1 : size(STUDY.setind,2)   % go over the sets from the first condition (if there are some)
            sind = find(STUDY.cluster(cluster_ind).sets(1,:) == k); % the component indices that belong to 
                                                                    % the dataset in the cluster 
            if isempty(sind)
                succompind{k} = 0;
            else
                succompind{k} = STUDY.cluster(cluster_ind).comps(sind);
            end
        end
    end

    if ~isempty(components_ind) % if there is a component selection
        if isstr(components_ind) % if input is a .mat file, load component indices
            try, 
                eval( [ 'load ' components_ind ]);
            catch
                error('The compind argument is not a valid filename')
            end
            if isempty('clustcomp')
                error('The compind .mat file must have a cell array variable ''clustcomp''');
            end
            components_ind = clustcomp;
        end
        if length(components_ind) ~= size(STUDY.setind,2)
            error('Size of cell array of component indices (''compsind'') ~= (subjects * sessions) in STUDY.');
        end;
        if ~isempty(cluster_ind) % components to cluster on must be both part of the specified cluster 
                                 % and selected components
            for ind = 1:size(STUDY.setind,2)
                if isempty(components_ind{ind})
                    seti = STUDY.datasetinfo(STUDY.setind(1,ind)).index;
                    components_ind{ind} = 1:size(ALLEEG(seti).icawinv,2);
                end
                succompind{ind} = intersect(components_ind{ind}, succompind{ind});
                if isempty(succompind{ind})
                    succompind{ind} = 0;
                end
            end
        else
            succompind = components_ind;
        end
    else % no component selection 
        if ~exist('succompind')
            for ind = 1:size(STUDY.setind,2), succompind{ind} = []; end;
        end
    end
    
    for ind = 1:size(STUDY.setind,2)
        if isempty(succompind{ind})
            seti = STUDY.datasetinfo(STUDY.setind(1,ind)).index;
            succompind{ind} = 1:size(ALLEEG(seti).icawinv,2);
        else
            succompind{ind} = succompind{ind}(find(succompind{ind}));% remove zeros
            succompind{ind} = sort(succompind{ind}); % sort the components
        end
    end;
    
    % Scan which component to remove (no dipole info)
    % -----------------------------------------------
    update_flag = 0;
    clear rv
    for index = 1:length(varargin) % scan commands
        strcom = varargin{index}{1};
        if strcmpi(strcom(1:3),'dip')
            update_flag = 1;
            if strcmpi(strcom,'dipoles') & ~exist('rv')
                rv = 1;
            elseif strcmpi(strcom,'dipselect')
                rv = varargin{index}{3};
                varargin(index) = []; %remove this command
                break
            end
        end        
    end
    if update_flag % dipole information is used to select components
        for si = 1:size(STUDY.setind,2)% scan datasets that are part of STUDY
            idat = STUDY.datasetinfo(STUDY.setind(1,si)).index;
            if isfield(ALLEEG(idat).dipfit, 'rv')
                indrm = []; % components that will be removed from the clustering
                indleft = []; % components that are left in clustering
                for icomp = succompind{si} % scan components
                    if (ALLEEG(idat).dipfit.model(icomp).rv >= rv) | isnan(ALLEEG(idat).dipfit.model(icomp).rv) 
                        % Components have rv bigger than asked for 
                        % fprintf('Component %d of dataset %d has no dipole info and was removed\n', icomp, idat)
                        indrm = [indrm icomp];
                    else
                        indleft = [indleft icomp];
                    end;
                end;
                if ~isempty(indrm)
                    succompind{si} = indleft;
                end
            else
                disp('No dipole information found in the data, using all dipoles')
            end
        end;
        % Create a cluster of removed components 
        % --------------------------------------
        update_flag = 0;
        if isempty(cluster_ind) % first step clustering
            for si = 1:size(STUDY.setind,2)
                idat = STUDY.datasetinfo(STUDY.setind(1,si)).index;
                % Check if not all components in this dataset are part of the clustering
                rmind = setdiff([1:size(ALLEEG(idat).icawinv,2)], succompind{si});
                if ~isempty(rmind) 
                    rmcomp{si} = rmind;
                    update_flag = 1;
                end
            end
        else % cluster on a specific cluster components
            for si = 1:size(STUDY.setind,2)
                % Check if not all components in the cluster are part of the clustering
                compind = find(STUDY.cluster(cluster_ind).sets(1,:) == si);
                rmind = setdiff(STUDY.cluster(cluster_ind).comps(compind), succompind{si});  
                if ~isempty(rmind) 
                    rmcomp{si} = rmind;
                    update_flag = 1;
                end
            end
        end
        if update_flag % Update STUDY with new components
            if isempty(cluster_ind)
                a =  ['% A subset of components with dipole information were readied for clustering.'...
                        'Other components were placed in a separate cluster'];
                [STUDY] = cls_createclust(STUDY, ALLEEG, 'Notclust');
                STUDY.cluster(end).parent{end} = STUDY.cluster(1).name; 
                STUDY.cluster(1).child{end+1} = STUDY.cluster(end).name;
            else
                a =  ['% A subset of components from ' STUDY.cluster(cluster_ind).name  'with dipole information, were readied for clustering.' ...
                    'Other components were placed in a separate cluster'];    
                [STUDY] = cls_createclust(STUDY, ALLEEG, [ 'Notclust ' num2str(cluster_ind) ] );
                STUDY.cluster(end).parent{end} = STUDY.cluster(cluster_ind).name;
                STUDY.cluster(cluster_ind).child{end+1} = STUDY.cluster(end).name;
            end
            for k = 1: length(rmcomp)
                if ~isempty(rmcomp{k})
                    STUDY.cluster(end).sets = [STUDY.cluster(end).sets k*ones(1,length(rmcomp{k}))];
                    STUDY.cluster(end).comps = [STUDY.cluster(end).comps rmcomp{k}];
                end
            end
            if Ncond > 1
                tmp = ones(Ncond, length(STUDY.cluster(end).sets));
                for l = 1:Ncond
                    tmp(l,:) = STUDY.cluster(end).sets + (l-1)*size(STUDY.setind,2);
                end
                STUDY.cluster(end).sets = tmp;
                clear tmp
			end

            STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
        end;
    end;
    
    % scan all commands
    % -----------------
    clustdata = [];
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
        cycles = 0;
        padratio  = 1;
        alpha = 0.01;
        fun_arg = [];
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
                case 'cycles'
                    cycles = varargin{index}{subind+1};
                case 'alpha'
                    alpha = varargin{index}{subind+1};
                case 'padratio'
                    padratio = varargin{index}{subind+1};
                otherwise 
                    fun_arg{length(fun_arg)+1} =  varargin{index}{subind+1};
            end
        end
        
        % If ersp / itc values already exist in some of the datasets,
        % make sure they are the same as the requested parameters.
        % -----------------------------------------------------------
        if (strcmpi(strcom,'ersp')  | strcmpi(strcom,'itc') )  
            overwrite = 0; 
            params = [];
            for si = 1:size(STUDY.setind,2) % scan datasets that are part of STUDY
                for cond = 1:Ncond
                    if overwrite == 0
                        idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;
                        if isfield(ALLEEG(idat).etc, 'icaerspparams')
                            params = ALLEEG(idat).etc.icaerspparams;
                        end
                        if ~isempty(params)
                            overwrite = 2;
                            if  sum(params.cycles ~= cycles) | sum(params.freqrange ~= freqrange)  | ...
                                    (padratio ~= params.padratio) | (alpha~= params.alpha) 
                                set_yes =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_no''), ''value'', 0);'];
                                set_no =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_yes''), ''value'', 0);' ];
                                ersp_ans = inputgui({[1] [1] [1] [1] [1 1] [1]}, ...
                                       { {'style' 'text' 'string' [upper(strcom) ' info exists for dataset: '  ALLEEG(idat).filename '. It does not fit requested values.' ] } ...
                                         {'style' 'text' 'string' ['Existing values are: frequency range - [' num2str(params.freqrange) '], wavelet cycles - [' num2str(params.cycles) ...
                                                 '], padratio - ' num2str(params.padratio) ', and bootstrap significance - ' num2str(params.alpha) ] } {} ...
                                         {'style' 'text' 'string' ['Would you like to recalculate ' upper(strcom) ' and overwrite those values?' ]} ...
                                         {'style' 'checkbox' 'tag' 'ersp_yes' 'string' 'Yes' 'value' 1 'Callback' set_yes }  ...
                                         {'style' 'checkbox' 'tag' 'ersp_no' 'string' 'Use existing ERSP info' 'value' 0 'Callback' set_no } {} },...
                                       [], ['Recalculate ' upper(strcom) ' parameters -- part of cls_ersp()']); 
                                if find(celltomat(ersp_ans))  == 2 % use existing ERSP info from this dataset
                                    overwrite = 2;
                                    cycles = params.cycles;
                                    freqrange = params.freqrange;
                                    alpha = params.alpha;
                                    padratio = params.padratio;
                               else % Over write data in dataset
                                    overwrite = 1;
                                end
                            end
                        end
                    end
                end
            end
        end

        
        % scan datasets
        % -------------
        data = [];
        for si = 1:size(STUDY.setind,2) 
            switch strcom
             
             % select ica ERP
             % --------------
             case 'erp'    , 
                  if ~isempty(succompind{si})
                      for cond = 1 : Ncond
                          idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                          if ALLEEG(idat).trials == 1
                              error('No epochs in dataset: ERP information has no meaning');
                         end
                         if ~isempty('timewindow')
                             [tmp, X, t] = cls_erp(ALLEEG(idat), succompind{si}, timewindow);
                             if ~isempty(tmp)
                                 ALLEEG(idat).etc = tmp;
                             end
                         else
                             [tmp, X, t] = cls_erp(ALLEEG(idat), succompind{si});
                             if ~isempty(tmp)
                                 ALLEEG(idat).etc = tmp;
                             end
                         end
                         if cond == 1
                             con_data = abs(X);
                             con_t = t;
                         else % concatenate across conditions
                             con_t = [con_t; t];
                             con_data = [con_data abs(X)];
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
             case 'scalp' , % NB: scalp maps are identical across conditions (within session)
                 idat = STUDY.datasetinfo(STUDY.setind(1,si)).index; 
                 if ~isempty(succompind{si})
                     [tmp, X] = cls_scalp(ALLEEG(idat), succompind{si});
                     if ~isempty(tmp)
                         ALLEEG(idat).etc = tmp;
                     end
                     if abso % absolute values
                         data = [ data; abs(X) ];
                     else
                         data = [ data; X ];
                     end
                     clear X tmp;
                 end       
             % select Laplacian ica comp scalp maps
             % ------------------------------------
             case 'scalpLaplac' , 
                 idat = STUDY.datasetinfo(STUDY.setind(1,si)).index; 
                 if ~isempty(succompind{si})
                     [tmp, X] = cls_scalpL(ALLEEG(idat), succompind{si}); 
                     if ~isempty(tmp)
                         ALLEEG(idat).etc = tmp;
                     end
                     if abso
                         data = [ data; abs(X)];
                     else
                         data = [ data; X];
                     end
                     clear X tmp;
                 end
             % select Gradient ica comp scalp maps
             % -----------------------------------
             case 'scalpGrad'   , 
                 idat = STUDY.datasetinfo(STUDY.setind(1,si)).index; 
                 if ~isempty(succompind{si})
                     [tmp, X] = cls_scalpG(ALLEEG(idat), succompind{si}); 
                     if ~isempty(tmp)
                         ALLEEG(idat).etc = tmp;
                     end
                     if abso
                         data = [ data; abs(X)];
                     else
                         data = [ data; X];
                     end
                     clear X tmp;
                 end 
             % select ica comp spectra
             % -----------------------
             case 'spec'   , 
                 if si == 1
                    overwrite = 0;
                end
                if ~isempty(succompind{si})
                    for cond = 1 : Ncond 
                         idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                         if ~isempty('freqrange')
                             [tmp, X, f,overwrite] = cls_spec(ALLEEG(idat),succompind{si}, ...
                                                                  freqrange, fun_arg,overwrite);
                             if ~isempty(tmp)
                                 ALLEEG(idat).etc = tmp;
                             end
                         else
                             [tmp, X, f,overwrite] = cls_spec(ALLEEG(idat),succompind{si}, ...
                                                                         [], fun_arg,overwrite);
                             if ~isempty(tmp)
                                 ALLEEG(idat).etc = tmp;
                             end
                         end
                         if cond == 1
                             con_data = X;
                             con_f = f;
                         else % concatenate across conditions
                             con_f = [con_f; f];
                             con_data = [con_data X];
                         end
                     end
                     if isempty(data)
                          frequencies = con_f;
                      else
                          frequencies = [ frequencies; con_f];
                      end
                      data = [ data; con_data ];
                      clear f X con_f con_data tmp;  
                  end               
             % select dipole information
             % -------------------------
             case 'dipoles' % NB: dipoles are identical across conditions (within session)
                            % (no need to scan across conditions)
              idat = STUDY.datasetinfo(STUDY.setind(1,si)).index; 
              count = size(data,1)+1;
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
             % cluster on ica ersp / itc values
             % --------------------------------
             case  {'ersp', 'itc'}
                type =  strcom;            
                if ~isempty(succompind{si})
                    idattot = [];
                    for cond = 1 : Ncond 
                         idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                         idattot = [idattot idat];
                         % compute ERSP/ ITC, if doesn't exist.
                         [tmp] = cls_ersp(ALLEEG(idat),succompind{si}, freqrange, timewindow, ...
                                                                     cycles, padratio, alpha, type);
                         if ~isempty(tmp)
                             ALLEEG(idat).etc = tmp;
                         end
                    end
                    % prepare ERSP / ITC data for clustering (select requested components, 
                    % mask and change to a common base if multiple conditions).
                    % --------------------------------------------------------------------
                    X = cls_ersptocluster(ALLEEG(idattot),succompind{si}, timewindow, freqrange, type); 
                    data = [data; X];
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
                        clustdatatmp = runpca( data.', npca, 1);
                        dsflag = 0;
                    catch,
                        % downsample frequency by 2 and times by 2
                        % ----------------------------------------
                        idat = STUDY.datasetinfo(STUDY.setind(1)).index; 
                        times = ALLEEG(idat).etc.icaerspparams.times;
                        freqs = ALLEEG(idat).etc.icaerspparams.logfreqs;
                        [data, freqs, times] = erspdownsample(data,4, freqs,times,Ncond); 
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
            otherwise
                clustdatatmp = runpca( data.', npca, 1);
                clustdatatmp = clustdatatmp.';
        end
        
        clear data;
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
        
    end;
    
    % Compute a second PCA of the already PCA'd data if there are too many PCA dimensions.
    % ------------------------------------------------------------------------------------
    if size(clustdata,2) > 10
        fprintf('Performing second-level PCA: reducing dimension from %d to 10 \n', size(clustdata,2));
        clustdata = runpca( clustdata.', 10, 1);
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
