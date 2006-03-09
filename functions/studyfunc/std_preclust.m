% std_preclust() - prepare STUDY component location and activity measures for later clustering.
%                  Selected measures (one or more from options: ERPs, dipole locations, spectra,
%                  scalp maps, ERSPs, and ITCs) are computed for each dataset in the STUDY 
%                  set, unless they already present. After all requested measures are computed 
%                  and saved in the STUDY datasets, each feature dimension is reduced by computing 
%                  a PCA  decomposition. These PCA matrices (one per measure) are concatenated and 
%                  used as input to the clustering  algorithm in pop_clust(). std_preclust() allows 
%                  selection of a subset of components to use in the clustering. This subset 
%                  may be a user-specified component subset, components with dipole model residual 
%                  variance lower than a defined threshold (see dipfit()), or components from 
%                  an already existing cluster (for hierarchical clustering). The EEG datasets
%                  in the ALLEEG structure are updated. If new measures are added, the updated 
%                  EEG sets are also saved to disk. Called by pop_preclust(). Follow with 
%                  eeg_clust() or pop_clust().
% Usage:    
%                >> [ALLEEG,STUDY] = std_preclust(ALLEEG,STUDY); % cluster all comps in all sets
%                >> [ALLEEG,STUDY] = std_preclust(ALLEEG,STUDY,clustind,compind, preproc1,...);
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
%                    * 'dipselect' = select components to cluster that have residual 
%                                  dipole variance less than a specified threshold. 
%                    * 'spec'    = cluster on the component log activity spectra (in dB)
%                                  (with the baseline mean dB spectrum subtracted).
%                    * 'scalp'   = cluster on component (topoplot()) scalp maps 
%                                  (or on their absolute values),
%                    * 'scalpLaplac' = cluster on component (topoplot()) laplacian scalp maps
%                                  (or on their absolute values),
%                    * 'scalpGrad' = cluster on the (topoplot()) scalp map gradients 
%                                  (or on their absolute values),
%                    * 'ersp'    = cluster on components ERSP. (requires: 'cycles', 
%                                  'freqrange', 'padratio', 'timewindow', 'alpha').
%                    * 'itc'     = cluster on components ITC.(requires: 'cycles', 
%                                  'freqrange', 'padratio', 'timewindow', 'alpha').
%                    * 'finaldim' = final number of dimensions. Enables second-level PCA. 
%                                  By default this command is not used (see Example below).
%
%                  'key'   optional inputs used in computing  the specified measures:
%                    * 'npca'    =  [integer] number of principal components (PCA dimension) of 
%                                   the selected measures to retain for clustering. {default: 5}
%                    * 'norm'    =  [0|1] 1 -> normalize the PCA components so the variance of 
%                                   first principal component is 1 (useful when using several 
%                                   clustering measures - 'ersp','scalp',...). {default: 1}
%                    * 'weight'  =  [integer] weight with respect to other clustering measures.
%                    * 'freqrange'  = [min max] frequency range (in Hz) to include in activity 
%                                   spectrum, 'ersp', and 'itc' measures.  
%                    * 'timewindow' = [min max] time window (in sec) to include in 'erp',
%                                   'ersp', and 'itc' measures.  
%                    * 'abso'    =  [0|1] 1 = use absolute values of topoplot(), gradient, or 
%                                   Laplacian maps {default: 1}
%                    * 'cycles'  =  [0| cycles_factor] for ERSP and ITC (see >> timef details) 
%                                   {default: 0 (=> FFT method)}
%                    * 'padratio'=  [integer] for ERSP and ITC (see >> timef details) {default:1}
%                    * 'alpha'   =  [integer] bootstrap probability significance threshold for 
%                                   masking component mean ERSP and ITC measures 
%                                   (>> timef details) {default: 0.01}
%                    * 'funarg'  =  [cell array] optional function arguments for mean spectrum 
%                                   calculation (>> help spectopo) {default: none}
%                    * 'rv'      =  [number < 1] for dipole locations ('dipselect'), max component 
%                                   model residual variance. Only components with a lower residual 
%                                   variance (rv) will be clustered {default: 0 (all components)}
% Outputs:
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures, modified by adding preprocessing 
%                  data as pointers to float files that hold ERSPs, spectra, and/or other measures.
%   STUDY        - the input STUDY set with pre-clustering data added, for use by pop_clust() 
%
% Example:
%   >> [ALLEEG  STUDY] = std_preclust(ALLEEG, STUDY, [], [] , { 'dipselect'  'rv'  0.15  } ,...
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
%                        % This prepares, for initial clustering, all components in the STUDY
%                        % datasets except components with dipole model residual variance above 0.15.
%                        % Clustering will be based on the components' mean spectra in the [3 25] Hz 
%                        % frequency range, on the components' ERPs in the [350 500] ms time window, on
%                        % the (absolute-value) component scalp maps, on the equivalent dipole locations,
%                        % and on the mean ERSP and ITC images. 
%                        % The final keyword asks for a second level PCA dimension reduction (to 10
%                        % principal dimensions) to be performed. See the web tutorial for details.
%
% Authors: Arnaud Delorme, Hilit Serby & Scott Makeig, SCCN, INC, UCSD, May 13, 2004

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

% $Log: not supported by cvs2svn $
% Revision 1.27  2006/03/08 21:09:26  arno
% checking study
%
% Revision 1.26  2006/03/08 20:30:41  arno
% rename func
%
% Revision 1.25  2006/03/07 22:35:14  arno
% catch error when scanning for dipoles
%
% Revision 1.24  2006/03/06 23:47:47  arno
% computing gradient or laplacian
%
% Revision 1.23  2006/03/03 23:02:37  arno
% convert to doulbe before runing PCA
%
% Revision 1.22  2006/02/23 19:14:21  scott
% help msg
%
% Revision 1.21  2006/02/22 21:17:01  arno
% same
%
% Revision 1.20  2006/02/22 21:16:11  arno
% better decoding of arguments
%
% Revision 1.19  2006/02/22 21:05:02  arno
% update call to checkstudy
%
% Revision 1.18  2006/02/22 21:02:25  arno
% Removing cluster information
%
% Revision 1.17  2006/02/22 20:03:20  arno
% fixing dimension reduction
%
% Revision 1.16  2006/02/22 19:55:36  arno
% second level pca
%
% Revision 1.15  2006/02/18 01:01:38  arno
% eeg_preclust -> std_preclust
%
% Revision 1.14  2006/02/18 00:56:34  arno
% showing more warnings
%
% Revision 1.13  2006/02/16 22:10:16  arno
% scalpL -> scalpl
%
% Revision 1.12  2006/02/16 22:09:43  arno
% scalpG -> scalpg
%
% Revision 1.11  2006/02/11 00:29:44  arno
% final fix
%

function [ STUDY, ALLEEG ] = std_preclust(STUDY, ALLEEG, cluster_ind, components_ind, varargin)
    
    if nargin < 2
        help std_preclust;
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
    else % no component selection (use pre-selected components)
        if ~isfield(STUDY.datasetinfo, 'comps')
            STUDY.datasetinfo(1).comps = [];
        end;
        if ~exist('succompind')
            for ind = 1:size(STUDY.setind,2)
                succompind{ind} = STUDY.datasetinfo(STUDY.setind(1,ind)).comps; 
            end;
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
    nodip = 0;
    if update_flag % dipole information is used to select components
        
        % remove previous clusters
        % ------------------------
        STUDY.cluster = [];
        STUDY = std_checkset(STUDY, ALLEEG);
        
        % find dipoles of interest
        % ------------------------
        for si = 1:size(STUDY.setind,2)% scan datasets that are part of STUDY
            idat = STUDY.datasetinfo(STUDY.setind(1,si)).index;
            if isfield(ALLEEG(idat).dipfit, 'model')
                fprintf('Selecting dipole with less than %2.1f residual variance in dataset ''%s''\n', 100*rv, ALLEEG(idat).setname)
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
                fprintf('No dipole information found in ''%s'' dataset, using all components\n', ALLEEG.setname)
                nodip = 1;
            end
        end;
        if nodip
            error('Some dipole information is missing; thus, dipole information may not be used for clustering\n');
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
                [STUDY] = std_createclust(STUDY, ALLEEG, 'Notclust');
                STUDY.cluster(end).parent{end} = STUDY.cluster(1).name; 
                STUDY.cluster(1).child{end+1} = STUDY.cluster(end).name;
            else
                a =  ['% A subset of components from ' STUDY.cluster(cluster_ind).name  'with dipole information, were readied for clustering.' ...
                    'Other components were placed in a separate cluster'];    
                [STUDY] = std_createclust(STUDY, ALLEEG, [ 'Notclust ' num2str(cluster_ind) ] );
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
            params = [];
            for si = 1:size(STUDY.setind,2) % scan datasets that are part of STUDY
                for cond = 1:Ncond
                    idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;
                    
                    % read spectrum paramters from file
                    % ---------------------------------
                    filename = fullfile( ALLEEG(idat).filepath,[ ALLEEG(idat).filename(1:end-3) 'icaersp']);
                    if exist(filename)
                        tmp = load('-mat', filename, 'parameters');
                        params = struct( tmp.parameters{:} );
                        if  sum(params.cycles ~= cycles) | (padratio ~= params.padratio) | (alpha~= params.alpha) 
                            set_yes = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_no''), ''value'', 0);'];
                            set_no  = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_yes''), ''value'', 0);' ];
                            ersp_ans = inputgui('geometry', {[1] [1] [1] [1] [1 1] [1]}, 'uilist', ...
                                                { {'style' 'text' 'string' [upper(strcom) ' info exists for dataset: ' ...
                                                ALLEEG(idat).filename '. It does not fit requested values.' ] } ...
                                                {'style' 'text' 'string' ['Existing values are: wavelet cycles - [' num2str(params.cycles) ...
                                                '], padratio - ' num2str(params.padratio) ', and bootstrap significance - ' num2str(params.alpha) ] } {} ...
                                                {'style' 'text' 'string' ['Would you like to recalculate ' upper(strcom) ' and overwrite those values?' ]} ...
                                                {'style' 'checkbox' 'tag' 'ersp_yes' 'string' 'Yes' 'value' 1 'Callback' set_yes }  ...
                                                {'style' 'checkbox' 'tag' 'ersp_no' 'string' 'Use existing ERSP info' 'value' 0 'Callback' set_no } {} },...
                                                'helpcom', '', 'title', ['Recalculate ' upper(strcom) ' parameters -- part of std_ersp()']); 
                            if find(celltomat(ersp_ans))  == 2 % use existing ERSP info from this dataset
                                cycles   = params.cycles;
                                alpha    = params.alpha;
                                padratio = params.padratio;
                                else % Over write data in dataset
                                    delete(filename);
                            end
                        else
                            disp('Using existing ERSP information...');
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
                             [X, t] = std_erp(ALLEEG(idat), succompind{si}, timewindow);
                         else
                             [X, t] = std_erp(ALLEEG(idat), succompind{si});
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
                     X = std_topo(ALLEEG(idat), succompind{si});
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
                     X = std_topo(ALLEEG(idat), succompind{si}, 'laplacian'); 
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
                     X = std_topo(ALLEEG(idat), succompind{si}, 'gradient'); 
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
                         [X, f,overwrite] = std_spec(ALLEEG(idat),succompind{si}, ...
                                                     freqrange, fun_arg,overwrite);
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
              try,
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
              catch,
                  error('Some dipole information is missing');
              end;
              
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
                        std_ersp(ALLEEG(idat),succompind{si}, freqrange, timewindow, ...
                                 cycles, padratio, alpha, type);
                    end
                    STUDY.preclust.erspclustfreqs = freqrange;
                    STUDY.preclust.erspclusttimes = timewindow;
                    
                    % prepare ERSP / ITC data for clustering (select requested components, 
                    % mask and change to a common base if multiple conditions).
                    % --------------------------------------------------------------------
                    for k = 1:length(succompind{si})
                        if strcmpi(type, 'ersp')
                            tmp = std_readersp(ALLEEG, idattot, succompind{si}(k), timewindow, freqrange); 
                        else
                            tmp = std_readitc( ALLEEG, idattot, succompind{si}(k), timewindow, freqrange); 
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
                clustdatatmp = runpca( double(data.'), npca, 1);
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
