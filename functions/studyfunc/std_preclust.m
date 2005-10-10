% eeg_preclust() - select EEGLAB dataset ICA measures (i.e. ERP, dipole, scalp maps, ...) 
%                          for clustering the ICA components in the several
%                          datasets that are part of the STUDY structure.
%
% Usage:    
%   >> [ALLEEG, STUDY] = eeg_preclust(ALLEEG, STUDY, clustind, compind, preproc1, preproc2, ...);
%
% Inputs:
%   ALLEEG       - ALLEEG data structure, can also be one EEG dataset structure.
%   STUDY        - an eeglab STUDY set (that contains some of these EEGs structures)
%   clustind      - a cluster index for further clustering- hierarchical clustering. 
%                        For example, if there is a mu cluster you wish to cluster into left mu
%                        and right mu. When left empty it is a first stage clustering (all STUDY clustering). 
%   compind    - [cell array | file name| []] component indices to cluster.
%                       Cell array of component indices for each dataset, empty brackets stands for   
%                       all the components in a set, 0 for none (dataset will not be included in pre-clustering).
%                       For instance  { [2:5] [10:15] [0] [] } will consider components 2 to 5 from the first 
%                       dataset, components 10 to 15 from the second dataset, not include dataset 3, and all the
%                       components of the fourth dataset.
%                       Or the name of a .mat file containing such a cell array variable named clustcomp.
%                       [] => use all components from all datasets.
%   preprocX   - {'comnd' 'key1' val1 'key2' val2 ...} 
%                'comnd'   ['erp'|'spec'|'scalp'|'dipoles'|'itc'|'ersp'|'scalpLaplac'|'scalpGrad'] 
%                              enter preprocessing command to apply to the data.
%                  * 'erp'     cluster on the ICA component ERPs,
%                  * 'dipoles'  cluster on the component [X Y Z] dipole locations
%                  * 'dipselect'  select components to cluster on that have residual 
%                                  dipole variance less than a specified threshold. 
%                  * 'spec'    cluster on the ICA component activity spectra 
%                                 (with the baseline removed).
%                  * 'scalp'  cluster on the ICA topoplot scalp maps (or their absolute values),
%                  * 'scalpLaplac'  cluster on the ICA topoplot scalp map Laplacians (or their absolute values),
%                  * 'scalpGrad'  cluster on the ICA topoplot scalp map Gradients (or their absolute values),
%                  * 'ersp' cluster on ICA components ERSP. (requires:
%                              'cycles', 'freqrange', 'padratio', 'timewindow', 'alpha').
%                  * 'itc' cluster on ICA components ITC.(requires:
%                              'cycles', 'freqrange', 'padratio', 'timewindow', 'alpha').
%
%                'key*'   optional inputs to accompany the 'comnd'   
%                  * 'npca'    = [integer] number of PCA components (PCA dimension) of the selected data
%                                    to retain for clustering. {default: 5}
%                  * 'norm'    =   [0|1] 1 normalizes the PCA components so the variance of first PCA 
%                                    component is 1 (good when using several preprocessing commands - 
%                                    'ersp','scalp',...). {default: 1}
%                  * 'weight'  =   [integer] weight with respect to other preprocessing commands.
%                  * 'freqrange' = [min max] frequency range in Hz for spectrum commands and 'ersp', 'itc' options.  
%                  * 'timewindow' = [min max] time window in seconds for 'erp' 'ersp' 'itc' options.  
%                  * 'abso'    =   [0|1] 1 take absolute values of topoplot map, map Gradient or map Laplacian. {default: 1}
%                  * 'cycles'  =   [0| wavecycles factor] for ERSP/ ITC (see timef() for details). {default: 0 (FFT)}
%                  * 'padratio'  =   [integer] for ERSP/ ITC (see timef() for details). {default: 1}
%                  * 'alpha'  =   [integer] bootstrap significance for ERSP/ ITC (see timef() for details). {default: 0.01}
%                  * 'funarg'  =  [cell array] optional function arguments for spectrum calculation (see spectopo() for details). {default: none}
%                  * 'rv'  =   [number < 1] for 'dipselect', a threshold on the component residual variance, 
%                               only components that have a lower rv will be clustered. {default: 0 (all components)}
%
% Outputs:
%   ALLEEG       - the input ALLEEG data structure modified with the preprocessing data 
%                      (pointers to float files that holds the ERSP, spectrum etc. for future use 
%                       by the clustering functions)
%   STUDY        - the input STUDY modified with the pre-clustering data, for use by the pop_clust() function.
%
% Example:
%   >> [ALLEEG  STUDY] = eeg_preclust(ALLEEG, STUDY, [], [] , { 'dipselect'  'rv'  0.15  } ,...
%                        { 'spec'  'npca'   10   'norm'  1   'weight'  1  'freqrange' [ 3  25  ] } , ...
%                        { 'erp'  'npca'   10   'norm'  1   'weight'  2  'timewindow' [ 350 500  ] } ,...
%                        { 'scalp'  'npca'  10  'norm'  1  'weight'  2  'abso'  1  } , { 'dipoles'   'norm'   1   'weight'  15  } , ...
%                        { 'ersp'  'npca'  10  'freqrange' [ 3  25  ] 'cycles' [ 3 0.5 ] 'alpha'  0.01  'padratio'  4  'timewindow' [ -1600  1495 ] 'norm'  1  'weight'  1  } ,...
%                        { 'itc'  'npca'  10  'freqrange' [ 3  25  ] 'cycles' [ 3 0.5 ] 'alpha'  0.01  'padratio'  4  'timewindow' [ -1600  1495 ] 'norm'  1  'weight'  1  });
%   % This is a first stage clustering done on all the components in the
%   %datasets, except for components with residual dipole variance higher
%   % than 0.15. The clustering will be based on the components spectra in
%   % the [3 25] Hz frequency range, the components 'erp' in the [350 500]ms
%   % time window, the absolute component scalp maps, the dipole location,
%   % and the ersp and itc information. 
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
    
    if nargin < 3
        help eeg_preclust;
        return;
    end;
    
    Ncond = length(STUDY.condition);
    if Ncond == 0
        Ncond = 1;
    end
    
    % Get component indices that are part of the cluster 
    if ~isempty(cluster_ind)
        if length(cluster_ind) ~= 1
            error('Only one cluster can be clustered at a time, if you wish to cluster more than one cluster merge clusters first.');
        end
        for k = 1 : size(STUDY.setind,2) % go over the sets from the first condition (if there are some)
            sind = find(STUDY.cluster(cluster_ind).sets(1,:) == k);%the component indices that belong to the dataset in the cluster 
            if isempty(sind)
                succompind{k} = 0;
            else
                succompind{k} = STUDY.cluster(cluster_ind).comps(sind);
            end
        end
    end

    if ~isempty(components_ind) % there is a component selection
        if isstr(components_ind) % if input is a .mat file load component indices
            try, 
                eval( [ 'load ' components_ind ]);
            catch
                error('eeg_preclust: compsind input is not a valid file')
            end
            if isempty('clustcomp')
                error('eeg_preclust: compsind input .mat file must have a string argument named clustcomp');
            end
            components_ind = clustcomp;
        end
        if length(components_ind) ~= size(STUDY.setind,2)
            error('eeg_preclust: there must be as many elements in the cell array of component indices compsind input as there are subjects x sessions in STUDY.');
        end;
        if ~isempty(cluster_ind) % components to cluster on must be both part of the specified cluster and slected components
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
            succompind{ind} = succompind{ind}(find(succompind{ind}));%remove zeros
            succompind{ind} = sort(succompind{ind}); %sort the components
        end
    end;
    
    % Scan which component to remove (no dipole info)
    % ------------------------------
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
    if update_flag % dipole infromation is used to select components
        for si = 1:size(STUDY.setind,2)% scan datasets that are part of STUDY
            idat = STUDY.datasetinfo(STUDY.setind(1,si)).index;
            if ~isfield(ALLEEG(idat), 'dipfit')
                error('No dipole information exists in data')
            end
            indrm = []; % components that will be removed from the clustering
            indleft = []; % components that are left in clustering
            for icomp = succompind{si} % scan components
                if (ALLEEG(idat).dipfit.model(icomp).rv >= rv) | isnan(ALLEEG(idat).dipfit.model(icomp).rv) % Components have rv bigger than asked for 
                    %fprintf('Component %d of dataset %d has no dipole info and was removed\n', icomp, idat)
                    indrm = [indrm icomp];
                else
                    indleft = [indleft icomp];
                end;
            end;
             if ~isempty(indrm)
                 succompind{si} = indleft;
             end
        end;
        % Create a cluster with removed components 
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
                a =  ['% A subset of components, which have dipole information, were taken for clustering.'...
                        'The rest of the components are placed in a separate cluster'];
                [STUDY] = cls_createclust(STUDY, ALLEEG, 'Notclust');
                STUDY.cluster(end).parent{end} = STUDY.cluster(1).name; 
                STUDY.cluster(1).child{end+1} = STUDY.cluster(end).name;
            else
                a =  ['% A subset of components from ' STUDY.cluster(cluster_ind).name  ', which have dipole information, were taken for clustering.' ...
                    'The rest of the components are placed in a separate cluster'];    
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
        
        % In case ersp / itc values already exist in some of the datasets,
        % make sure they are the same as the requested parameters.
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
                              error('No epochs in dataset, ERP information has no meaning');
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
                         else %concatenate across conditions
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
             case 'scalp' , %scalp maps are identicle across conditions
                 idat = STUDY.datasetinfo(STUDY.setind(1,si)).index; 
                 if ~isempty(succompind{si})
                     [tmp, X] = cls_scalp(ALLEEG(idat), succompind{si});
                     if ~isempty(tmp)
                         ALLEEG(idat).etc = tmp;
                     end
                     if abso %absolute values
                         data = [ data; abs(X) ];
                     else
                         data = [ data; X ];
                     end
                     clear X tmp;
                 end       
             % select ica Laplacian scalp maps
             % -----------------------------------
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
             % select ica Gradient scalp maps
             % ------------------------------
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
             % select ica spectrum
             % --------------------------
             case 'spec'   , 
                 if si == 1
                    overwrite = 0;
                end
                if ~isempty(succompind{si})
                    for cond = 1 : Ncond 
                         idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                         if ~isempty('freqrange')
                             [tmp, X, f,overwrite] = cls_spec(ALLEEG(idat),succompind{si}, freqrange, fun_arg,overwrite);
                             if ~isempty(tmp)
                                 ALLEEG(idat).etc = tmp;
                             end
                         else
                             [tmp, X, f,overwrite] = cls_spec(ALLEEG(idat),succompind{si}, [], fun_arg,overwrite);
                             if ~isempty(tmp)
                                 ALLEEG(idat).etc = tmp;
                             end
                         end
                         if cond == 1
                             con_data = X;
                             con_f = f;
                         else %concatenate across conditions
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
             case 'dipoles' %dipoles are identicle across conditions (no need to scan across condition)
              idat = STUDY.datasetinfo(STUDY.setind(1,si)).index; 
              count = size(data,1)+1;
              for icomp = succompind{si}
                  % select between 3 suboptions
                  % ---------------------------
                  if ~isempty(succompind{si})
                      ldip = 1;
                      if size(ALLEEG(idat).dipfit.model(icomp).posxyz,1) == 2 % two dipoles model
                          if any(ALLEEG(idat).dipfit.model(icomp).posxyz(1,:)) & any(ALLEEG(idat).dipfit.model(icomp).posxyz(2,:)) %both dipoles exist
                              [garb ldip] = max(ALLEEG(idat).dipfit.model(icomp).posxyz(:,2)); %find the left most dipole
                          elseif any(ALLEEG(idat).dipfit.model(icomp).posxyz(2,:)) 
                              ldip = 2; %left most dipole is the only one that exists
                          end
                      end
                      data(count,:) = ALLEEG(idat).dipfit.model(icomp).posxyz(ldip,:);
                      count = count+1;
                  end
              end                     
             % cluster on ica ersp / itc values
             % --------------------------
             case  {'ersp', 'itc'}
                type =  strcom;            
                if ~isempty(succompind{si})
                    idattot = [];
                    for cond = 1 : Ncond 
                         idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                         idattot = [idattot idat];
                         % compute ERSP/ ITC, if doesn't exist.
                         [tmp] = cls_ersp(ALLEEG(idat),succompind{si}, freqrange, timewindow, cycles, padratio, alpha, type);
                         if ~isempty(tmp)
                             ALLEEG(idat).etc = tmp;
                         end
                    end
                    % prepare ERSP / ITC data for clustering (selecet requested components, mask 
                    % and change ro a common base if multiple conditions).
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
                              disp('Warning: wrong number of components (rows) generated by string command');
                          end;
                      end;
                  end
            end;

            
        end; % end scan datasets
        
        % adjust number of PCA values
        % ---------------------------
        if isnan(npca), npca = 5; end; % default number of components
        if npca >= size(data,2)
            % no need to run PCA, just copy the data
            % but we will still run it to "normalize" coordinates
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
            fprintf('Retaining the three dimensional location of the dipoles (normalize:%s; weight:%d)\n', ...
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
                        [data, freqs, times] = erspdownsample(data,4, freqs,times,Ncond); %down sample frequency by 2 and times by 2
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
                %normalize each cordinate by the std of the radius
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
    
    % Compute a second PCA of the already PCA data if there are too many
    % PCA dimentions.
    if size(clustdata,2) > 10
        fprintf('Performing a second level of PCA. The dimension is reduced from %d to 10 \n', size(clustdata,2));
        clustdata = runpca( clustdata.', 10, 1);
        clustdata = clustdata.';
    end
    
    STUDY.etc.preclust.preclustdata = clustdata;
    STUDY.etc.preclust.preclustparams = varargin;
    STUDY.etc.preclust.preclustcomps = succompind;
    % The preclustering level equals to the parent cluster that the components belong to.
    if ~isempty(cluster_ind)
        STUDY.etc.preclust.clustlevel = cluster_ind; 
    else
        STUDY.etc.preclust.clustlevel = 1; % No parent cluster (cluster on all components in STUDY).
    end

return
