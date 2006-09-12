% std_precomp() - Precompute measures (ERP, spectrum, ERSP, ITC) for channels in a study. 
%                 If channels are interpolated before computing the measures, the updated 
%                 EEG datasets are also saved to disk. Called by pop_precomp(). Follow with 
%                 pop_plotstudy(). See Example below.
% Usage:    
% >> [ALLEEG,STUDY] = std_precomp(ALLEEG, STUDY, chanlist, 'key', 'val', ...);
%
% Required inputs:
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   chanlist     - [cell array] Channel names for which to precompute the
%                  selected measures. Note that the name of the channel is
%                  not case-sensitive.
% Optional inputs:
%  'erp'     - ['on'|'off'] pre-compute ERPs for each dataset.
%  'spec'    - ['on'|'off'] pre-compute spectrum for each dataset.
%              Use 'specparams' to set spectrum parameters.
%  'ersp'    - ['on'|'off'] pre-compute ERSP for each dataset.
%              Use 'erspparams' to set time/frequency parameters.
%  'itc'     - ['on'|'off'] pre-compute ITC for each dataset.
%              Use 'erspparams' to set time/frequency parameters.
%  'specparams' - [cell array] Parameters for the spectopo function are given as 
%              optional arguments. Note that it is advised to compute spectrum 
%              over all frequencies since plotting function can always reduce
%              the range of plotted frequencies.
%  'erspparams' - [cell array] Optional arguments are 'cycles', 'freqrange',
%              'padratio', 'winsize', 'alpha' (see newtimef()). Note that it 
%              is adivised to select the largest frequency range and time window
%              as plotting function are capable of plotting subranges of
%              these. An important optional parameter that is
%                    'savetrials' = ['on'|'off'] save single-trials ERSP.
%                                   Requires a lot of disk space (dataset
%                                   space on disk times 10) but allow for
%                                   refined single-trial statistics.
%
% Outputs:
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures, modified by adding preprocessing 
%                  data as pointers to Matlab files that hold the pre-clustering component measures.
%   STUDY        - the input STUDY set with pre-clustering data added, for use by pop_clust() 
%
% Example:
%   >> [ALLEEG STUDY] = std_precomp(ALLEEG, STUDY, { 'cz' 'oz' }, 'interpolate', 'on', 'erp', 'on', ...
%          'spec', 'on', 'ersp', 'on', 'erspparams', { 'cycles' [ 3 0.5 ], 'alpha', 0.01, 'padratio' 1 });
%                          
%           % This prepares, channels 'cz' and 'oz' in the STUDY datasets.
%           % If a data channel is missing in one dataset, it will be
%           % interpolated (see eeg_interp()). The ERP, spectrum, ERSP, and 
%           % ITC for each dataset is then computed. 
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2006-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2006, arno@sccn.ucsd.edu
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

function [ STUDY, ALLEEG ] = std_precomp(STUDY, ALLEEG, chanlist, varargin)
    
    if nargin < 2
        help std_precomp;
        return;
    end;
    if nargin == 2
        chanlist = []; % default to clustering the whole STUDY 
    end    
    Ncond = length(STUDY.condition);
    if Ncond == 0
        Ncond = 1;
    end

    g = finputcheck(varargin, { 'erp'         'string'  { 'on' 'off' }     'off';
                                'interpolate' 'string'  { 'on' 'off' }     'off';
                                'ersp'        'string'  { 'on' 'off' }     'off';
                                'spec'        'string'  { 'on' 'off' }     'off';
                                'itc'         'string'  { 'on' 'off' }     'off';
                                'specparams'        'cell'    {}                 {};
                                'erspparams'        'cell'    {}                 {}}, 'std_precomp');
    if isstr(g), error(g); end;
    
    % union of all channel structures
    % -------------------------------
    if isempty(chanlist)
        alllocs = ALLEEG(STUDY.datasetinfo(1).index).chanlocs;
        alllabs = { alllocs.labels };
        for index = 2:length(STUDY.datasetinfo)
           tmplocs = ALLEEG(STUDY.datasetinfo(index).index).chanlocs;
           alllocs = eeg_mergechan(alllocs, tmplocs);
        end;
        chanlist = { alllocs.labels };
    end;
    
    % Interpolate all datasets first
    % ------------------------------
    if strcmpi(g.interpolate, 'on')
        fprintf('Interpolation of data channels\n');
        fprintf('------------------------------\n');
        [ STUDY, ALLEEG ] = std_interp(STUDY, ALLEEG, chanlist);
    end;
    
    % compute ERPs
    % ------------
    if strcmpi(g.erp, 'on')
        for index = 1:length(STUDY.datasetinfo)
            std_erp(ALLEEG(STUDY.datasetinfo(index).index), 'channels', chanlist);
        end;
    end;

    % compute spectrum
    % ----------------
    if strcmpi(g.spec, 'on')
        for index = 1:length(STUDY.datasetinfo)
            std_spec(ALLEEG(STUDY.datasetinfo(index).index), 'channels', chanlist, g.specparams{:});
        end;
    end;

    % compute ERSP
    % ------------
    if strcmpi(g.ersp, 'on') | strcmpi(g.itc, 'on')
        if strcmpi(g.ersp, 'on') & strcmpi(g.itc, 'on'), type = 'both';
        elseif strcmpi(g.ersp, 'on')                   , type = 'ersp';
        else                                             type = 'itc';
        end;
        for index = 1:length(STUDY.datasetinfo)
            std_ersp(ALLEEG(STUDY.datasetinfo(index).index), 'channels', chanlist, 'type', type, g.erspparams{:});
        end;
    end;
    return;
    
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
        alpha = NaN;
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
                    if ~isnan(STUDY.setind(cond,si))
                        idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;
                    
                        % read spectrum paramters from file
                        % ---------------------------------
                        filename = fullfile( ALLEEG(idat).filepath,[ ALLEEG(idat).filename(1:end-3) 'icaersp']);
                        if exist(filename)
                            tmp = load('-mat', filename, 'parameters');
                            params = struct( tmp.parameters{:} );
                            if  sum(params.cycles ~= cycles) | (padratio ~= params.padratio) | ...
                                ( (alpha~= params.alpha) & ~( isnan(alpha) & isnan(params.alpha)) )
                                set_yes = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_no''), ''value'', 0);'];
                                set_no  = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_yes''), ''value'', 0);' ];
                                ersp_ans = inputgui('geometry', {[1] [1] [1] [1] [0.5 1] [1]}, 'uilist', ...
                                                { {'style' 'text' 'string' [upper(strcom) ' info exists for dataset: ' ...
                                                ALLEEG(idat).filename '. It does not fit requested values. Existing values used these parameters: ' ] } ...
                                                {'style' 'text' 'string' ['wavelet cycles - [' num2str(params.cycles(1)) ' ' num2str(params.cycles(2)) ...
                                                '] instead of [' num2str(cycles(1)) ' ' num2str(cycles(2)) '], padratio - ' num2str(params.padratio) ' instead of ' num2str(padratio) ...
                                                ', and bootstrap significance - ' num2str(params.alpha) ' instead of ' num2str(alpha) ] } {} ...
                                                {'style' 'text' 'string' ['Would you like to recompute ' upper(strcom) ' and overwrite those values?' ]} ...
                                                {'style' 'checkbox' 'tag' 'ersp_yes' 'string' 'Yes, recompute' 'value' 1 'Callback' set_yes }  ...
                                                {'style' 'checkbox' 'tag' 'ersp_no' 'string' 'No, use data file parameters (and existing ERSP info)' 'value' 0 'Callback' set_no } {} },...
                                                'helpcom', '', 'title', ['Recalculate ' upper(strcom) ' parameters -- part of std_ersp()']); 
                                if isempty(ersp_ans), return; end;
                                if find(celltomat(ersp_ans))  == 2 % use existing ERSP info from this dataset
                                    cycles   = params.cycles;
                                    alpha    = params.alpha;
                                    padratio = params.padratio;
                                else % Over write data in dataset
                                    delete(filename);
                                end
                            else
                                fprintf('Using %s information available on disk for dataset %d...\n', upper(strcom), idat);
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
             case 'erp', 
                  for kk = 1:length(STUDY.cluster)
                      if isfield(STUDY.cluster(kk).centroid, 'erp')
                          STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'erp');
                          STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'erp_times');
                      end;
                  end;
                  if ~isempty(succompind{si})
                      for cond = 1 : Ncond
                         if ~isnan(STUDY.setind(cond,si))
                            idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                            if ALLEEG(idat).trials == 1
                               error('No epochs in dataset: ERP information has no meaning');
                            end
                            con_t = []; con_data = [];
                            [X, t] = std_erp(ALLEEG(idat), succompind{si}, timewindow);
                            STUDY.precomp.erpclusttimes = timewindow;
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
             case 'scalp' , % NB: scalp maps are identical across conditions (within session)
                 for cond = 1 : Ncond   % Find first nonNaN index
                    if ~isnan(STUDY.setind(cond,si)), break; end
                 end       
                 idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;
                 fprintf('Computing/loading interpolated scalp maps for dataset %d...\n', idat);
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
             % ------------------------------------si
             case 'scalpLaplac'
                 for cond = 1 : Ncond   % Find first nonNaN index
                    if ~isnan(STUDY.setind(cond,si)), break; end
                 end       
                 idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
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
             case 'scalpGrad'
                 for cond = 1 : Ncond   % Find first nonNaN index
                    if ~isnan(STUDY.setind(cond,si)), break; end
                 end
                 idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
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
             case 'spec', 
                 for kk = 1:length(STUDY.cluster)
                     if isfield(STUDY.cluster(kk).centroid, 'spec')
                         STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'spec');
                         STUDY.cluster(kk).centroid = rmfield(STUDY.cluster(kk).centroid, 'spec_freqs');
                     end;
                 end;
                 if si == 1
                     overwrite = 0;
                 end
                 if ~isempty(succompind{si})
                     for cond = 1 : Ncond 
                         if ~isnan(STUDY.setind(cond,si))
                            idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index; 
                            con_f = []; con_data = [];
                            [X, f,overwrite] = std_spec(ALLEEG(idat),succompind{si}, ...
                                                     freqrange, fun_arg,overwrite);
                            STUDY.precomp.specclustfreqs = freqrange;
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
              for cond = 1 : Ncond  % Scan for first nonNaN index.
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
                error('Some dipole information is missing');
              end
              
             % cluster on ica ersp / itc values
             % --------------------------------
             case  {'ersp', 'itc'}
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
                type =  strcom;            
                if ~isempty(succompind{si})
                    idattot = [];
                    for cond = 1 : Ncond 
                        if ~isnan(STUDY.setind(cond,si))
                            idat = STUDY.datasetinfo(STUDY.setind(cond,si)).index;  
                            idattot = [idattot idat];
                            % compute ERSP/ ITC, if doesn't exist.
                            std_ersp(ALLEEG(idat),succompind{si}, freqrange, timewindow, ...
                                 cycles, padratio, alpha, type);
                        end
                    end
                    STUDY.precomp.erspclustfreqs = freqrange;
                    STUDY.precomp.erspclusttimes = timewindow;
                    
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
        
    end
    
    % Compute a second PCA of the already PCA'd data if there are too many PCA dimensions.
    % ------------------------------------------------------------------------------------
    if size(clustdata,2) > secondlevpca
        fprintf('Performing second-level PCA: reducing dimension from %d to %d \n', ...
                size(clustdata,2), secondlevpca);
        clustdata = runpca( double(clustdata.'), secondlevpca, 1);
        clustdata = clustdata.';
    end
    
    STUDY.etc.precomp.precompdata = clustdata;
    STUDY.etc.precomp.precompparams = varargin;
    STUDY.etc.precomp.precompcomps = succompind;

    % The precompering level is equal to the parent cluster that the components belong to.
    if ~isempty(cluster_ind)
        STUDY.etc.precomp.clustlevel = cluster_ind; 
    else
        STUDY.etc.precomp.clustlevel = 1; % No parent cluster (cluster on all components in STUDY).
    end

return

% erspdownsample() - down samples component ERSP/ITC images if the
%        PCA operation in the clustering feature reduction step fails.
%        This is a helper function called by eeg_precomp().

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
