% std_dipoleclusters - Plots clusters of ICs as colored dipoles in MRI
%           images (side, rear, top  and oblique angles possible)
%
% std_dipoleclusters(STUDY,ALLEEG,'key1',value1, 'key2',value2, ... );
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
% 'clusters' - [vector of numbers] list of clusters to plot in same head space
% 'title'    - [string] figure title
% 'viewnum'  - [vector] list of views to plot: 1=top, 2=side, 3=rear, 4 is an oblique view; 
%                     length(viewnum) gives the number of subplots that will be produced and the 
%                     values within the vector tell the orientation and order of views
% 'rowcolplace' - [rows cols subplot] If plotting into an existing figure, specify the number of rows, 
%                columns and the subplot number to start plotting dipole panels.
% 'colors'   - [vector or matrix] if 1 x 3 vector of RGB values, this will plot all dipoles as the
%            same color. ex. [1 0 0] is red, [0 0 1] is blue, [0 1 0] is green.
%            If a matrix, should be n x 3, with the number of rows equal to the number 
%            of clusters to be plotted and the columns should be RGB values for each. 
%            If [], will plot clusters as 'jet' colorscale from the first to the last cluster
%            requested (therefore an alternate way to control dipole color is to input a specific
%            order of clusters).
%            [] will assign colors from hsv color scale.
% 'centroid' - ['only', 'add', 'off'] 'only' will plot only cluster centroids, 'add' will superimpose
%            centroids over cluster dipoles, 'off' will skip centroid plotting and only plot 
%            cluster-member dipoles.
%
% Authors: Julie Onton, SCCN/INC UCSD, June 2009

% Copyright (C) Julie Onton, SCCN/INC/UCSD, October 11, 2009
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

function std_dipoleclusters(STUDY,ALLEEG, varargin);
    
if nargin < 2
    help std_dipoleclusters;
    return;
end

% Set default values------------------------------------------------------------------------
    if strcmp(STUDY.cluster(2),'outlier 2') % don't plot outlier cluster #2
        clustvec = [3:length(STUDY.cluster)]; % plot all clusters in STUDY
    else
        clustvec = [2:length(STUDY.cluster)]; % plot all clusters in STUDY
    end
    onecolor = [];
    colvec = [];
    centroid = 'off';
    viewnum = [1:4]; % plot all views and oblique
    rowcolplace = [2 2 1]; % 2 X 2 figure starting in subplot 1
    figureon = 1; % plot on a new figure
    ttl = ''; % no title
    %---------------------------------------------------------------------------------
    for k = 3:2:nargin
        switch varargin{k-2}
         case 'clusters'
          clusters = varargin{k-1}; % redefine from all to specified clusters
         case 'title'
          ttl = varargin{k-1};
         case 'viewnum', 
          viewnum = varargin{k-1};
         case 'rowcolplace' %, mode = varargin{k-1}; % what is mode? JO
          rowcolplace = varargin{k-1};
          if length(rowcolplace) < 3
              fprintf('\nThe variable ''rowcolplace'' must contain 3 values.\n');
              return;
          end
          row = rowcolplace(1);
          col = rowcolplace(2);
          place = rowcolplace(3);          
          figureon = 0; % don't pop a new figure if plotting into existing fig
         case 'colors'
          colvec = varargin{k-1};
         case 'centroid'
          centroid = varargin{k-1};
        end
    end    
    % adjust color matrix for dipoles:---------------
    if isempty(colvec)
        cols = jet(length(clusters));% default colors
    else
        cols = colvec; % input RGB colors
    end;    

    % extract IC cluster and data path info from STUDY structure
    clear clustcps fullpaths gdcomps
    x = cell(1,length(unique({STUDY.datasetinfo.subject})));
    subjs = unique_bc({STUDY.datasetinfo.subject});
    origlist = cell(1,length(unique({STUDY.datasetinfo.subject})));
    sets = cell(1,length(unique({STUDY.datasetinfo.subject})));
    for clust = 1:length(STUDY.cluster)
        clustcps{clust} = x;
        for st = 1:size(STUDY.cluster(clust).sets,2)
            currset = STUDY.cluster(clust).sets(1,st);            
            currcomp = STUDY.cluster(clust).comps(1,st);
            subjidx = strmatch(STUDY.datasetinfo(currset).subject,subjs);
            clustcps{clust}{subjidx}(end+1) = currcomp;
            origlist{subjidx} = [origlist{subjidx} currcomp];
            [origlist{subjidx} idx] = unique_bc(origlist{subjidx});
            sets{subjidx} = currset;
        end;    
    end
    %-----------------------------------------------------------
    % extract dipole info for ALL ICs to be plotted subj by subj
    for nx = 1:length(origlist)
        dipsources = [];
        if ~isempty(origlist{nx})
            EEG = ALLEEG(sets{nx}); % call in a dataset from subj
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end
            dipsources.posxyz = EEG.dipfit.model(origlist{nx}(1)).posxyz;
            dipsources.momxyz = EEG.dipfit.model(origlist{nx}(1)).momxyz;
            dipsources.rv = EEG.dipfit.model(origlist{nx}(1)).rv;p=1;
            
            for w = 1:length(origlist{nx})                
                dipsources(1,p).posxyz = EEG.dipfit.model(origlist{nx}(w)).posxyz;
                dipsources(1,p).momxyz = EEG.dipfit.model(origlist{nx}(w)).momxyz;
                dipsources(1,p).rv = EEG.dipfit.model(origlist{nx}(w)).rv;  
                p=p+1;
            end
            allbesa1{nx} = dipsources; new = 0;
        end
    end;    

    %-----------------------------------------------------------
    % collect cluster dipole info from extracted dipole infos (above)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    new = 1;     pp=1;  bic = 1;
    centrstr = struct('posxyz',[0 0 0],'momxyz',[0 0 0],'rv',0);
    for clst =1:length(clusters)
        clust = clusters(clst);
        centr = [];
        centr2 = [];
        for nx = 1:length(clustcps{clust})
            if ~isempty(clustcps{clust}{nx})
                for k = 1:length(clustcps{clust}{nx})  
                    if new == 1
                        allbesa = allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k)));
                        centr = [centr; allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(1,:)];
                        if size(allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz,1) > 1& allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,1) ~= 0 % actual values, not zero
                            if allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,2) > 0 % on the wrong side, switch with centr1
                                centr2 = [centr2;allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                                centr2(end,2) = centr2(end,2)*-1; centr1(end,2) = centr1(end,2)*-1;
                            else
                                centr2 = [centr2;allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                            end
                        end
                        new = 0;
                    else
                        allbesa(1,end+1) = allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k)));
                        centr = [centr; allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(1,:)];                       
                        if size(allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz,1) > 1 && allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,1) ~= 0 % actual values, not zero
                            if allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,2) > 0 % on the wrong side, switch with centr1
                                centr2 = [centr2; allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                                centr2(end,2) = centr2(end,2)*-1; centr1(end,2) = centr1(end,2)*-1;
                            else
                                centr2 = [centr2;allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                            end
                       end
                    end;                    
                    colset{pp} = cols(clst,:); pp = pp+1;                 
                end
            end
        end
        if length(allbesa) > 1
            centr = mean(centr,1);
            centr2 = mean(centr2,1);
            centrstr(clst).posxyz = centr;
            centrstr(clst).momxyz = allbesa(end).momxyz(1,:);
            centrstr(clst).rv = 2;
            centcols{clst} = cols(clst,:);
            centcols2{clst} = cols(clst,:)/5;
            if ~isempty(find(centr2))
                centrstr2(bic).posxyz = centr2;
                centrstr2(bic).momxyz = allbesa(end).momxyz(1,:);
                centrstr2(bic).rv = 2;
                bic = bic + 1; % separate count for bilaterals
            end
        end
    end
    if figureon == 1
        figure; row = 2; col = 2; place= 1;
    end
    %-------------------------------------------
    % PLOT the clusster dipoles:
    
    
    if length(allbesa) > 1
        for sbpt = 1:length(viewnum)
            if sbpt < 4
                prjimg = 'off';
            else
                prjimg = 'on';
            end
            sbplot(row,col,place)
            if strcmp(centroid,'only')
                dipplot(centrstr,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',centcols,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);      
                if ~isempty(find(centrstr2(1).posxyz)) % only if there were bilaterals
                    dipplot(centrstr2,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',centcols,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0); camzoom(.8)   
                else
                    camzoom(1)  
                end
            elseif strcmp(centroid,'add')
                dipplot(allbesa,'image','mri','gui','off','dipolelength',0,'dipolesize',25,'normlen','on','spheres','on','color',colset,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);     
                dipplot(centrstr,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',centcols2,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0); camzoom(.8)  
                if ~isempty(find(centrstr2(1).posxyz)) % only if there were bilaterals
                    dipplot(centrstr2,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',centcols2,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);camzoom(.8)   
                else
                    camzoom(1)  
                end
            else
                dipplot(allbesa,'image','mri','gui','off','dipolelength',0,'dipolesize',25,'normlen','on','spheres','on','color',colset,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);  camzoom(1.1)                   
            end
            if viewnum(sbpt) == 3
                view(0,0)
            elseif viewnum(sbpt) == 1
                view(0,90)
            elseif viewnum(sbpt) == 4
                view(63,22);
            end
            place = place+1;
        end;  
        if ~isempty(ttl)
            if sbpt == 4 % if oblique
                ph = text(-75,-75,125,ttl); set(ph,'color','r');
            elseif sbpt == 1 % 2d image:
                ph = text(-50,110,125,ttl); set(ph,'color','r');
            elseif sbpt == 2 % 2d image:
                ph = text(-75,-75,125,ttl); set(ph,'color','r');
            elseif sbpt == 3 % 2d image:
                ph = text(-100,-50,130,ttl); set(ph,'color','r');
            end
        end
    end
    
