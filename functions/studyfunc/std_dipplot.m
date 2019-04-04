% std_dipplot() - Commandline function to plot cluster component dipoles. Dipoles for each
%                 named cluster is displayed in a separate figure. To view all the clustered 
%                 components in the STUDY on the same figure (in a separate subplot), all 
%                 STUDY clusters must be requested.
%                 To visualize dipoles, they first must be stored in the EEG dataset structures
%                 using dipfit(). Only components that have dipole locations will be displayed,
%                 along with the cluster mean dipole (in red). 
% Usage:    
%                 >> [STUDY] = std_dipplot(STUDY, ALLEEG, clusters);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in 
%                the STUDY. ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [numeric vector | 'all']  -> specific cluster numbers to plot.
%                'all'  -> plot all clusters in STUDY.
%                {default: 'all'}.
%   'comps'    - [numeric vector]  -> indices of the cluster components to plot.
%                'all'  -> plot all the components in the cluster 
%                {default: 'all'}.
%   'mode'     - ['together'|'apart'|'multicolor'] Display all requested cluster on one 
%                figure ('together') or separate figures ('apart'). 
%                'together'-> plot all 'clusters' individuall in one multi-pane figure (without the gui).
%                'apart'   -> plot each cluster in a separate figure. 
%                'multicolor' -> plot all clusters in one figure, 
%                Note that this parameter has no effect if the 'comps' option (above) is used.
%                {default: 'together'}
%   'figure'   - ['on'|'off'] plots on a new figure ('on')  or plots on current
%                figure ('off'). If 'figure','off' does not display gui controls,
%                Useful for incomporating a cluster dipplot into a complex figure. 
%                {default: 'on'}. 
%   'groups'   - ['on'|'off'] use different colors for different groups.
%                {default: 'off'}.
%   'dipcolor' - [cell vector] color for dipoles in each cluster. (multicolor mode)
%   'dipsize'  - [numeric vector] size for each cluster. (multicolor mode)
%                {default if unspecified: will automatically color/size each
%                cluster}
% Outputs:
%   STUDY      - the input STUDY set structure modified with plotted cluster 
%                mean dipole, to allow quick replotting (unless cluster means 
%                already exists in the STUDY).  
%   Example:
%   >> [STUDY] = std_dipplot(STUDY,ALLEEG, 'clusters', 5, 'mode', 'apart', 'figure', 'off');
%                % Plot cluster-5 component dipoles (in blue), plus ther mean dipole (in red), 
%                % on an exisiting (gui-less) figure. 
%
%  See also  pop_clustedit(), dipplot()        
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005
%          'groups' added by Makoto Miyakoshi on June 2012.
%          'multicolor' mode added by John Iversen to draw all clusters on
%           a single panel, with each cluster indicated by different color/size dipoles.

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 08, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_dipplot(STUDY, ALLEEG, varargin)

% Set default values
cls = []; % plot all clusters in STUDY
figureon = 1; % plot on a new figure
mode = 'apart';

STUDY = pop_dipparams(STUDY, 'default');

opt_dipplot = {'projlines',STUDY.etc.dipparams.projlines, 'axistight', STUDY.etc.dipparams.axistight, 'projimg', STUDY.etc.dipparams.projimg, 'spheres', 'on', 'dipolelength', 0, 'density', STUDY.etc.dipparams.density};

dipcolor = [];
dipsize = [];
%, 'spheres', 'on'
groupval = 'off';
for k = 3:2:nargin
    switch varargin{k-2}
        case 'clusters'
            if isnumeric(varargin{k-1})
                cls = varargin{k-1};
                if isempty(cls)
                    cls = 2:length(STUDY.cluster);
                end
            else
                if ischar(varargin{k-1}) && strcmpi(varargin{k-1}, 'all')
                    cls = 2:length(STUDY.cluster);
                else
                    error('std_dipplot: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
            if length(cls) == 1, mode = 'apart'; else mode = 'together'; end
        case 'comps'
            if strcmpi(STUDY.etc.dipparams.density, 'on')
                disp('Single dipole should not be plotted using dipole density, reverting to dipole plotting');
                opt_dipplot{end} = 'off';
            end
            STUDY = std_plotcompdip(STUDY, ALLEEG,  cls, varargin{k-1}, opt_dipplot{:});
            return;
        case 'plotsubjects', % do nothing
        case 'mode', mode = varargin{k-1};
        case 'groups', groupval = varargin{k-1};
        case 'figure'
            if strcmpi(varargin{k-1},'off') 
                if strcmpi(STUDY.etc.dipparams.density, 'on')
                    disp('Cannot plot dipole density within figure, reverting to dipole plotting');
                    opt_dipplot{end} = 'off';
                end
                opt_dipplot{end + 1} = 'gui';
                opt_dipplot{end + 1} = 'off';
                figureon = 0;
            end
      case 'dipcolor'
        dipcolor = varargin{k-1};
      case 'dipsize'
        dipsize = varargin{k-1};
    end
end

if strcmpi(mode, 'together')
    if strcmpi(STUDY.etc.dipparams.density, 'on')
        disp('Cannot plot dipole density within figure, reverting to dipole plotting');
        opt_dipplot{end} = 'off';
    end
end

% select clusters to plot
% -----------------------
if isempty(cls)
    tmp =[];
    cls = 2:length(STUDY.cluster); % plot all clusters in STUDY
    for k = 1: length(cls)
        % don't include 'Notclust' clusters
        if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) && ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)
            tmp = [tmp cls(k)];
        end
    end
    cls = tmp;
end

if strcmpi(mode, 'apart')  % case each cluster on a separate figure
    for clus = 1: length(cls) % For each cluster requested
        if length(STUDY.cluster(cls(clus)).comps) > 0  % check there are comps in cluster
            max_r = 0;
            clear cluster_dip_models;
            len = length(STUDY.cluster(cls(clus)).comps);
            ndip = 0;
            dip_ind = [];
            if ~isfield(STUDY.cluster(cls(clus)),'dipole')
                STUDY = std_centroid(STUDY,ALLEEG, cls(clus) , 'dipole');
            elseif isempty(STUDY.cluster(cls(clus)).dipole)
                STUDY = std_centroid(STUDY,ALLEEG, cls(clus) , 'dipole');
            end
            for k = 1:len
                abset   = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).index;
                subject = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).subject;
               if ~isfield(ALLEEG(abset), 'dipfit')
                   warndlg2(['No dipole information available in dataset ' ALLEEG(abset).filename ' , abort plotting'], 'Aborting plot dipoles');
                   return;
               end
               comp = STUDY.cluster(cls(clus)).comps(k); 
               cluster_dip_models(k).posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
               cluster_dip_models(k).momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
               cluster_dip_models(k).rv = ALLEEG(abset).dipfit.model(comp).rv;
               if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                   if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                       load('-mat', ALLEEG(abset).dipfit.hdmfile);
                       max_r = max(max_r, max(vol.r));
                   else % old version of dipfit
                       max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
                   end
               end
               comp_to_disp{k} = [subject  ', ' 'IC' num2str(comp) ];
               if ~isempty(cluster_dip_models(k).posxyz)
                   ndip = ndip +1;
                   dip_ind = [dip_ind k];
               end
            end % finished going over cluster comps
            
            STUDY.cluster(cls(clus)).dipole = computecentroid(cluster_dip_models);
            cluster_dip_models(end + 1) = STUDY.cluster(cls(clus)).dipole;
           
           % additional options
           % ------------------
           dip_color = cell(1,ndip+1);
           dip_color(1:ndip) = {'b'};
           dip_color(end) = {'r'};
           options = opt_dipplot;
           options{end+1} =  'mri';
           options{end+1} =  ALLEEG(abset).dipfit.mrifile;
           options{end+1} =  'coordformat';
           options{end+1} =  ALLEEG(abset).dipfit.coordformat;
           options{end+1} =  'dipnames';
           options{end+1} = {comp_to_disp{dip_ind } [STUDY.cluster(cls(clus)).name ' mean']};
           options{end+1} = 'color';
           options{end+1} = dip_color;
           
           % if 'groups'==1, overwrite cluster_dip_models, dip_color and dipnames in option -makoto
           if strcmpi(groupval, 'on')
               [cluster_dip_models, options] = dipgroups(ALLEEG, STUDY, cls, comp_to_disp, cluster_dip_models, options);
               break
           end
  
           if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
               options{end+1} = 'sphere';
               options{end+1} = max_r;
           else
               options{end+1} = 'meshdata';
               options{end+1} = ALLEEG(abset).dipfit.hdmfile;
           end
           %enable both lines and images to be projected; increase #dipole limits
           if ndip < 20 && strcmpi(options{1}, 'projlines') && length(cls) == 1
               options{2} = 'on';
           end
           if ndip < 20 && strcmpi(options{5}, 'projimg') && length(cls) == 1
             options{6} = 'on';
           end
           
           % Dealing with projection lines
           if strcmpi(options{2},'on')
               projlinesvect = ones(1,length(cluster_dip_models));
           elseif strcmpi(options{2},'off') && strcmp(STUDY.etc.dipparams.centrline,'on')
               projlinesvect = zeros(1,length(cluster_dip_models));
               projlinesvect(end) = 1;
           elseif strcmpi(options{2},'off')
               projlinesvect = zeros(1,length(cluster_dip_models));
           end
           options{2} = projlinesvect;
           
           if figureon
               dipplot(cluster_dip_models, options{:});
               fig_h = gcf;
               set(fig_h,'Name', [STUDY.cluster(cls(clus)).name ' - ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) ...
                       ' sets - ' num2str(length(STUDY.cluster(cls(clus)).comps)) ' components (' num2str(ndip) ' dipoles)' ],'NumberTitle','off');
           else
               dipplot(cluster_dip_models, options{:},'view', [0.5 -0.5 0.5]);
               for gind = 1:length(options) % remove the 'gui' 'off' option
                   if ischar(options{gind}) 
                       if strfind(options{gind}, 'gui')
                           break;
                       end
                   end
               end
               options(gind:gind+1) = [];
               dipinfo.dipmod =  cluster_dip_models;
               dipinfo.op = options;
               diptitle = [STUDY.cluster(cls(clus)).name ', ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) ' sets -' ...
                   num2str(length(STUDY.cluster(cls(clus)).comps)) ' components (' num2str(ndip) ' dipoles)' ];
               dipinfo.title = diptitle;
               set(gcf, 'UserData', dipinfo);
               set(gca,'UserData', dipinfo);
               rotate3d off;
               axcopy(gca, ['dipinfo = get(gca, ''''UserData''''); dipplot(dipinfo.dipmod, dipinfo.op{:}); set(gcf, ''''Name'''', dipinfo.title,''''NumberTitle'''',''''off''''); ']);
           end
        end % finished the if condition that cluster isn't empty
    end % finished going over requested clusters
end 

if strcmpi(mode, 'together')  % case all clusters are plotted in the same figure (must be a new figure)
    N = length(cls);
    rowcols(2) = ceil(sqrt(N)); % Number of rows in the subplot figure.
    rowcols(1) = ceil(N/rowcols(2));
    fig_h = figure; 
    orient tall
    set(fig_h,'Color', 'black');
    set(fig_h,'Name', 'All clusters dipoles','NumberTitle','off');
    set(fig_h, 'resize','off');
    for l = 1:N
        len = length(STUDY.cluster(cls(l)).comps);
        max_r = 0;
        clear cluster_dip_models;
        if ~isfield(STUDY.cluster(cls(l)),'dipole')
            STUDY = std_centroid(STUDY,ALLEEG, cls(l), 'dipole');
        elseif isempty(STUDY.cluster(cls(l)).dipole)
            STUDY = std_centroid(STUDY,ALLEEG, cls(l), 'dipole');
        end
        for k = 1: len
            abset = STUDY.datasetinfo(STUDY.cluster(cls(l)).sets(1,k)).index;
           if ~isfield(ALLEEG(abset), 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(abset) ' , abort plotting'], 'Aborting plot dipoles');
               return;
           end
           comp = STUDY.cluster(cls(l)).comps(k);
           cluster_dip_models(k).posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
           cluster_dip_models(k).momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
           cluster_dip_models(k).rv = ALLEEG(abset).dipfit.model(comp).rv;
           if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                   load('-mat', ALLEEG(abset).dipfit.hdmfile);
                   max_r = max(max_r, max(vol.r));
               else % old version of dipfit
                   max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
               end
           end
        end % finished going over cluster comps
        STUDY.cluster(cls(l)).dipole = computecentroid(cluster_dip_models);
        cluster_dip_models(end + 1) = STUDY.cluster(cls(l)).dipole;
        dip_color = cell(1,length(cluster_dip_models));
        dip_color(1:end-1) = {'b'};
        dip_color(end) = {'r'};
        options = opt_dipplot;
        options{end + 1} =  'gui';
        options{end + 1} =  'off';
        options{end+1} =  'mri';
        options{end+1} =  ALLEEG(abset).dipfit.mrifile;
        options{end+1} =  'coordformat';
        options{end+1} =  ALLEEG(abset).dipfit.coordformat;
        options{end+1} = 'color';
        options{end+1} = dip_color;
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
            options{end+1} = 'sphere';
            options{end+1} = max_r;
        else
            options{end+1} = 'meshdata';
            options{end+1} = ALLEEG(abset).dipfit.hdmfile;
        end
        
        % Dealing with projection lines
        if strcmpi(options{2},'on')
            projlinesvect = ones(1,length(cluster_dip_models));
        elseif strcmpi(options{2},'off') && strcmp(STUDY.etc.dipparams.centrline,'on')
            projlinesvect = zeros(1,length(cluster_dip_models));
            projlinesvect(end) = 1;
        elseif strcmpi(options{2},'off')
            projlinesvect = zeros(1,length(cluster_dip_models));
        end
           options{2} = projlinesvect;
        
        hsbplot = subplot(rowcols(1),rowcols(2),l);
        % Creating new axis for title only
        set(hsbplot,'visible','off');
        htitle = title([ STUDY.cluster(cls(l)).name ' (' num2str(length(unique(STUDY.cluster(cls(l)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(l)).comps)),' ICs)'],'color','white','Visible', 'on');
        if rowcols(1)> 1, set(htitle,'Position', get(htitle,'Position')+ [0 0.05 0].*get(htitle,'Position') ); end
        axes('Position', get(hsbplot,'Position'),'Color', 'black');
        dipplot(cluster_dip_models, options{:});

        %diptitle = [STUDY.cluster(cls(l)).name ', ' num2str(length(unique(STUDY.cluster(cls(l)).sets(1,:)))) 'Ss'];
        %title(diptitle, 'Color', 'white');
        % Complex axcopy
        %if l == 1
        %    for gind = 1:length(options) % remove the 'gui' 'off' option
        %        if ischar(options{gind}) 
        %            if strfind(options{gind}, 'gui')
        %                break;
        %            end
        %        end
        %    end
        %    options(gind:gind+1) = [];
        %end
        %dipinfo.dipmod =  cluster_dip_models;
        %dipinfo.op = options;
        %dipinfo.title = diptitle;
        %set(gcf, 'UserData', dipinfo);
        %set(gca,'UserData', dipinfo);
        %axcopy(gcf, ['dipinfo = get(gca, ''''UserData''''); dipplot(dipinfo.dipmod, dipinfo.op{:}); set(gcf, ''''Name'''', dipinfo.title,''''NumberTitle'''',''''off'''');']);
   end %finished going over all clusters
   set(fig_h, 'resize','on');
end % finished case of 'all' clusters

% ========================================================================================
% multicolor mode
%   all clusters are plotted in the same axis, with each cluster indicated by color/size
%   also enable 'data cursor' so that in data tip mode, clicking on a dipole
%    will display its name
if strcmpi(mode, 'multicolor')
  N = length(cls);
  %%%%%%%%%%%%%%%%%%%%% color list %%%%%%%%%%%%%%%%%%%%%
  % This color list was developped for std_envtopo
  % modified from dipgroups below
  colors{1}  = [1 1 1];            % White
  colors{2}  = [1 1 0];            % Yellow
  colors{3}  = [1 0 1];            % Fuchsia
  colors{4}  = [1 0 0];            % Red
  colors{5}  = [0.875 0.875 0.875]; % Silver
  colors{6}  = [0.5 0.5 0.5];      % Gray
  colors{7}  = [0.5 0.5 0];        % Olive
  colors{8}  = [0.5 0 0.5];        % Purple
  colors{9}  = [0.5 0 0];          % Maroon
  colors{10} = [0 1 1];            % Aqua
  colors{11} = [0 1 0];            % Lime
  colors{12} = [0 0.5 0.5];        % Teal
  colors{13} = [0 0.5 0];          % Green
  colors{14} = [0 0 1];            % Blue
  colors{15} = [0 0 0.5];          % Navy
  colors{16} = [0 0 0];            % Black
  % Choosing and sorting 13 colors for clusters: Red, Green, Blue,
  % Fuchsia, Lime, Aqua, Maroon, Olive, Purple, Teal, Navy, Gray, and White
  colors = colors([4 13 14 3 11 10 9 7 8 12 15 6 1 ]);
  fig_h = figure;
  orient tall
  set(fig_h,'Color', 'black');
  set(fig_h,'Name', 'All clusters dipoles','NumberTitle','off');
  set(fig_h, 'resize','off');
  
  idx = 0; %cumulative dipole index
  centroidIdx = [];
  clear cluster_dip_models;
  for l = 1:N %loop over clusters
    len = length(STUDY.cluster(cls(l)).comps);
    max_r = 0;
    
    %color, size for every cluster can be passed as arguments, or
    %automatically iterate through a set of colors/sizes
    if ~isempty(dipcolor)
      clusterColors{l} = dipcolor{idx};
    else
      colorIndex = mod(l-1, length(colors))+1;
      clusterColors{l} = colors{colorIndex};
    end
    
    if ~isempty(dipsize)
      clusterSizes(l) = dipsize(idx);
    else
      %after rotating through color list once, change size of dipole
      if l <= length(colors)
        clusterSizes(l) = 30;
      elseif l <= 2*length(colors)
        clusterSizes(l) = 20;
      else
        clusterSizes(l) = 15;
      end
    end
    
    if ~isfield(STUDY.cluster(cls(l)),'dipole')
      STUDY = std_centroid(STUDY,ALLEEG, cls(l), 'dipole');
    elseif isempty(STUDY.cluster(cls(l)).dipole)
      STUDY = std_centroid(STUDY,ALLEEG, cls(l), 'dipole');
    end
    clustStartIdx = idx + 1;
    
    for k = 1:len %loop over components within a cluster
      idx = idx + 1;
      abset = STUDY.datasetinfo(STUDY.cluster(cls(l)).sets(1,k)).index;
      subjname = STUDY.datasetinfo(STUDY.cluster(cls(l)).sets(1,k)).subject;
      if ~isfield(ALLEEG(abset), 'dipfit')
        warndlg2(['No dipole information available in dataset ' num2str(abset) ' , abort plotting'], 'Aborting plot dipoles');
        return;
      end
      comp = STUDY.cluster(cls(l)).comps(k);
      cluster_dip_models(idx).posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
      cluster_dip_models(idx).momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
      cluster_dip_models(idx).rv = ALLEEG(abset).dipfit.model(comp).rv;
      if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
        if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
          load('-mat', ALLEEG(abset).dipfit.hdmfile);
          max_r = max(max_r, max(vol.r));
        else % old version of dipfit
          max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
        end
      end
      
      dip_color{idx} = clusterColors{l};
      dip_size(idx) = clusterSizes(l);
      dip_label{idx} = sprintf('Cls %d (%s IC%d)',cls(l), subjname, comp);
    end % finished going over cluster comps
    
    %add the cluster centroid
    clustEndIdx = idx;
    STUDY.cluster(cls(l)).dipole = computecentroid(cluster_dip_models(clustStartIdx:clustEndIdx));
    idx = idx + 1;
    centroidIdx(end+1) = idx;
    cluster_dip_models(idx) = STUDY.cluster(cls(l)).dipole;
    dip_color(idx) = {'k'};
    dip_size(idx) = 10;
    dip_label{idx} = sprintf('Cls %d centroid',cls(l));
    clusterLabels{l} = sprintf('Cls %d',cls(l));
    
  end %loop over clusters
  
  options = opt_dipplot;
  options{end + 1} =  'gui';
  options{end + 1} =  'off';
  options{end+1} =  'mri';
  options{end+1} =  ALLEEG(abset).dipfit.mrifile;
  options{end+1} =  'coordformat';
  options{end+1} =  ALLEEG(abset).dipfit.coordformat;
  options{end+1} = 'color';
  options{end+1} = dip_color;
  options{end+1} = 'dipolesize';
  options{end+1} = dip_size;
  options{end+1} = 'dipnames';
  options{end+1} = dip_label;
  if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
    options{end+1} = 'sphere';
    options{end+1} = max_r;
  else
    options{end+1} = 'meshdata';
    options{end+1} = ALLEEG(abset).dipfit.hdmfile;
  end
  
  % Dealing with projection lines
  if strcmpi(options{2},'on')
      projlinesvect = ones(1,length(cluster_dip_models));
  elseif strcmpi(options{2},'off') && strcmp(STUDY.etc.dipparams.centrline,'on')
      projlinesvect = zeros(1,length(cluster_dip_models));
      projlinesvect(centroidIdx) = 1;
  elseif strcmpi(options{2},'off')
      projlinesvect = zeros(1,length(cluster_dip_models));
  end
  options{2} = projlinesvect;
  
  dipplot(cluster_dip_models, options{:});
  set(fig_h, 'resize','on');
  
  %data cursor will show component label
  datacursormode on
  dcm = datacursormode(gcf);
  set(dcm,'Enable','on', 'UpdateFcn', @componentText)
  
end % multicolor. Supporting functions at end of file

% ========================================================================================

% std_plotcompdip() - Commandline function, to visualizing cluster components dipoles. 
%                   Displays the dipoles of specified cluster components with the cluster mean 
%                   dipole on separate figures. 
%                   To visualize dipoles they first must be stored in the EEG dataset structures
%                   using dipfit(). Only components that have a dipole locations will be displayed,
%                   along with the cluster mean dipole in red. 
% Usage:    
%                   >> [STUDY] = std_plotcompdip(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster {default: 'all'}. 
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster
%                     dipole mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= 1;  
%                         >> [STUDY] = std_plotcompdip(STUDY,ALLEEG, cluster, comps);
%                    Plots component 1 dipole in blue with the cluster 4 mean dipole in red. 
%
%  See also  pop_clustedit, dipfit, std_dipplot         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 08, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_plotcompdip(STUDY, ALLEEG, cls, comp_ind, varargin)
if ~exist('cls')
    error('std_plotcompdip: you must provide a cluster number as an input.');
end
if isempty(cls)
   error('std_plotcompdip: you must provide a cluster number as an input.');
end
if nargin == 3 % no components indices were given
    % Default plot all components of the cluster
    [STUDY] = std_dipplot(STUDY, ALLEEG, 'clusters', cls);
    return
end

for ci = 1:length(comp_ind)
    abset = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).index;
    comp = STUDY.cluster(cls).comps(comp_ind(ci));
    subject = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).subject;
    if ~isfield(ALLEEG(abset), 'dipfit')
        warndlg2(['No dipole information available in dataset ' num2str(abset) ' , abort plotting'], 'Aborting plot dipoles');
        return;
    end
    if length(comp_ind) == 1 && isempty(ALLEEG(abset).dipfit.model(comp).posxyz)
        warndlg2(strvcat('There is no dipole information available in', ...
                       [ 'dataset ' num2str(abset) ' for this component, abort plotting']), 'Aborting plot dipoles');
        return;
    end
    if ~isfield(STUDY.cluster(cls),'dipole')
        STUDY = std_centroid(STUDY,ALLEEG, cls , 'dipole');
    elseif isempty(STUDY.cluster(cls).dipole)
        STUDY = std_centroid(STUDY,ALLEEG, cls , 'dipole');
    end
    comp_to_disp = [subject  ' / ' 'IC' num2str(comp) ];
    cluster_dip_models.posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
    cluster_dip_models.momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
    cluster_dip_models.rv = ALLEEG(abset).dipfit.model(comp).rv;
    cluster_dip_models(2) = STUDY.cluster(cls).dipole;
    if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
        if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
            load('-mat', ALLEEG(abset).dipfit.hdmfile);
            max_r = max(vol.r);
        else
            max_r = max(ALLEEG(abset).dipfit.vol.r);
        end
        dipplot(cluster_dip_models, 'sphere', max_r, 'mri', ALLEEG(abset).dipfit.mrifile,'coordformat', ALLEEG(abset).dipfit.coordformat , ...
           'normlen' ,'on', 'pointout' ,'on','color', {'b', 'r'}, 'dipnames', {comp_to_disp [ STUDY.cluster(cls).name ' mean' ] },...
            'spheres', 'on', 'verbose', 'off', varargin{:});
    else
       dipplot(cluster_dip_models, 'meshdata', ALLEEG(abset).dipfit.hdmfile, 'mri', ALLEEG(abset).dipfit.mrifile,'coordformat', ALLEEG(abset).dipfit.coordformat , ...
          'normlen' ,'on', 'pointout' ,'on','color', {'b', 'r'}, 'dipnames', {comp_to_disp [STUDY.cluster(cls).name ' mean']}, ...
          'spheres', 'on', 'verbose', 'off', varargin{:});
    end
    fig_h = gcf;
    set(fig_h,'Name', [subject ' / ' 'IC' num2str(comp) ', ' STUDY.cluster(cls).name],'NumberTitle','off');
end
        
% -----------------------
% load all dipoles and
% compute dipole centroid
% DEVELOPMENT: this function
% should be the only one to
% access dipole information
% -----------------------
function STUDY = std_centroid(STUDY,ALLEEG, clsind, tmp);

    for clust = 1:length(clsind)
        max_r = 0;
        len = length(STUDY.cluster(clsind(clust)).comps);
        tmppos = [ 0 0 0 ];
        tmpmom = [ 0 0 0 ];
        tmprv = 0;
        ndip = 0;
        for k = 1:len 
            fprintf('.');
            comp  = STUDY.cluster(clsind(clust)).comps(k);
            abset = STUDY.cluster(clsind(clust)).sets(1,k);
            if ~isfield(ALLEEG(abset), 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(abset) ], 'Aborting compute centroid dipole');
               return;
            end
            if ~isempty(ALLEEG(abset).dipfit.model(comp).posxyz)
                ndip   = ndip +1;
                posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
                momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
                if size(posxyz,1) == 2
                    if all(posxyz(2,:) == [ 0 0 0 ])
                        posxyz(2,:) = [];
                        momxyz(2,:) = [];
                    end
                end
                tmppos = tmppos + mean(posxyz,1);
                tmpmom = tmpmom + mean(momxyz,1);
                tmprv = tmprv + ALLEEG(abset).dipfit.model(comp).rv;
                if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                   if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                       load('-mat', ALLEEG(abset).dipfit.hdmfile);
                       max_r = max(max_r, max(vol.r));
                   else % old version of dipfit
                       max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
                   end
               end
            end
        end
        centroid{clust}.dipole.posxyz =  tmppos/ndip;
        centroid{clust}.dipole.momxyz =  tmpmom/ndip;
        centroid{clust}.dipole.rv =  tmprv/ndip;
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical') && (~isfield(ALLEEG(abset).dipfit, 'hdmfile')) %old dipfit
            centroid{clust}.dipole.maxr = max_r;
        end
        STUDY.cluster(clsind(clust)).dipole = centroid{clust}.dipole;
    end
    fprintf('\n');

% --------------------------------
% new function to compute centroid
% was programmed to debug the function
% above but is now used in the code
% --------------------------------
function dipole = computecentroid(alldipoles)

        max_r = 0;
        len = length(alldipoles);
        dipole.posxyz = [ 0 0 0 ];
        dipole.momxyz = [ 0 0 0 ];
        dipole.rv = 0;
        ndip = 0;
        count = 0;
        warningon = 1;
        for k = 1:len 
            if size(alldipoles(k).posxyz,1) == 2
                if all(alldipoles(k).posxyz(2,:) == [ 0 0 0 ])
                    alldipoles(k).posxyz(2,:) = [];
                    alldipoles(k).momxyz(2,:) = [];
                end
            end
            if ~isempty(alldipoles(k).posxyz)
                dipole.posxyz = dipole.posxyz + mean(alldipoles(k).posxyz,1);
                dipole.momxyz = dipole.momxyz + mean(alldipoles(k).momxyz,1);
                dipole.rv     = dipole.rv     + alldipoles(k).rv;
                count = count+1;
            elseif warningon
                disp('Some components do not have dipole information');
                warningon = 0;
            end
        end
        dipole.posxyz = dipole.posxyz/count;
        dipole.momxyz = dipole.momxyz/count;
        dipole.rv     = dipole.rv/count;
        if isfield(alldipoles, 'maxr')
            dipole.maxr = alldipoles(1).max_r;
        end
        
function [cluster_dip_models, options] = dipgroups(ALLEEG, STUDY, cls, comp_to_disp, cluster_dip_models, options);

    % first, extract the subject number
    for n = 1:length(comp_to_disp)
        subjectnum(n,1) = str2num(comp_to_disp{n}(1:3));
    end

    % second, extract group info
    for n = 1:length(subjectnum)
        subj_group{n,1} = ALLEEG(1,subjectnum(n)).group;
    end

    % third, replace the group names with numbers
    for n = 1:length(subj_group)
        for m = 1:length(STUDY.group)
            if strcmp(subj_group{n,1}, STUDY.group{1,m})
                subj_groupnum(n,1) = m;
                break
            end
        end
    end

    % fourth, compute centroid for each group
    for n = 1:length(STUDY.group)
        samegroupIC = find(subj_groupnum==n);
        cluster_dip_models(1,length(subj_groupnum)+n) = computecentroid(cluster_dip_models(1, samegroupIC));
    end

    % fifth, use subj_groupnum as a type of dipole color

        %%%%%%%%%%%%%%%%%%%%% color list %%%%%%%%%%%%%%%%%%%%%
        % This color list was developped for std_envtopo
        % 16 colors names officially supported by W3C specification for HTML
        colors{1,1}  = [1 1 1];            % White
        colors{2,1}  = [1 1 0];            % Yellow
        colors{3,1}  = [1 0 1];            % Fuchsia
        colors{4,1}  = [1 0 0];            % Red
        colors{5,1}  = [0.75  0.75  0.75]; % Silver
        colors{6,1}  = [0.5 0.5 0.5];      % Gray
        colors{7,1}  = [0.5 0.5 0];        % Olive
        colors{8,1}  = [0.5 0 0.5];        % Purple
        colors{9,1}  = [0.5 0 0];          % Maroon
        colors{10,1} = [0 1 1];            % Aqua
        colors{11,1} = [0 1 0];            % Lime
        colors{12,1} = [0 0.5 0.5];        % Teal
        colors{13,1} = [0 0.5 0];          % Green
        colors{14,1} = [0 0 1];            % Blue
        colors{15,1} = [0 0 0.5];          % Navy
        colors{16,1} = [0 0 0];            % Black
        % Silver is twice brighter because silver is used for a background color
        colors{5,1} = [0.875 0.875 0.875];
        % Choosing and sorting 12 colors for line plot, namely Red, Blue, Green, Fuchsia, Lime, Aqua, Maroon, Olive, Purple, Teal, Navy, and Gray
        selectedcolors = colors([4 13 14 3 11 10 9 7 8 12 15 6]);

    % determine the new dip colors
    for n = 1:length(subj_groupnum)
        dip_color{1,n}=selectedcolors{subj_groupnum(n,1)+1};
    end
    for n = 1:length(STUDY.group)
        dip_color{1,end+1}= selectedcolors{n+1};
    end
    
    for n = 1:length(options)
        if      strcmp(options{1,n}, 'color')
            options{1,n+1} = dip_color;
        elseif  strcmp(options{1,n}, 'dipnames')
            dipnames = options{1,n+1};
            for m = 1:length(STUDY.group)
                dipnames{1,length(subj_groupnum)+m}= [STUDY.group{1,m} ' mean'];
            end
            options{1,n+1} = dipnames;
        end
    end
    
% ========================================================================================
% multicolor mode support functions
function str = componentText(~,obj)
% look up component name--used in data cursor callback in 'multicolor' mode
try
  str = obj.Target.UserData.name;
  h=findall(gcf,'type','hggroup');
  h(1).FontSize = 14;
catch
  str = '';
end



