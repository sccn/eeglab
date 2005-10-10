% cls_centroid() - compute clusters centroid in EEGLAB dataset STUDY.
%                          Compute and store the averege (centroid) for
%                          specific clusters or for all clusters in STUDY.
%                          You can compute one or some of the 6 optional
%                          avereges: Scalp map, ERP, Spectra, ERSP , ITC, dipole. 
%
% Usage:    
%   >> [STUDY, centroid] = cls_centroid(STUDY, ALLEEG, clsind, ctr1, ctr2, ...);
%
% Inputs:
%   STUDY        - an eeglab STUDY set (that contains some of these EEGs structures)
%   ALLEEG       - ALLEEG data structure, can also be one EEG dataset structure.
%   clsind         - [vector], cluster indices. Compute centroids only for
%                        indicated clusters. 
%                        If left empty compute centroids for all clusters in STUDY.
%   ctrX            - ['erp'|'spec'|'scalp'|'dipole'|'itc'|'ersp'].   
%                      * The field/s to calculate the mean of each cluster.
%                      * 'erp'   -  the mean ERP of each cluster.
%                      * 'dipole'  - the mean dipole of each cluster.
%                      * 'spec' - the mean spectra of each cluster (with the baseline removed).
%                      * 'scalp' - the mean topoplot scalp map of each cluster.
%                      * 'ersp'  - the mean ERSP of each cluster. 
%                      * 'itc'  - the mean ITC of each cluster. 
%                      If left empty re-compute the centroid for whatever mean centroid was already in the cluster.  
%
% Outputs:
%   STUDY        - the STUDY structure updated with centroids. If centroids
%                        already exist, will overwite them. 
%   centroid     - a cell array of centroid structures, each cell corraspond to a different cluster requested in clsind.
%
% Examples:
%   >> [STUDY, centroid] = cls_centroid(STUDY, ALLEEG,[], 'scalp'); 
%   % For each of the clusters in STUDY an average scalp map is computed. 
%   % The centroids are saved in the STUDY structure under
%   % STUDY.cluster(k).centroid.scalp. The centroids are also returned in a
%   % cell array the size of the clusters (stored in: centroid(k).scalp).
%
%   >> [STUDY, centroid] = cls_centroid(STUDY, ALLEEG,5,'spec','scalp'); 
%   % Same as before but now compute the centroid of cluster 5 only. The centroid has 2 fields one for 
%   % spectrum the other for scalp map.
%
% Authors: Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, Feb 03, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Feb 03, 2005, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

function [STUDY, centroid] = cls_centroid(STUDY,ALLEEG, clsind, varargin);

 if nargin < 3
     help cls_centroid;
     return
 end
 
 if isempty(clsind)
     for k = 2: length(STUDY.cluster) %don't include the ParentCluster
         if ~strncmpi('Notclust',STUDY.cluster(k).name,8) % don't include 'Notclust' clusters
             clsind = [clsind k];
         end
     end
 end
 %defult values
erpC =0;
specC =0 ;
scalpC = 0;
dipoleC = 0;
itcC = 0;
erspC = 0;
 
commands = {};
if isempty(varargin)
    if isfield(STUDY.cluster(clsind(1)).centroid,'scalp')
        commands{end+1} = 'scalp';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'spec')
        commands{end+1} = 'spec';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'erp')
        commands{end+1} = 'erp';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'ersp')
        commands{end+1} = 'ersp';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'itc')
        commands{end+1} = 'itc';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'dipole')
        commands{end+1} = 'dipole';
    end
else
    commands = varargin;
end

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end 
centroid = cell(length(clsind),1);
fprintf('Computing ');
for k = 1:length(clsind)
    for l = 1:Ncond 
        for ind = 1:length(commands)
            ctr = commands{ind};
            switch ctr
                case 'scalp'
                    centroid{k}.scalp = 0; 
                    scalpC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('scalp ');
                    end
                case 'erp'
                    centroid{k}.erp{l} = 0; 
                    erpC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('erp ');
                    end
                case 'spec'
                    centroid{k}.spec{l} = 0; 
                    specC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('spectrum ');
                    end
                case 'ersp'
                    centroid{k}.ersp{l} = 0; 
                    centroid{k}.ersp_limits{l} = 0;
                    erspC =1;
                    if (l ==1) & (k ==1)
                        fprintf('ersp ');
                    end
                case 'itc'
                    centroid{k}.itc{l} = 0; 
                    centroid{k}.itc_limits{l} = 0;
                    itcC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('itc ');
                    end                    
                case 'dipole'
                    dipoleC =1;
                    if (l ==1) & (k ==1)
                        fprintf('dipole ');
                    end
            end
        end
    end
end   
fprintf('centroid \n');
if itcC | erspC | specC | erpC | scalpC
    for cond = 1:Ncond %compute for all conditions
        for clust = 1:length(clsind) %go over all requested clusters
            for ind = 1:size(STUDY.setind,2) %go through all datasets in STUDY
                abset = STUDY.datasetinfo(STUDY.setind(cond,ind)).index;
                compind = find(STUDY.cluster(clsind(clust)).sets(1,:) == ind);
                for k = 1:length(compind) % go through all components
                    comp = STUDY.cluster(clsind(clust)).comps(compind(k));
                    if scalpC & cond == 1  %scalp centroid, does not depend on condition 
                        grid = cls_readscalp(ALLEEG, abset, comp);
                        if isempty(grid)
                            return;
                        end
                        centroid{clust}.scalp = centroid{clust}.scalp + grid;
                    end
                    if erpC %erp centroid
                        [erp, t] = cls_readerp(ALLEEG, abset, comp);
                        fprintf('.');
                        if isempty(erp)
                            return;
                        end
                        if (k==1) & (ind == STUDY.cluster(clsind(clust)).sets(1,1))
                            all_erp = zeros(length(erp),length(STUDY.cluster(clsind(clust)).comps));
                        end
                        all_erp(:,compind(k)) = erp;
                        if (k == length(compind) ) &  (ind == STUDY.cluster(clsind(clust)).sets(1,end)) 
                            [all_erp pol] = comppol(all_erp);
                            centroid{clust}.erp{cond} = mean(all_erp,2);
                            centroid{clust}.erp_t = t;
                        end
                    end
                    if specC %spec centroid
                        [spec, f] = cls_readspec(ALLEEG, abset, comp);
                        fprintf('.');
                        if isempty(spec)
                            return;
                        end
                        centroid{clust}.spec{cond} = centroid{clust}.spec{cond} + spec;
                        centroid{clust}.spec_f = f;
                    end
                    if erspC %ersp centroid
                        fprintf('.');
                        if cond == 1
                            abset = [STUDY.datasetinfo(STUDY.setind(:,ind)).index];
                            [ersp, logfreqs] = cls_readersp(ALLEEG, abset, comp);
                            if isempty(ersp)
                                return;
                            end
                            for m = 1:Ncond
                                centroid{clust}.ersp{m} = centroid{clust}.ersp{m} + ersp(:,:,m);
                                centroid{clust}.ersp_limits{m} = max(floor(max(max(abs(ersp(:,:,m))))), centroid{clust}.ersp_limits{m});
                            end
                            centroid{clust}.ersp_logf = logfreqs;
                            abset = STUDY.datasetinfo(STUDY.setind(cond,ind)).index; %return  default value
                        end
                    end
                    if itcC %itc centroid
                        fprintf('.');
                        [itc, logfreqs] = cls_readitc(ALLEEG, abset, comp);
                        if isempty(itc)
                            return;
                        end
                        centroid{clust}.itc{cond} = centroid{clust}.itc{cond} + itc;
                        centroid{clust}.itc_limits{cond} = max(floor(max(max(abs(itc)))), centroid{clust}.itc_limits{cond}); %ersp image limits 
                        centroid{clust}.itc_logf = logfreqs;
                    end
                end
            end
            fprintf('\n');
        end
	end
end

if dipoleC %dipole centroid
    for clust = 1:length(clsind)
        max_r = 0;
        len = length(STUDY.cluster(clsind(clust)).comps);
        tmppos = 0;
        tmpmom = 0;
        tmprv = 0;
        ndip = 0;
        for k = 1:len 
            fprintf('.');
            abset = STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(clsind(clust)).sets(1,k))).index;
            comp = STUDY.cluster(clsind(clust)).comps(k);
            if ~isfield(ALLEEG(abset), 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(abset) ], 'Aborting compute centroid dipole');
               return;
            end
            if ~isempty(ALLEEG(abset).dipfit.model(comp).posxyz)
                ndip = ndip +1;
                tmppos = tmppos + ALLEEG(abset).dipfit.model(comp).posxyz;
                tmpmom = tmpmom + ALLEEG(abset).dipfit.model(comp).momxyz;
                tmprv = tmprv + ALLEEG(abset).dipfit.model(comp).rv;
                if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                   if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                       eval(['load ' ALLEEG(abset).dipfit.hdmfile]);
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
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical') & (~isfield(ALLEEG(abset).dipfit, 'hdmfile')) %old dipfit
            centroid{clust}.dipole.maxr = max_r;
        end
        STUDY.cluster(clsind(clust)).centroid.dipole = centroid{clust}.dipole;
   end
end

%updat STUDY
for clust =  1:length(clsind) %go over all requested clusters
    for cond  = 1:Ncond
        ncomp = length(STUDY.cluster(clsind(clust)).comps);
        if scalpC & cond == 1%scalp centroid
            centroid{clust}.scalp  = centroid{clust}.scalp/ncomp;
            STUDY.cluster(clsind(clust)).centroid.scalp = centroid{clust}.scalp ;
        end
        if erpC
            STUDY.cluster(clsind(clust)).centroid.erp{cond} = centroid{clust}.erp{cond};
		    STUDY.cluster(clsind(clust)).centroid.erp_t = centroid{clust}.erp_t;
        end
		if specC
            centroid{clust}.spec{cond} = centroid{clust}.spec{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.spec{cond} = centroid{clust}.spec{cond};
		    STUDY.cluster(clsind(clust)).centroid.spec_f = centroid{clust}.spec_f;
        end
        if erspC %ersp centroid
            centroid{clust}.ersp{cond} = centroid{clust}.ersp{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.ersp{cond} = centroid{clust}.ersp{cond};
            STUDY.cluster(clsind(clust)).centroid.ersp_limits{cond} = floor(0.75*centroid{clust}.ersp_limits{cond}); 
            %[round(0.9*min(cell2mat({centroid{clust}.ersp_limits{cond,:}})))  round(0.9*max(cell2mat({centroid{clust}.ersp_limits{cond,:}})))];
            STUDY.cluster(clsind(clust)).centroid.ersp_logf = centroid{clust}.ersp_logf;
        end
        if itcC
            centroid{clust}.itc{cond} = centroid{clust}.itc{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.itc{cond} = centroid{clust}.itc{cond} ;
            STUDY.cluster(clsind(clust)).centroid.itc_limits{cond} = floor(0.75*centroid{clust}.itc_limits{cond});%round(0.9*max(cell2mat({centroid{clust}.itc_limits{cond,:}})));
            STUDY.cluster(clsind(clust)).centroid.itc_logf = centroid{clust}.itc_logf;
        end
        
    end
end
fprintf('\n');