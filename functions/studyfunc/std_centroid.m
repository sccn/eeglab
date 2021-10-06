% std_centroid() - compute cluster centroid in EEGLAB dataset STUDY.
%                  Compute and store the centroid(s) (i.e., mean(s)) 
%                  for some combination of six measures on specified
%                  clusters in a STUDY. Possible measures include: scalp
%                  maps, ERPs, spectra, ERSPs, ITCs, dipole_locations
% Usage:    
%        >> [STUDY, centroid] = std_centroid(STUDY, ALLEEG, ...
%                                              clusters, measure1, measure2, ...);
%
% Inputs:
%   STUDY        - STUDY set 
%   ALLEEG       - ALLEEG dataset vector (else an EEG dataset) containing the STUDY
%                  datasets, typically created using load_ALLEEG().
%   clusters     - [vector] of cluster indices. Computes measure means for the 
%                  specified clusters. {deffault|[]: compute means for all 
%                  STUDY clusters} 
%   measure(s)   - ['erp'|'spec'|'scalp'|'dipole'|'itc'|'ersp'].   
%                  The measures(s) for which to calculate the cluster centroid(s):
%                     'erp'    ->  mean ERP of each cluster.
%                     'dipole' ->  mean dipole of each cluster.
%                     'spec'   ->  mean spectrum of each cluster (baseline removed).
%                     'scalp'  ->  mean topoplot scalp map of each cluster.
%                     'ersp'   ->  mean ERSP of each cluster. 
%                     'itc'    ->  mean ITC of each cluster. 
%                  If [], re-compute the centroid for whichever centroids 
%                  have previously been computed.
% Outputs:
%   STUDY        - input STUDY structure with computed centroids added. 
%                  If the requested centroids already exist, overwrites them. 
%   centroid     - cell array of centroid structures, each cell corrasponding 
%                  to a different cluster requested in 'clusters' (above).
%                  fields of 'centroid' may include centroid.erp, centroid.dipole,
%                  etc. (as above). The structure is similar as the output
%                  of the std_readdata() function (with some fields
%                  about the cluster name and index missing).
% Examples:
%
%   >> [STUDY, centroid] = std_centroid(STUDY, ALLEEG,[], 'scalp'); 
%   % For each of the clusters in STUDY, compute a mean scalp map.
%   % The centroids are saved in the STUDY structure as entries in array
%   % STUDY.cluster(k).centroid.scalp. The centroids are also returned in 
%   % a cell array the size of the clusters (i.e., in: centroid(k).scalp).
%
%   >> [STUDY, centroid] = std_centroid(STUDY, ALLEEG,5,'spec','scalp'); 
%   % Same as above, but now compute only two centroids for Cluster 5. 
%   % The returned 'centroid' has two fields: centroid.scalp and centroid.spec
%
% Authors: Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, Feb 03, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Feb 03, 2005, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

function [STUDY, centroid] = std_centroid(STUDY,ALLEEG, clsind, varargin);

 if nargin < 3
     help std_centroid;
     return
 end
 
 if isempty(clsind)
     for k = 2: length(STUDY.cluster) %don't include the ParentCluster
         if ~strncmpi('Notclust',STUDY.cluster(k).name,8) 
             % don't include 'Notclust' clusters
             clsind = [clsind k];
         end
     end
 end
 %default values
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
                    if (l ==1) && (k ==1)
                        fprintf('scalp ');
                    end
                case 'erp'
                    centroid{k}.erp{l} = 0; 
                    erpC = 1;
                    if (l ==1) && (k ==1)
                        fprintf('erp ');
                    end
                case 'spec'
                    centroid{k}.spec{l} = 0; 
                    specC = 1;
                    if (l ==1) && (k ==1)
                        fprintf('spectrum ');
                    end
                case 'ersp'
                    centroid{k}.ersp{l} = 0; 
                    centroid{k}.ersp_limits{l} = 0;
                    erspC =1;
                    if (l ==1) && (k ==1)
                        fprintf('ersp ');
                    end
                case 'itc'
                    centroid{k}.itc{l} = 0; 
                    centroid{k}.itc_limits{l} = 0;
                    itcC = 1;
                    if (l ==1) && (k ==1)
                        fprintf('itc ');
                    end                    
                case 'dipole'
                    dipoleC =1;
                    if (l ==1) && (k ==1)
                        fprintf('dipole ');
                    end
            end
        end
    end
end   
fprintf('centroid (only done once)\n');
if itcC || erspC || specC || erpC || scalpC
    for clust = 1:length(clsind) %go over all requested clusters
        for cond = 1:Ncond %compute for all conditions
            for k = 1:length(STUDY.cluster(clsind(clust)).comps) % go through all components
                comp  = STUDY.cluster(clsind(clust)).comps(k);
                abset = STUDY.cluster(clsind(clust)).sets(cond,k);
                if scalpC && cond == 1  %scalp centroid, does not depend on condition 
                    grid = std_readtopo(ALLEEG, abset, comp);
                    if isempty(grid)
                        return;
                    end
                    centroid{clust}.scalp = centroid{clust}.scalp + grid;
                end
                if erpC %erp centroid
                    [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                    fprintf('.');
                    if isempty(erp)
                        return;
                    end
                    if (cond==1) && (k==1)
                        all_erp = zeros(length(erp),length(STUDY.cluster(clsind(clust)).comps));
                    end
                    all_erp(:,k) = erp';
                    if k == length(STUDY.cluster(clsind(clust)).comps)
                        [all_erp pol] = std_comppol(all_erp);
                        centroid{clust}.erp{cond} = mean(all_erp,2);
                        centroid{clust}.erp_times = t;
                    end
                end
                if specC %spec centroid
                    [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                    fprintf('.');
                    if isempty(spec)
                        return;
                    end
                    centroid{clust}.spec{cond} = centroid{clust}.spec{cond} + spec;
                    centroid{clust}.spec_freqs = f;
                end
                if erspC %ersp centroid
                    fprintf('.');
                    if cond == 1
                        tmpabset = STUDY.cluster(clsind(clust)).sets(:,k);
                        [ersp, logfreqs, timevals] = std_readersp(ALLEEG, tmpabset, comp, STUDY.preclust.erspclusttimes, ...
                                                                STUDY.preclust.erspclustfreqs );
                        if isempty(ersp)
                            return;
                        end
                        for m = 1:Ncond
                            centroid{clust}.ersp{m} = centroid{clust}.ersp{m} + ersp(:,:,m);
                            centroid{clust}.ersp_limits{m} = max(floor(max(max(abs(ersp(:,:,m))))), centroid{clust}.ersp_limits{m});
                        end
                        centroid{clust}.ersp_freqs  = logfreqs;
                        centroid{clust}.ersp_times = timevals;
                    end
                end
                if itcC %itc centroid
                    fprintf('.');
                    [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, STUDY.preclust.erspclusttimes, ...
                                                                STUDY.preclust.erspclustfreqs );
                    if isempty(itc)
                        return;
                    end
                    centroid{clust}.itc{cond} = centroid{clust}.itc{cond} + itc;
                    centroid{clust}.itc_limits{cond} = max(floor(max(max(abs(itc)))), centroid{clust}.itc_limits{cond}); %ersp image limits 
                    centroid{clust}.itc_freqs  = logfreqs;
                    centroid{clust}.itc_times = timevals;
                end
            end
        end
        if ~scalpC
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
            comp  = STUDY.cluster(clsind(clust)).comps(k);
            abset = STUDY.cluster(clsind(clust)).sets(1,k);
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
        STUDY.cluster(clsind(clust)).centroid.dipole = centroid{clust}.dipole;
   end
end

%update STUDY
for clust =  1:length(clsind) %go over all requested clusters
    for cond  = 1:Ncond
        ncomp = length(STUDY.cluster(clsind(clust)).comps);
        if scalpC && cond == 1%scalp centroid
            centroid{clust}.scalp  = centroid{clust}.scalp/ncomp;
            STUDY.cluster(clsind(clust)).centroid.scalp = centroid{clust}.scalp ;
        end
        if erpC
            STUDY.cluster(clsind(clust)).centroid.erp{cond} = centroid{clust}.erp{cond};
		    STUDY.cluster(clsind(clust)).centroid.erp_times = centroid{clust}.erp_times;
        end
		if specC
            centroid{clust}.spec{cond} = centroid{clust}.spec{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.spec{cond} = centroid{clust}.spec{cond};
		    STUDY.cluster(clsind(clust)).centroid.spec_freqs = centroid{clust}.spec_freqs;
        end
        if erspC %ersp centroid
            centroid{clust}.ersp{cond} = centroid{clust}.ersp{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.ersp{cond} = centroid{clust}.ersp{cond};
            STUDY.cluster(clsind(clust)).centroid.ersp_limits{cond} = floor(0.75*centroid{clust}.ersp_limits{cond}); 
            %[round(0.9*min(cell2mat({centroid{clust}.ersp_limits{cond,:}})))  round(0.9*max(cell2mat({centroid{clust}.ersp_limits{cond,:}})))];
            STUDY.cluster(clsind(clust)).centroid.ersp_freqs = centroid{clust}.ersp_freqs;
            STUDY.cluster(clsind(clust)).centroid.ersp_times = centroid{clust}.ersp_times;
        end
        if itcC
            centroid{clust}.itc{cond} = centroid{clust}.itc{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.itc{cond} = centroid{clust}.itc{cond} ;
            STUDY.cluster(clsind(clust)).centroid.itc_limits{cond} = floor(0.75*centroid{clust}.itc_limits{cond});%round(0.9*max(cell2mat({centroid{clust}.itc_limits{cond,:}})));
            STUDY.cluster(clsind(clust)).centroid.itc_freqs = centroid{clust}.itc_freqs;
            STUDY.cluster(clsind(clust)).centroid.itc_times = centroid{clust}.itc_times;
        end
        
    end
end
fprintf('\n');
